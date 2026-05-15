#!/usr/bin/env python3
"""
Multi-beam RFI cleaning via cross-beam covariance eigendecomposition.

Takes N SIGPROC filterbank files (one per beam), produces N cleaned ones.

Algorithm (per (time-window, freq-window) tile, looped over all tiles):
  1. Whiten each beam-channel by its running robust scale.
  2. Form the BxB cross-beam covariance C.
  3. Eigendecompose C.
  4. Flag eigenpair (lambda_k, u_k) as RFI iff
        lambda_k > MP_edge(B, N) * safety        (matched-filter null)
        AND participation_ratio(u_k) >= min_support  (>= N_min beams)
  5. Subtract that subspace:   x_clean = x - U_rfi U_rfi^T x
  6. Re-quantize to the input bit-depth and write.

Reads everything (nchans, nbits, fch1, foff, tsamp, tstart, ...) from
the SIGPROC headers; only requires that all input files agree on the
header fields that define the data layout (nchans, nbits, nifs,
fch1, foff, tsamp). Tstart is required to be within a configurable
tolerance.

Dependencies: numpy only (plus Python >= 3.8 stdlib).

Memory model (used at peak):
  - one float32 (B, nchans, chunk_samps) buffer for the working chunk
  - one float32 (B, nchans) scratch each for running mean and sigma
  - one float32 (n_freq_batch_tiles * n_time_tiles, B, Fw * Tw) scratch
    + one same-shape scratch for the projection delta
  - per-beam transient buffers during unpack and pack (a few hundred MB
    each, freed immediately)

So peak RAM is ~(B * nchans * chunk_samps * 4 + ~hundreds of MB).
For 72 beams x 2048 chans x 50us tsamp at --chunk-sec 1.0 that's
~12 GB. Drop --chunk-sec to fit smaller nodes.

Usage:
    multibeam_clean.py beam0.fil beam1.fil ... \\
        --outdir cleaned/ \\
        --cores 20 \\
        --min-support 20 \\
        --chunk-sec 1.0 \\
        --tile-time-ms 100 \\
        --tile-freq-chans 4

All tunable parameters listed under --help.
"""

# --------------------------------------------------------------------------
# IMPORTANT: parse --cores BEFORE importing numpy, so that BLAS thread pools
# pick up the right env vars.
# --------------------------------------------------------------------------
import os
import sys


def _early_parse_cores(argv):
    for i, a in enumerate(argv):
        if a in ("--cores", "-j") and i + 1 < len(argv):
            try:
                return int(argv[i + 1])
            except ValueError:
                pass
        if a.startswith("--cores="):
            try:
                return int(a.split("=", 1)[1])
            except ValueError:
                pass
    return None


_cores_env = _early_parse_cores(sys.argv)
if _cores_env is not None and _cores_env > 0:
    os.environ["OMP_NUM_THREADS"]      = str(_cores_env)
    os.environ["OPENBLAS_NUM_THREADS"] = str(_cores_env)
    os.environ["MKL_NUM_THREADS"]      = str(_cores_env)
    os.environ["NUMEXPR_NUM_THREADS"]  = str(_cores_env)
    os.environ["BLIS_NUM_THREADS"]     = str(_cores_env)

import argparse
import json
import math
import struct
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# threadpoolctl is a thin pure-Python lib that lets us cap BLAS thread pools
# *inside a context manager*. Conda/pip numpy stacks ship it as a transitive
# dep almost universally. Soft import: if missing, we just don't parallelize
# BLAS-heavy work (parallel workers would oversubscribe).
try:
    from threadpoolctl import threadpool_limits as _threadpool_limits
    _HAS_TPCTL = True
except Exception:
    _threadpool_limits = None
    _HAS_TPCTL = False


# ==========================================================================
# SIGPROC header I/O (preserves field order so the writer round-trips)
# ==========================================================================

_TAG_INT    = {"machine_id", "telescope_id", "data_type", "nchans", "nbits",
               "nifs", "nbeams", "ibeam", "barycentric", "pulsarcentric",
               "nbins", "nsamples"}
_TAG_DOUBLE = {"tstart", "tsamp", "fch1", "foff", "refdm",
               "src_raj", "src_dej", "az_start", "za_start",
               "period", "fchannel"}
_TAG_INT64  = {"npuls"}
_TAG_BYTE   = {"signed"}
_TAG_STRING = {"rawdatafile", "source_name"}


def _read_str(f) -> Optional[str]:
    head = f.read(4)
    if len(head) < 4:
        return None
    (n,) = struct.unpack("<i", head)
    if not (1 <= n <= 80):
        raise ValueError(f"sigproc string length out of range: {n}")
    s = f.read(n)
    if len(s) < n:
        raise IOError("truncated sigproc string")
    return s.decode("ascii", errors="strict")


def _write_str(f, s: str) -> None:
    b = s.encode("ascii")
    if not (1 <= len(b) <= 80):
        raise ValueError(f"sigproc string too long/empty: {s!r}")
    f.write(struct.pack("<i", len(b)))
    f.write(b)


@dataclass
class SigprocHeader:
    """Ordered list of (tag, type, value) tuples + byte length."""
    fields: List[Tuple[str, str, object]] = field(default_factory=list)
    header_bytes: int = 0

    def get(self, tag: str, default=None):
        for t, _, v in self.fields:
            if t == tag:
                return v
        return default

    def set(self, tag: str, value) -> None:
        for i, (t, ty, _) in enumerate(self.fields):
            if t == tag:
                self.fields[i] = (t, ty, value)
                return
        if tag in _TAG_INT:    self.fields.append((tag, "int",    int(value)))
        elif tag in _TAG_DOUBLE: self.fields.append((tag, "double", float(value)))
        elif tag in _TAG_INT64:  self.fields.append((tag, "int64",  int(value)))
        elif tag in _TAG_BYTE:   self.fields.append((tag, "byte",   int(value)))
        elif tag in _TAG_STRING: self.fields.append((tag, "string", str(value)))
        else: raise ValueError(f"unknown sigproc tag {tag!r}")


def read_sigproc_header(f) -> SigprocHeader:
    tag = _read_str(f)
    if tag != "HEADER_START":
        raise ValueError(f"not a sigproc filterbank: first tag {tag!r}")
    fields: List[Tuple[str, str, object]] = []
    while True:
        tag = _read_str(f)
        if tag is None:
            raise IOError("truncated sigproc header (no HEADER_END)")
        if tag == "HEADER_END":
            break
        if tag in _TAG_INT:
            (v,) = struct.unpack("<i", f.read(4))
            fields.append((tag, "int", v))
        elif tag in _TAG_DOUBLE:
            (v,) = struct.unpack("<d", f.read(8))
            fields.append((tag, "double", v))
        elif tag in _TAG_INT64:
            (v,) = struct.unpack("<q", f.read(8))
            fields.append((tag, "int64", v))
        elif tag in _TAG_BYTE:
            (v,) = struct.unpack("<b", f.read(1))
            fields.append((tag, "byte", v))
        elif tag in _TAG_STRING:
            v = _read_str(f)
            fields.append((tag, "string", v))
        else:
            raise ValueError(f"unknown sigproc tag {tag!r}")
    return SigprocHeader(fields=fields, header_bytes=f.tell())


def write_sigproc_header(f, hdr: SigprocHeader) -> int:
    start = f.tell()
    _write_str(f, "HEADER_START")
    for tag, ty, value in hdr.fields:
        _write_str(f, tag)
        if ty == "int":    f.write(struct.pack("<i", int(value)))
        elif ty == "double": f.write(struct.pack("<d", float(value)))
        elif ty == "int64":  f.write(struct.pack("<q", int(value)))
        elif ty == "byte":   f.write(struct.pack("<b", int(value)))
        elif ty == "string": _write_str(f, str(value))
        else: raise ValueError(f"bad field type {ty!r} for tag {tag!r}")
    _write_str(f, "HEADER_END")
    return f.tell() - start


# ==========================================================================
# Pack / unpack between nbits-packed bytes and float32 samples
# All decoders write into a caller-supplied `out` slice when possible so we
# don't allocate (chunk_samps * nchans * 4) bytes per beam per chunk.
# ==========================================================================

def unpack_bits_to_float32(packed: np.ndarray, nbits: int,
                           nsamps_total: int,
                           out: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Decode `packed` uint8 bytes to `nsamps_total` float32 samples. If `out`
    is provided (length >= nsamps_total, float32), decode in place into it
    and return a slice view; otherwise allocate.
    SIGPROC sample order; for nbits=2, 4 samples per byte, low bits first.
    """
    if out is None:
        out = np.empty(nsamps_total, dtype=np.float32)
    o = out[:nsamps_total]

    if nbits == 8:
        np.copyto(o, packed[:nsamps_total].astype(np.float32, copy=False))
        return o
    if nbits == 16:
        np.copyto(o, packed.view(np.int16)[:nsamps_total].astype(np.float32, copy=False))
        return o
    if nbits == 32:
        np.copyto(o, packed.view(np.float32)[:nsamps_total])
        return o
    if nbits == 4:
        n_pairs = (nsamps_total + 1) // 2
        b = packed[:n_pairs]
        o2 = out[:2 * n_pairs].reshape(n_pairs, 2)
        o2[:, 0] = b & 0x0F
        o2[:, 1] = (b >> 4) & 0x0F
        return out[:nsamps_total]
    if nbits == 2:
        # one byte -> 4 samples. Process in 64 MiB packed slabs to keep
        # transient intermediates small.
        n_full = (nsamps_total + 3) // 4
        b = packed[:n_full]
        slab_bytes = 1 << 26  # 64 MiB packed -> 256 MiB float32
        for s0 in range(0, n_full, slab_bytes):
            s1 = min(s0 + slab_bytes, n_full)
            bb = b[s0:s1]
            o4 = out[4 * s0:4 * s1].reshape(-1, 4)
            o4[:, 0] = bb & 0x3
            o4[:, 1] = (bb >> 2) & 0x3
            o4[:, 2] = (bb >> 4) & 0x3
            o4[:, 3] = (bb >> 6) & 0x3
        return out[:nsamps_total]
    if nbits == 1:
        n_full = (nsamps_total + 7) // 8
        b = packed[:n_full]
        o8 = out[:8 * n_full].reshape(-1, 8)
        for k in range(8):
            o8[:, k] = (b >> k) & 0x1
        return out[:nsamps_total]
    raise ValueError(f"unsupported nbits: {nbits}")


def pack_float32_to_bits(samples: np.ndarray, nbits: int,
                         out: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Re-quantize float32 samples (already mapped to the correct integer range)
    to nbits-packed uint8 bytes. SIGPROC sample order; for nbits<8 the input
    length must be a multiple of (8//nbits).
    """
    n = samples.size
    if nbits == 8:
        return np.clip(samples, 0, 255, out=out.astype(np.float32) if out is not None else None).astype(np.uint8)
    if nbits == 16:
        return np.clip(samples, -32768, 32767).astype(np.int16).view(np.uint8)
    if nbits == 32:
        return samples.astype(np.float32).view(np.uint8)
    if nbits == 4:
        if n % 2: raise ValueError("4-bit pack needs even length")
        q = np.clip(np.rint(samples), 0, 15).astype(np.uint8)
        return (q[0::2] | (q[1::2] << 4)).astype(np.uint8)
    if nbits == 2:
        if n % 4: raise ValueError("2-bit pack needs length % 4 == 0")
        # Process in 64 MiB float slabs to bound the rint/clip transients.
        n_out = n // 4
        if out is None:
            out = np.empty(n_out, dtype=np.uint8)
        slab = 1 << 24  # 16 M samples -> 64 MiB float32
        for s0 in range(0, n, slab):
            s1 = min(s0 + slab, n)
            x = np.rint(samples[s0:s1])
            np.clip(x, 0, 3, out=x)
            q = x.astype(np.uint8)
            o0 = s0 // 4
            o1 = s1 // 4
            np.copyto(out[o0:o1],
                       q[0::4] | (q[1::4] << 2) | (q[2::4] << 4) | (q[3::4] << 6))
        return out
    if nbits == 1:
        if n % 8: raise ValueError("1-bit pack needs length % 8 == 0")
        q = np.clip(np.rint(samples), 0, 1).astype(np.uint8)
        out = np.zeros(q.size // 8, dtype=np.uint8)
        for k in range(8):
            out |= (q[k::8] << k)
        return out
    raise ValueError(f"unsupported nbits: {nbits}")


# ==========================================================================
# Per-beam whitening stats (low-memory, robust)
# ==========================================================================

def _beam_stats_trimmed(x: np.ndarray, clip_sigma: float = 4.0,
                        n_iter: int = 2) -> Tuple[np.ndarray, np.ndarray]:
    """
    Robust per-channel mean/sigma using sigma-clipped iterative estimation.
    x: (C, T) float32. Returns (mean (C,), sigma (C,)) as float32.

    This is much cheaper than np.median + np.median-of-abs: each iter is
    just two reductions + a mask, all SIMD-vectorized in numpy. Memory
    transient is bounded by the boolean mask (~C*T bytes) which we reuse.
    """
    C, T = x.shape
    mu = x.mean(axis=1)
    sigma = x.std(axis=1)
    if T < 16 or n_iter <= 0:
        return mu.astype(np.float32), np.maximum(sigma, 1e-3).astype(np.float32)
    for _ in range(n_iter):
        lo = (mu - clip_sigma * sigma)[:, None]
        hi = (mu + clip_sigma * sigma)[:, None]
        # boolean mask, ~C*T bytes
        mask = (x >= lo) & (x <= hi)
        cnt = mask.sum(axis=1).astype(np.float32)
        # masked sums via where -> 0 outside, no extra allocation beyond mask
        mu_new = (np.where(mask, x, 0.0)).sum(axis=1) / np.maximum(cnt, 1.0)
        dev = np.where(mask, x - mu_new[:, None], 0.0)
        var_new = (dev * dev).sum(axis=1) / np.maximum(cnt - 1.0, 1.0)
        sigma_new = np.sqrt(np.maximum(var_new, 0.0))
        mu, sigma = mu_new, sigma_new
    return mu.astype(np.float32), np.maximum(sigma, 1e-3).astype(np.float32)


def _beam_stats_mad(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Median + 1.4826*MAD per channel. Allocates ~x.size bytes for |x-med|."""
    med = np.median(x, axis=1)
    mad = np.median(np.abs(x - med[:, None]), axis=1)
    sigma = 1.4826 * mad
    bad = sigma <= 0
    if bad.any():
        sigma[bad] = np.maximum(x[bad].std(axis=1), 1e-6)
    return med.astype(np.float32), sigma.astype(np.float32)


def _beam_stats_plain(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Plain per-channel mean + std (axis=1). Cheapest option; not robust to
    outliers, but for cleaning across beams the eigen analysis only needs
    consistency, not absolute correctness."""
    mu = x.mean(axis=1)
    sigma = x.std(axis=1)
    return mu.astype(np.float32), np.maximum(sigma, 1e-3).astype(np.float32)


# ==========================================================================
# Parallel helpers (memory-bound numpy ufuncs are single-threaded by default,
# which is the single biggest perf footgun on large-chunk multi-core nodes).
# ==========================================================================

def _split_workload(n: int, n_workers: int) -> List[Tuple[int, int]]:
    """Return list of (start, stop) index ranges covering [0, n)."""
    if n_workers <= 1 or n <= 1:
        return [(0, n)]
    parts = np.array_split(np.arange(n), min(n_workers, n))
    return [(int(p[0]), int(p[-1]) + 1) for p in parts if len(p)]


def _parallel_whiten(chunk: np.ndarray, mean: np.ndarray, sigma: np.ndarray,
                     pool: ThreadPoolExecutor, n_workers: int) -> None:
    """In place: chunk = (chunk - mean[..., None]) / sigma[..., None]
    parallelized across the leading (beam) axis."""
    ranges = _split_workload(chunk.shape[0], n_workers)
    def _work(rng):
        i0, i1 = rng
        chunk[i0:i1] -= mean[i0:i1, :, None]
        chunk[i0:i1] /= sigma[i0:i1, :, None]
    if len(ranges) == 1:
        _work(ranges[0]); return
    list(pool.map(_work, ranges))


def _parallel_dewhiten(chunk: np.ndarray, mean: np.ndarray, sigma: np.ndarray,
                       pool: ThreadPoolExecutor, n_workers: int) -> None:
    """Inverse of _parallel_whiten."""
    ranges = _split_workload(chunk.shape[0], n_workers)
    def _work(rng):
        i0, i1 = rng
        chunk[i0:i1] *= sigma[i0:i1, :, None]
        chunk[i0:i1] += mean[i0:i1, :, None]
    if len(ranges) == 1:
        _work(ranges[0]); return
    list(pool.map(_work, ranges))


# ==========================================================================
# Eigen-cleaning kernel (in-place, freq-batch streaming)
# ==========================================================================

@dataclass
class CleanParams:
    tile_time_samp: int
    tile_freq_chans: int
    min_support: float
    safety: float
    max_rank: int
    freq_batch_tiles: int = 8     # number of freq tiles to batch per einsum
    prefilter: bool = True        # skip eigh on tiles that can't possibly host RFI


@dataclass
class CleanDiag:
    n_tiles: int = 0
    n_rfi_eigvals: int = 0
    n_tiles_with_rfi: int = 0
    n_tiles_skipped: int = 0     # passed pre-filter, no eigh done
    max_eigval: float = 0.0
    mp_edge: float = 0.0
    safety_threshold: float = 0.0
    t_cov: float = 0.0
    t_eigh: float = 0.0
    t_project: float = 0.0
    t_reshape: float = 0.0


def _clean_one_batch(whitened, c0, c1, nFb, Fw, nTt, Tw, B, N,
                     threshold, mr, min_support, prefilter, tr_skip,
                     timers):
    """Clean a single freq-tile batch in place inside `whitened`.
    Returns (n_rfi_eigvals, n_tiles_with_rfi, n_tiles_skipped, max_eigval)."""
    t0 = time.perf_counter()
    view = whitened[:, c0:c1, :nTt * Tw]
    tiles = np.ascontiguousarray(
        view.reshape(B, nFb, Fw, nTt, Tw)
            .transpose(1, 3, 0, 2, 4)
            .reshape(nFb, nTt, B, N)
    )
    timers["reshape"] += time.perf_counter() - t0

    t0 = time.perf_counter()
    cov = tiles @ tiles.swapaxes(-2, -1)
    cov *= (1.0 / float(N))
    timers["cov"] += time.perf_counter() - t0

    if prefilter:
        tr = np.einsum("FTbb->FT", cov)
        active = tr > tr_skip
    else:
        active = np.ones((nFb, nTt), dtype=bool)

    n_active = int(active.sum())
    n_total = nFb * nTt

    if n_active == 0:
        t0 = time.perf_counter()
        whitened[:, c0:c1, :nTt * Tw] = (
            tiles.reshape(nFb, nTt, B, Fw, Tw)
                 .transpose(2, 0, 3, 1, 4)
                 .reshape(B, nFb * Fw, nTt * Tw)
        )
        timers["reshape"] += time.perf_counter() - t0
        return (0, 0, n_total, 0.0)

    idx_active = np.argwhere(active)
    cov_active = cov[active]
    t0 = time.perf_counter()
    eigvals_a, eigvecs_a = np.linalg.eigh(cov_active)
    timers["eigh"] += time.perf_counter() - t0

    eigvals_a = eigvals_a[..., ::-1]
    eigvecs_a = eigvecs_a[..., ::-1]
    top_eigvals = eigvals_a[..., :mr]
    top_eigvecs = eigvecs_a[..., :, :mr]

    significant = top_eigvals > threshold
    u_sq = top_eigvecs * top_eigvecs
    u4_sum = (u_sq * u_sq).sum(axis=-2)
    pr = 1.0 / np.maximum(u4_sum, 1e-30)
    enough_support = pr >= float(min_support)
    is_rfi = significant & enough_support

    max_eigval = float(top_eigvals[..., 0].max()) if top_eigvals.size else 0.0
    rfi_per_tile = is_rfi.any(axis=-1)
    n_rfi_here = int(is_rfi.sum())
    n_tiles_rfi = int(rfi_per_tile.sum())
    n_skipped = n_total - n_active

    if n_rfi_here:
        sub_mask = rfi_per_tile
        sub_top_eigvecs = top_eigvecs[sub_mask]
        sub_is_rfi = is_rfi[sub_mask]
        sub_idx = idx_active[sub_mask]
        U_rfi = sub_top_eigvecs * sub_is_rfi[..., None, :].astype(np.float32)
        sub_tiles = tiles[sub_idx[:, 0], sub_idx[:, 1]]
        t0 = time.perf_counter()
        v = U_rfi.swapaxes(-2, -1) @ sub_tiles
        sub_tiles -= U_rfi @ v
        timers["project"] += time.perf_counter() - t0
        tiles[sub_idx[:, 0], sub_idx[:, 1]] = sub_tiles

    t0 = time.perf_counter()
    whitened[:, c0:c1, :nTt * Tw] = (
        tiles.reshape(nFb, nTt, B, Fw, Tw)
             .transpose(2, 0, 3, 1, 4)
             .reshape(B, nFb * Fw, nTt * Tw)
    )
    timers["reshape"] += time.perf_counter() - t0
    return (n_rfi_here, n_tiles_rfi, n_skipped, max_eigval)


def clean_chunk(whitened: np.ndarray, params: CleanParams,
                pool: Optional[ThreadPoolExecutor] = None,
                n_workers: int = 1) -> CleanDiag:
    """
    In-place clean of a whitened (B, C, T) chunk.

    If `pool` is provided and n_workers > 1, freq-tile batches are processed
    in parallel; the gather/scatter reshapes (which dominate wall time on
    large chunks) parallelize cleanly because batches touch disjoint slices
    of `whitened`.
    """
    B, C, T = whitened.shape
    Fw = max(1, int(params.tile_freq_chans))
    Tw = max(1, int(params.tile_time_samp))
    nFt = C // Fw
    nTt = T // Tw
    C_use = nFt * Fw
    T_use = nTt * Tw

    diag = CleanDiag()
    if nFt == 0 or nTt == 0:
        return diag

    N = Fw * Tw
    q = float(B) / float(N)
    mp_edge = (1.0 + math.sqrt(q)) ** 2
    threshold = mp_edge * params.safety
    mr = params.max_rank if params.max_rank > 0 else B
    mr = min(mr, B)
    fbatch = max(1, int(params.freq_batch_tiles))
    tr_skip = threshold + (B - 1)

    # Build the work list of (f0, f1) batches.
    batches = []
    for f0 in range(0, nFt, fbatch):
        f1 = min(f0 + fbatch, nFt)
        batches.append((f0, f1))

    timers = {"reshape": 0.0, "cov": 0.0, "eigh": 0.0, "project": 0.0}

    # Parallel path: each worker handles one freq-tile batch with BLAS locked
    # to a single thread so the 20-way worker pool doesn't oversubscribe
    # against MKL's internal 20-thread BLAS pool. Memory-bound reshapes
    # parallelize linearly; small BLAS calls (per-tile project) are also
    # better in many-1-thread mode than few-many-thread mode.
    use_parallel = (pool is not None and n_workers > 1 and len(batches) > 1
                    and _HAS_TPCTL)

    if use_parallel:
        thread_timers: List[Dict[str, float]] = [
            {"reshape": 0.0, "cov": 0.0, "eigh": 0.0, "project": 0.0}
            for _ in batches
        ]
        def _do(args):
            i, (f0, f1) = args
            nFb = f1 - f0
            with _threadpool_limits(limits=1):
                return _clean_one_batch(
                    whitened, f0 * Fw, f1 * Fw,
                    nFb, Fw, nTt, Tw, B, N,
                    threshold, mr, params.min_support,
                    params.prefilter, tr_skip,
                    thread_timers[i],
                )
        results = list(pool.map(_do, list(enumerate(batches))))
        for tt in thread_timers:
            for k in timers: timers[k] += tt[k]
    else:
        results = []
        for (f0, f1) in batches:
            nFb = f1 - f0
            results.append(_clean_one_batch(
                whitened, f0 * Fw, f1 * Fw,
                nFb, Fw, nTt, Tw, B, N,
                threshold, mr, params.min_support,
                params.prefilter, tr_skip,
                timers,
            ))

    total_rfi = sum(r[0] for r in results)
    total_rfi_tiles = sum(r[1] for r in results)
    n_skipped = sum(r[2] for r in results)
    max_eigval = max((r[3] for r in results), default=0.0)
    total_tiles = sum((f1 - f0) for f0, f1 in batches) * nTt

    diag.n_tiles = total_tiles
    diag.n_rfi_eigvals = total_rfi
    diag.n_tiles_with_rfi = total_rfi_tiles
    diag.n_tiles_skipped = n_skipped
    diag.max_eigval = max_eigval
    diag.mp_edge = mp_edge
    diag.safety_threshold = threshold
    diag.t_cov = timers["cov"]
    diag.t_eigh = timers["eigh"]
    diag.t_project = timers["project"]
    diag.t_reshape = timers["reshape"]
    return diag


# ==========================================================================
# Multi-beam I/O + per-beam helpers
# ==========================================================================

@dataclass
class BeamFile:
    path: Path
    header: SigprocHeader
    fp_in: object = None
    fp_out: object = None
    data_bytes_total: int = 0
    nsamples_total: int = 0
    out_path: Optional[Path] = None


def open_input(path: Path) -> BeamFile:
    fp = open(path, "rb")
    hdr = read_sigproc_header(fp)
    size = path.stat().st_size
    data_bytes = size - hdr.header_bytes
    nchans = hdr.get("nchans"); nbits = hdr.get("nbits"); nifs = hdr.get("nifs") or 1
    bps = (nbits * nchans * nifs) / 8.0
    nsamples = int(data_bytes // bps) if bps > 0 else 0
    return BeamFile(path=path, header=hdr, fp_in=fp,
                    data_bytes_total=int(data_bytes), nsamples_total=nsamples)


def validate_compatible(beams: List[BeamFile], tstart_tol_sec: float) -> None:
    if not beams:
        raise ValueError("no input files")
    ref = beams[0].header
    for key in ("nchans", "nbits", "nifs"):
        ref_v = ref.get(key)
        for b in beams[1:]:
            v = b.header.get(key)
            if v != ref_v:
                raise ValueError(
                    f"{b.path.name}: {key}={v} differs from "
                    f"{beams[0].path.name}: {key}={ref_v}")
    for key in ("fch1", "foff", "tsamp"):
        ref_v = ref.get(key)
        for b in beams[1:]:
            v = b.header.get(key)
            if not (ref_v is not None and v is not None and
                    abs(float(v) - float(ref_v)) < 1e-9 * max(1.0, abs(float(ref_v)))):
                raise ValueError(
                    f"{b.path.name}: {key}={v} differs from "
                    f"{beams[0].path.name}: {key}={ref_v}")
    ref_tstart = float(ref.get("tstart", 0.0))
    tol_days = tstart_tol_sec / 86400.0
    for b in beams[1:]:
        ts = float(b.header.get("tstart", 0.0))
        if abs(ts - ref_tstart) > tol_days:
            raise ValueError(
                f"{b.path.name}: tstart {ts} differs from {ref_tstart} by "
                f">{tstart_tol_sec}s tolerance")
    if (ref.get("nifs") or 1) != 1:
        raise ValueError("multibeam_clean: only nifs=1 is supported")


def open_outputs(beams: List[BeamFile], outdir: Path, suffix: str,
                 out_nbits: Optional[int]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for b in beams:
        b.out_path = outdir / f"{b.path.stem}{suffix}.fil"
        b.fp_out = open(b.out_path, "wb")
        hdr_to_write = SigprocHeader(fields=list(b.header.fields))
        if out_nbits is not None and out_nbits != b.header.get("nbits"):
            hdr_to_write.set("nbits", int(out_nbits))
        if hdr_to_write.get("rawdatafile") is not None:
            hdr_to_write.set("rawdatafile", b.out_path.name)
        write_sigproc_header(b.fp_out, hdr_to_write)


# ==========================================================================
# Main pipeline (streaming, in-place where possible)
# ==========================================================================

@dataclass
class PipelineParams:
    chunk_sec: float
    whitening_tau_sec: float
    tile_time_ms: float
    tile_freq_chans: int
    min_support: float
    safety: float
    max_rank: int
    out_nbits: Optional[int]
    suffix: str
    tstart_tol_sec: float
    diag_path: Optional[Path]
    cores: int
    freq_batch_tiles: int
    robust_stats: str          # 'trimmed' (default) or 'mad'
    prefilter: bool = True     # cheap trace-based skip on RFI-free tiles


# ----- threading / environment introspection -----

def _affinity_cpus() -> Optional[int]:
    try:
        return len(os.sched_getaffinity(0))         # Linux: cgroup-aware
    except Exception:
        return None


def _threading_banner(requested_cores: int) -> List[str]:
    lines: List[str] = []
    aff = _affinity_cpus()
    lines.append(f"[env] python={sys.version.split()[0]}  numpy={np.__version__}")
    lines.append(f"[env] requested cores       : {requested_cores}")
    lines.append(f"[env] OMP_NUM_THREADS       : {os.environ.get('OMP_NUM_THREADS', '<unset>')}")
    lines.append(f"[env] MKL_NUM_THREADS       : {os.environ.get('MKL_NUM_THREADS', '<unset>')}")
    lines.append(f"[env] OPENBLAS_NUM_THREADS  : {os.environ.get('OPENBLAS_NUM_THREADS', '<unset>')}")
    lines.append(f"[env] os.cpu_count()        : {os.cpu_count()}")
    lines.append(f"[env] cgroup affinity CPUs  : {aff if aff is not None else 'n/a'}")
    if aff is not None and aff < requested_cores:
        lines.append(f"[env] WARNING: only {aff} CPU(s) actually granted by the kernel "
                     f"(you asked for {requested_cores}). BLAS will be limited.")
    lines.append(f"[env] threadpoolctl available: {_HAS_TPCTL}   "
                 f"(needed for the parallel clean_chunk path)")
    # Try a tiny GEMM and time it to confirm BLAS is actually threaded
    try:
        n = 1024
        A = np.random.randn(n, n).astype(np.float32)
        B = np.random.randn(n, n).astype(np.float32)
        # warmup
        _ = A @ B
        t0 = time.monotonic()
        for _ in range(3): _ = A @ B
        dt = (time.monotonic() - t0) / 3.0
        gflops = (2.0 * n * n * n) / dt / 1e9
        lines.append(f"[env] sgemm({n}x{n}) GFLOPS  : {gflops:7.1f}   "
                     f"({dt*1000:.1f} ms/call;  threaded BLAS should give "
                     f">>{20*requested_cores} GFLOPS)")
    except Exception as e:
        lines.append(f"[env] BLAS micro-benchmark failed: {e}")
    # Try a small eigh to see if LAPACK is multi-threaded
    try:
        n = 256
        S = np.random.randn(n, n).astype(np.float32)
        S = S @ S.T
        # warmup
        _ = np.linalg.eigh(S)
        t0 = time.monotonic()
        for _ in range(3): _ = np.linalg.eigh(S)
        dt = (time.monotonic() - t0) / 3.0
        lines.append(f"[env] eigh({n}x{n})           : {dt*1000:7.1f} ms/call   "
                     f"(if this is the same as 1-core LAPACK, your numpy is "
                     f"using netlib LAPACK)")
    except Exception as e:
        lines.append(f"[env] eigh micro-benchmark failed: {e}")
    return lines


def process(beams: List[BeamFile], params: PipelineParams,
            quiet: bool = False) -> Dict:
    nchans = beams[0].header.get("nchans")
    nbits  = beams[0].header.get("nbits")
    tsamp  = float(beams[0].header.get("tsamp"))
    B = len(beams)
    out_nbits = params.out_nbits if params.out_nbits is not None else nbits

    chunk_samps = int(round(params.chunk_sec / tsamp))
    tile_time_samps = max(1, int(round(params.tile_time_ms * 1e-3 / tsamp)))

    # Snap chunk to a multiple of tile_time_samps and to a byte boundary.
    if chunk_samps % tile_time_samps:
        chunk_samps = (chunk_samps // tile_time_samps) * tile_time_samps
    if chunk_samps == 0:
        raise ValueError("chunk_sec too short for tile_time_ms; pick larger chunk")
    if (chunk_samps * nchans * nbits) % 8:
        raise ValueError(
            "chunk size does not produce an integer number of bytes per beam; "
            "increase --tile-time-ms or --chunk-sec to align")
    bytes_per_chunk = (chunk_samps * nchans * nbits) // 8

    nsamps_avail = min(b.nsamples_total for b in beams)
    n_chunks = nsamps_avail // chunk_samps

    alpha = float(params.chunk_sec) / float(params.whitening_tau_sec)
    alpha = max(min(alpha, 1.0), 1e-3)

    clean_params = CleanParams(
        tile_time_samp = tile_time_samps,
        tile_freq_chans = max(1, params.tile_freq_chans),
        min_support = params.min_support,
        safety = params.safety,
        max_rank = params.max_rank,
        freq_batch_tiles = params.freq_batch_tiles,
        prefilter = params.prefilter,
    )

    if params.robust_stats == "mad":
        stats_fn = _beam_stats_mad
    elif params.robust_stats == "trimmed":
        stats_fn = _beam_stats_trimmed
    else:
        stats_fn = _beam_stats_plain

    if not quiet:
        for line in _threading_banner(params.cores):
            print(line)
        chunk_gb = (B * nchans * chunk_samps * 4) / (1 << 30)
        print(f"[multibeam_clean] beams={B}  nchans={nchans}  nbits={nbits}  "
              f"tsamp={tsamp*1e6:.3f} us")
        print(f"[multibeam_clean] usable samples per beam = {nsamps_avail}  "
              f"({nsamps_avail*tsamp:.2f} s = {nsamps_avail*tsamp/3600:.3f} hr)")
        print(f"[multibeam_clean] chunk = {chunk_samps} samples "
              f"({chunk_samps*tsamp:.4f} s);  {n_chunks} chunks total;  "
              f"chunk buffer = {chunk_gb:.2f} GB float32")
        print(f"[multibeam_clean] tile_time = {tile_time_samps} samples "
              f"({tile_time_samps*tsamp*1e3:.3f} ms);  tile_freq = "
              f"{clean_params.tile_freq_chans} chans;  freq_batch_tiles = "
              f"{params.freq_batch_tiles}")
        print(f"[multibeam_clean] whitening tau = {params.whitening_tau_sec:.3f} s "
              f"(alpha = {alpha:.4f});  robust_stats = {params.robust_stats}")
        print(f"[multibeam_clean] min_support = {params.min_support}  "
              f"safety = {params.safety}  max_rank = {params.max_rank}  "
              f"out_nbits = {out_nbits}")
        print(f"[multibeam_clean] cores = {params.cores}")

    # Pre-allocate the (B, C, T) chunk buffer once, reuse across iterations.
    chunk = np.empty((B, nchans, chunk_samps), dtype=np.float32)
    # Running whitening state.
    state_mean  = np.zeros((B, nchans), dtype=np.float32)
    state_sigma = np.ones ((B, nchans), dtype=np.float32)
    state_init  = False

    diag_fp = open(params.diag_path, "w") if params.diag_path else None
    io_pool = ThreadPoolExecutor(max_workers=min(B, max(2, params.cores)))

    t0 = time.monotonic()
    last_report = t0
    total_rfi = 0
    total_rfi_tiles = 0
    total_tiles = 0
    # Per-stage cumulative timers (reset every progress report)
    t_acc = {"read": 0.0, "unpack": 0.0, "stats": 0.0, "whiten": 0.0,
             "clean": 0.0, "dewhiten": 0.0, "pack": 0.0, "write": 0.0,
             "cov": 0.0, "eigh": 0.0, "project": 0.0, "reshape": 0.0}
    chunks_since_report = 0

    try:
        for ck in range(n_chunks):
            # --- 1. Read packed bytes per beam (parallel) ---
            ts = time.perf_counter()
            def _read(b):
                return b.fp_in.read(bytes_per_chunk)
            packed_list = list(io_pool.map(_read, beams))
            t_acc["read"] += time.perf_counter() - ts
            if any(len(p) != bytes_per_chunk for p in packed_list):
                if not quiet:
                    print(f"[multibeam_clean] short read at chunk {ck}, stopping",
                          file=sys.stderr)
                break

            # --- 2. Unpack into the per-beam slice of `chunk` (parallel).
            ts = time.perf_counter()
            def _unpack(arg):
                idx, packed = arg
                arr = np.frombuffer(packed, dtype=np.uint8)
                flat = unpack_bits_to_float32(arr, nbits, chunk_samps * nchans)
                chunk[idx, :, :] = flat.reshape(chunk_samps, nchans).T
            list(io_pool.map(_unpack, list(enumerate(packed_list))))
            t_acc["unpack"] += time.perf_counter() - ts

            # --- 3. Whitening stats per beam (parallel).
            ts = time.perf_counter()
            stats_out = [None] * B
            def _stats(bi):
                m, s = stats_fn(chunk[bi])
                stats_out[bi] = (m, s)
            list(io_pool.map(_stats, range(B)))
            for bi in range(B):
                m, s = stats_out[bi]
                if not state_init:
                    state_mean[bi]  = m
                    state_sigma[bi] = s
                else:
                    state_mean[bi]  = (1.0 - alpha) * state_mean[bi]  + alpha * m
                    state_sigma[bi] = (1.0 - alpha) * state_sigma[bi] + alpha * s
            state_init = True
            np.maximum(state_sigma, 1e-3, out=state_sigma)
            t_acc["stats"] += time.perf_counter() - ts

            # --- 4. Whiten in place (parallel across beams).
            ts = time.perf_counter()
            _parallel_whiten(chunk, state_mean, state_sigma, io_pool, params.cores)
            t_acc["whiten"] += time.perf_counter() - ts

            # --- 5. Eigen-clean in place (parallel freq batches).
            ts = time.perf_counter()
            diag = clean_chunk(chunk, clean_params,
                               pool=io_pool, n_workers=params.cores)
            t_acc["clean"]   += time.perf_counter() - ts
            t_acc["cov"]     += diag.t_cov
            t_acc["eigh"]    += diag.t_eigh
            t_acc["project"] += diag.t_project
            t_acc["reshape"] += diag.t_reshape

            total_rfi      += diag.n_rfi_eigvals
            total_rfi_tiles += diag.n_tiles_with_rfi
            total_tiles    += diag.n_tiles

            # --- 6. De-whiten in place (parallel).
            ts = time.perf_counter()
            _parallel_dewhiten(chunk, state_mean, state_sigma, io_pool, params.cores)
            t_acc["dewhiten"] += time.perf_counter() - ts

            # --- 7. Per-beam re-quantize + pack + write (parallel).
            packed_outs: List[Optional[np.ndarray]] = [None] * B
            ts = time.perf_counter()
            def _pack(idx):
                beam = chunk[idx]
                samples = np.ascontiguousarray(beam.T).reshape(-1)
                packed_outs[idx] = pack_float32_to_bits(samples, out_nbits)
            list(io_pool.map(_pack, range(B)))
            t_acc["pack"] += time.perf_counter() - ts

            ts = time.perf_counter()
            def _write(idx):
                beams[idx].fp_out.write(packed_outs[idx].tobytes())
            list(io_pool.map(_write, range(B)))
            t_acc["write"] += time.perf_counter() - ts

            chunks_since_report += 1

            # --- 8. Diagnostics + progress ---
            if diag_fp is not None:
                diag_fp.write(json.dumps({
                    "chunk": ck,
                    "t_sec": ck * chunk_samps * tsamp,
                    "n_tiles": diag.n_tiles,
                    "n_tiles_skipped": diag.n_tiles_skipped,
                    "n_rfi_eigvals": diag.n_rfi_eigvals,
                    "n_tiles_with_rfi": diag.n_tiles_with_rfi,
                    "max_eigval": diag.max_eigval,
                    "mp_edge": diag.mp_edge,
                    "safety_threshold": diag.safety_threshold,
                }) + "\n")

            now = time.monotonic()
            if (not quiet) and (now - last_report > 5.0 or ck == n_chunks - 1):
                done = (ck + 1) / n_chunks
                elapsed = now - t0
                eta = elapsed * (1.0 / max(done, 1e-9) - 1.0)
                tiles_rfi_pct = 100.0 * diag.n_tiles_with_rfi / max(1, diag.n_tiles)
                tiles_skip_pct = 100.0 * diag.n_tiles_skipped / max(1, diag.n_tiles)
                print(f"  chunk {ck+1:6d}/{n_chunks}  "
                      f"{100*done:5.1f}%  elapsed {elapsed/60:6.1f} min  "
                      f"ETA {eta/60:6.1f} min  "
                      f"skipped {tiles_skip_pct:5.1f}%  "
                      f"rfi {tiles_rfi_pct:5.1f}%  "
                      f"max_lam={diag.max_eigval:6.2f}", flush=True)
                # Per-stage timing breakdown, averaged since last report
                n = max(1, chunks_since_report)
                stages = ", ".join(f"{k}={t_acc[k]/n*1e3:6.1f}ms"
                                   for k in ("read", "unpack", "stats", "whiten",
                                             "clean", "dewhiten", "pack", "write"))
                substages = ", ".join(f"{k}={t_acc[k]/n*1e3:6.1f}ms"
                                      for k in ("cov", "eigh", "project", "reshape"))
                print(f"    per-chunk avg: {stages}", flush=True)
                print(f"    clean sub    : {substages}", flush=True)
                for k in t_acc:
                    t_acc[k] = 0.0
                chunks_since_report = 0
                last_report = now
    finally:
        io_pool.shutdown(wait=True)
        if diag_fp: diag_fp.close()
        for b in beams:
            try: b.fp_in.close()
            except Exception: pass
            try:
                if b.fp_out: b.fp_out.close()
            except Exception: pass

    return {
        "n_chunks": n_chunks,
        "nsamps_per_beam": n_chunks * chunk_samps,
        "tobs_sec": n_chunks * chunk_samps * tsamp,
        "total_tiles": total_tiles,
        "total_rfi_eigvals": total_rfi,
        "total_tiles_with_rfi": total_rfi_tiles,
        "elapsed_sec": time.monotonic() - t0,
    }


# ==========================================================================
# Benchmark mode: pure compute, no IO
# ==========================================================================

def _run_benchmark(beams: List[BeamFile], params: PipelineParams,
                   quiet: bool = False) -> int:
    """
    Run the in-memory compute path (unpack -> stats -> whiten -> clean ->
    de-whiten -> pack) for N_BENCH chunks on synthetic data, with no file
    IO whatsoever. Reports per-stage timing so cluster users can see whether
    they're CPU-bound or IO-bound.
    """
    nchans = beams[0].header.get("nchans")
    nbits  = beams[0].header.get("nbits")
    tsamp  = float(beams[0].header.get("tsamp"))
    B = len(beams)
    out_nbits = params.out_nbits if params.out_nbits is not None else nbits

    chunk_samps = int(round(params.chunk_sec / tsamp))
    tile_time_samps = max(1, int(round(params.tile_time_ms * 1e-3 / tsamp)))
    if chunk_samps % tile_time_samps:
        chunk_samps = (chunk_samps // tile_time_samps) * tile_time_samps
    chunk_samps = max(chunk_samps, tile_time_samps)
    if (chunk_samps * nchans * nbits) % 8 != 0:
        chunk_samps = ((chunk_samps * nchans * nbits + 7) // 8 * 8) // (nchans * nbits)

    bytes_per_chunk = (chunk_samps * nchans * nbits) // 8
    n_bench = 4   # run 4 chunks

    for line in _threading_banner(params.cores):
        print(line)
    print(f"[bench] B={B} nchans={nchans} nbits={nbits} tsamp={tsamp*1e6:.2f}us")
    print(f"[bench] chunk_samps={chunk_samps} ({chunk_samps*tsamp:.3f}s)")
    print(f"[bench] chunk buffer = {B*nchans*chunk_samps*4/(1<<30):.2f} GB float32")

    # Pre-allocate everything as the real pipeline would.
    chunk = np.empty((B, nchans, chunk_samps), dtype=np.float32)
    state_mean  = np.zeros((B, nchans), dtype=np.float32)
    state_sigma = np.ones ((B, nchans), dtype=np.float32)

    clean_params = CleanParams(
        tile_time_samp = tile_time_samps,
        tile_freq_chans = max(1, params.tile_freq_chans),
        min_support = params.min_support,
        safety = params.safety,
        max_rank = params.max_rank,
        freq_batch_tiles = params.freq_batch_tiles,
        prefilter = params.prefilter,
    )
    if params.robust_stats == "mad":
        stats_fn = _beam_stats_mad
    elif params.robust_stats == "trimmed":
        stats_fn = _beam_stats_trimmed
    else:
        stats_fn = _beam_stats_plain
    io_pool = ThreadPoolExecutor(max_workers=min(B, max(2, params.cores)))

    rng = np.random.default_rng(0)
    # synthetic packed bytes once, reuse
    packed_one = rng.integers(0, 256, size=bytes_per_chunk, dtype=np.uint8).tobytes()
    packed_list = [packed_one for _ in range(B)]

    t_acc = {"unpack": 0.0, "stats": 0.0, "whiten": 0.0,
             "clean": 0.0, "dewhiten": 0.0, "pack": 0.0,
             "cov": 0.0, "eigh": 0.0, "project": 0.0, "reshape": 0.0}
    t_total = 0.0
    for ck in range(n_bench):
        t_chunk = time.perf_counter()

        ts = time.perf_counter()
        def _unpack(arg):
            idx, packed = arg
            arr = np.frombuffer(packed, dtype=np.uint8)
            flat = unpack_bits_to_float32(arr, nbits, chunk_samps * nchans)
            chunk[idx, :, :] = flat.reshape(chunk_samps, nchans).T
        list(io_pool.map(_unpack, list(enumerate(packed_list))))
        t_acc["unpack"] += time.perf_counter() - ts

        ts = time.perf_counter()
        stats_out = [None] * B
        def _stats(bi):
            m, s = stats_fn(chunk[bi]); stats_out[bi] = (m, s)
        list(io_pool.map(_stats, range(B)))
        for bi in range(B):
            m, s = stats_out[bi]
            if ck == 0:
                state_mean[bi]  = m
                state_sigma[bi] = s
        np.maximum(state_sigma, 1e-3, out=state_sigma)
        t_acc["stats"] += time.perf_counter() - ts

        ts = time.perf_counter()
        _parallel_whiten(chunk, state_mean, state_sigma, io_pool, params.cores)
        t_acc["whiten"] += time.perf_counter() - ts

        ts = time.perf_counter()
        diag = clean_chunk(chunk, clean_params,
                           pool=io_pool, n_workers=params.cores)
        t_acc["clean"]   += time.perf_counter() - ts
        t_acc["cov"]     += diag.t_cov
        t_acc["eigh"]    += diag.t_eigh
        t_acc["project"] += diag.t_project
        t_acc["reshape"] += diag.t_reshape

        ts = time.perf_counter()
        _parallel_dewhiten(chunk, state_mean, state_sigma, io_pool, params.cores)
        t_acc["dewhiten"] += time.perf_counter() - ts

        ts = time.perf_counter()
        def _pack(idx):
            beam = chunk[idx]
            samples = np.ascontiguousarray(beam.T).reshape(-1)
            return pack_float32_to_bits(samples, out_nbits)
        list(io_pool.map(_pack, range(B)))
        t_acc["pack"] += time.perf_counter() - ts

        chunk_elapsed = time.perf_counter() - t_chunk
        t_total += chunk_elapsed
        print(f"[bench] chunk {ck+1}/{n_bench}: {chunk_elapsed*1000:.1f} ms  "
              f"(skipped {100.0*diag.n_tiles_skipped/max(1,diag.n_tiles):.1f}% of tiles)",
              flush=True)

    io_pool.shutdown(wait=True)
    print(f"\n[bench] per-chunk stage averages ({n_bench} chunks):")
    for k, v in t_acc.items():
        print(f"  {k:10s}: {v/n_bench*1000:7.1f} ms")
    print(f"  TOTAL     : {t_total/n_bench*1000:7.1f} ms")
    rt = (n_bench * chunk_samps * tsamp) / t_total
    print(f"\n[bench] pure-compute throughput: {rt:.2f}x real-time on {params.cores} cores")
    return 0


# ==========================================================================
# CLI
# ==========================================================================

def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="+",
                   help="SIGPROC filterbank files, one per beam (>=1)")
    p.add_argument("--outdir", type=Path, default=Path("cleaned"))
    p.add_argument("--suffix", default="_clean")
    p.add_argument("--out-nbits", type=int, default=None,
                   help="output bit depth (default: same as input)")
    p.add_argument("--cores", "-j", type=int, default=1,
                   help="number of CPU cores for BLAS + parallel IO (default: 1)")
    p.add_argument("--chunk-sec", type=float, default=1.0,
                   help="processing chunk length in seconds; the dominant "
                        "memory user (chunk buffer is B*nchans*chunk_samps*4 "
                        "bytes). Drop this if you OOM/swap. Default: 1.0")
    p.add_argument("--whitening-tau-sec", type=float, default=5.0,
                   help="EMA timescale for whitening stats (default: 5.0 s)")
    p.add_argument("--tile-time-ms", type=float, default=100.0,
                   help="time tile size in ms for eigen-clean (default: 100)")
    p.add_argument("--tile-freq-chans", type=int, default=4,
                   help="freq tile size in channels (default: 4)")
    p.add_argument("--freq-batch-tiles", type=int, default=8,
                   help="number of freq tiles batched per einsum (default: 8); "
                        "increase for fewer BLAS calls at the cost of more "
                        "scratch, decrease for lower memory")
    p.add_argument("--min-support", type=float, default=20.0,
                   help="participation-ratio threshold for eigenvector spatial "
                        "support (default: 20)")
    p.add_argument("--safety", type=float, default=1.10,
                   help="multiplicative safety factor on MP edge (default: 1.10)")
    p.add_argument("--max-rank", type=int, default=5,
                   help="cap on # eigenvectors subtracted per tile (default: 5)")
    p.add_argument("--robust-stats", choices=("none", "trimmed", "mad"),
                   default="none",
                   help="per-beam, per-channel whitening estimator. "
                        "'none' (default) is plain mean/std, cheapest. "
                        "'trimmed' iteratively sigma-clips. "
                        "'mad' is textbook MAD but heaviest. "
                        "The eigen analysis only needs whitening to be "
                        "consistent across beams, so 'none' is usually fine.")
    p.add_argument("--no-prefilter", action="store_true",
                   help="disable the cheap trace-based pre-filter that skips "
                        "tiles which cannot possibly host an RFI eigenvalue. "
                        "Pre-filter avoids ~99%% of eigh calls on typical data.")
    p.add_argument("--tstart-tol-sec", type=float, default=1e-3,
                   help="tolerance on tstart agreement across beams (default: 1 ms)")
    p.add_argument("--diag", type=Path, default=None,
                   help="optional path for per-chunk JSON-lines diagnostics")
    p.add_argument("--dry-run", action="store_true",
                   help="parse headers, validate compatibility, print plan and exit")
    p.add_argument("--benchmark", action="store_true",
                   help="ignore input files; generate B=72 (or len(inputs)) "
                        "synthetic beams in RAM and run a 4-second smoke through "
                        "the pipeline to time pure compute. Useful for "
                        "diagnosing slow IO vs slow BLAS on a cluster node.")
    p.add_argument("--quiet", action="store_true")
    args = p.parse_args(argv)

    if args.cores < 1:
        print("multibeam_clean: --cores must be >= 1", file=sys.stderr)
        return 2

    input_paths = [Path(s) for s in args.inputs]
    for ip in input_paths:
        if not ip.exists():
            print(f"multibeam_clean: input not found: {ip}", file=sys.stderr)
            return 2

    beams = [open_input(ip) for ip in input_paths]
    validate_compatible(beams, args.tstart_tol_sec)

    nchans = beams[0].header.get("nchans")
    nbits  = beams[0].header.get("nbits")
    tsamp  = float(beams[0].header.get("tsamp"))
    out_nbits = args.out_nbits if args.out_nbits is not None else nbits

    params = PipelineParams(
        chunk_sec = args.chunk_sec,
        whitening_tau_sec = args.whitening_tau_sec,
        tile_time_ms = args.tile_time_ms,
        tile_freq_chans = args.tile_freq_chans,
        min_support = args.min_support,
        safety = args.safety,
        max_rank = args.max_rank,
        out_nbits = out_nbits,
        suffix = args.suffix,
        tstart_tol_sec = args.tstart_tol_sec,
        diag_path = args.diag,
        cores = args.cores,
        freq_batch_tiles = args.freq_batch_tiles,
        robust_stats = args.robust_stats,
        prefilter = not args.no_prefilter,
    )

    if args.benchmark:
        return _run_benchmark(beams, params, quiet=args.quiet)

    if args.dry_run:
        chunk_samps = int(round(args.chunk_sec / tsamp))
        if chunk_samps == 0: chunk_samps = 1
        chunk_gb = (len(beams) * nchans * chunk_samps * 4) / (1 << 30)
        tobs = min(b.nsamples_total for b in beams) * tsamp
        print(f"[dry-run] beams              : {len(beams)}")
        for b in beams:
            print(f"           {b.path.name:40s}  "
                  f"{b.nsamples_total:>12d} samp  "
                  f"({b.nsamples_total*tsamp:8.2f} s)")
        print(f"[dry-run] nchans / nbits      : {nchans} / {nbits} "
              f"(out_nbits={out_nbits})")
        print(f"[dry-run] tsamp               : {tsamp*1e6:.3f} us")
        print(f"[dry-run] usable tobs         : {tobs:.2f} s ({tobs/3600:.3f} hr)")
        print(f"[dry-run] chunk_sec           : {args.chunk_sec}")
        print(f"[dry-run] chunk buffer        : {chunk_gb:.2f} GB float32")
        print(f"[dry-run] outdir              : {args.outdir}")
        print(f"[dry-run] cores               : {args.cores}")
        for b in beams: b.fp_in.close()
        return 0

    open_outputs(beams, args.outdir, args.suffix, args.out_nbits)
    summary = process(beams, params, quiet=args.quiet)

    if not args.quiet:
        print(f"\n[multibeam_clean] DONE. processed {summary['tobs_sec']:.2f} s "
              f"in {summary['elapsed_sec']/60:.2f} min "
              f"(rt-factor {summary['tobs_sec']/max(summary['elapsed_sec'],1e-9):.2f}x)")
        print(f"[multibeam_clean] tiles_with_rfi: "
              f"{summary['total_tiles_with_rfi']} / {summary['total_tiles']} "
              f"({100.0*summary['total_tiles_with_rfi']/max(1,summary['total_tiles']):.2f}%)")
        for b in beams:
            print(f"[multibeam_clean] wrote {b.out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
