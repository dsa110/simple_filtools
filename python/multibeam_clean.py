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
  5. Subtract that subspace:   x_clean = x - U_I U_I^T x
  6. Re-quantize to the input bit-depth and write.

Reads everything (nchans, nbits, fch1, foff, tsamp, tstart, ...) from
the SIGPROC headers; only requires that all input files agree on the
header fields that define the data layout (nchans, nbits, nifs,
fch1, foff, tsamp). Tstart is required to be within a configurable
tolerance.

Dependencies: numpy only (plus Python >= 3.8 stdlib).

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


# ==========================================================================
# SIGPROC header I/O (preserves field order so the writer round-trips)
# ==========================================================================

# Known tag types. Anything not here is rejected -- we don't want to silently
# pass through fields we can't faithfully re-serialize.
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
        # not present -- append with inferred type
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
    """Write the header back. Returns total bytes written (= new header_bytes)."""
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
# ==========================================================================

def unpack_bits_to_float32(packed: np.ndarray, nbits: int,
                           nsamps_total: int) -> np.ndarray:
    """
    `packed`: uint8 array of length ceil(nsamps_total * nbits / 8).
    Returns a float32 array of length `nsamps_total` containing the samples,
    in SIGPROC order (4 samples per byte for nbits=2, low bits first).
    """
    if nbits == 8:
        return packed[:nsamps_total].astype(np.float32, copy=True)
    if nbits == 16:
        v = packed.view(np.int16)[:nsamps_total].astype(np.float32)
        return v
    if nbits == 32:
        v = packed.view(np.float32)[:nsamps_total].astype(np.float32, copy=True)
        return v
    if nbits == 4:
        b = packed.astype(np.uint8, copy=False)
        lo = (b & 0x0F).astype(np.float32)
        hi = ((b >> 4) & 0x0F).astype(np.float32)
        out = np.empty(b.size * 2, dtype=np.float32)
        out[0::2] = lo; out[1::2] = hi
        return out[:nsamps_total]
    if nbits == 2:
        b = packed.astype(np.uint8, copy=False)
        s0 = (b      & 0x3).astype(np.float32)
        s1 = ((b>>2) & 0x3).astype(np.float32)
        s2 = ((b>>4) & 0x3).astype(np.float32)
        s3 = ((b>>6) & 0x3).astype(np.float32)
        out = np.empty(b.size * 4, dtype=np.float32)
        out[0::4] = s0; out[1::4] = s1; out[2::4] = s2; out[3::4] = s3
        return out[:nsamps_total]
    if nbits == 1:
        b = packed.astype(np.uint8, copy=False)
        out = np.empty(b.size * 8, dtype=np.float32)
        for k in range(8):
            out[k::8] = ((b >> k) & 0x1).astype(np.float32)
        return out[:nsamps_total]
    raise ValueError(f"unsupported nbits: {nbits}")


def pack_float32_to_bits(samples: np.ndarray, nbits: int) -> np.ndarray:
    """
    Re-quantize float32 samples (already mapped to the correct integer range)
    to `nbits`-packed uint8 bytes. SIGPROC sample order; for nbits<8 the
    length must be a multiple of (8//nbits).
    """
    if nbits == 8:
        clipped = np.clip(samples, 0, 255).astype(np.uint8)
        return clipped
    if nbits == 16:
        return np.clip(samples, -32768, 32767).astype(np.int16).view(np.uint8)
    if nbits == 32:
        return samples.astype(np.float32).view(np.uint8)
    if nbits == 4:
        q = np.clip(np.rint(samples), 0, 15).astype(np.uint8)
        if q.size % 2: raise ValueError("4-bit pack needs even length")
        return (q[0::2] | (q[1::2] << 4)).astype(np.uint8)
    if nbits == 2:
        q = np.clip(np.rint(samples), 0, 3).astype(np.uint8)
        if q.size % 4: raise ValueError("2-bit pack needs length % 4 == 0")
        return (q[0::4] | (q[1::4] << 2) | (q[2::4] << 4) | (q[3::4] << 6)).astype(np.uint8)
    if nbits == 1:
        q = np.clip(np.rint(samples), 0, 1).astype(np.uint8)
        if q.size % 8: raise ValueError("1-bit pack needs length % 8 == 0")
        out = np.zeros(q.size // 8, dtype=np.uint8)
        for k in range(8):
            out |= (q[k::8] << k)
        return out
    raise ValueError(f"unsupported nbits: {nbits}")


def bytes_per_sample(nbits: int, nchans: int) -> float:
    return (nbits * nchans) / 8.0


# ==========================================================================
# Whitening: running per-(beam, channel) robust mean / scale
# ==========================================================================

@dataclass
class WhiteningState:
    """One state per beam. Holds running mean and scale per channel."""
    mean: np.ndarray   # (nchans,) float32
    sigma: np.ndarray  # (nchans,) float32
    initialized: bool = False


def _robust_chunk_stats(x: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """x: (C, T). Returns (median per channel, MAD per channel) as float32."""
    med = np.median(x, axis=1)
    mad = np.median(np.abs(x - med[:, None]), axis=1)
    # 1.4826 maps MAD -> sigma for Gaussian. For 2-bit clipped data this is a
    # mild overestimate but consistent and well-defined; the eigen analysis
    # only needs whitening to be consistent across beams, not absolutely
    # correct.
    sigma = 1.4826 * mad
    # avoid divide-by-zero: replace sigma==0 by std fallback
    bad = sigma <= 0
    if bad.any():
        sigma[bad] = np.maximum(x[bad].std(axis=1), 1e-6)
    return med.astype(np.float32), sigma.astype(np.float32)


def update_whitening(state: WhiteningState, chunk: np.ndarray, alpha: float) -> None:
    """In-place EMA update of running mean / sigma using this chunk's robust stats."""
    med, sigma = _robust_chunk_stats(chunk)
    if not state.initialized:
        state.mean = med
        state.sigma = sigma
        state.initialized = True
    else:
        state.mean  = (1.0 - alpha) * state.mean  + alpha * med
        state.sigma = (1.0 - alpha) * state.sigma + alpha * sigma


# ==========================================================================
# Eigen-cleaning kernel (batched over all tiles in a chunk)
# ==========================================================================

@dataclass
class CleanParams:
    tile_time_samp: int     # samples per time tile
    tile_freq_chans: int    # channels per freq tile
    min_support: float      # eigenvector participation ratio threshold
    safety: float           # multiplicative safety on MP edge
    max_rank: int           # cap on # subtracted eigenvecs per tile


@dataclass
class CleanDiag:
    n_tiles: int = 0
    n_rfi_eigvals: int = 0       # total across all tiles in this chunk
    n_tiles_with_rfi: int = 0
    max_eigval: float = 0.0
    mp_edge: float = 0.0
    safety_threshold: float = 0.0


def clean_chunk(whitened: np.ndarray, params: CleanParams) -> Tuple[np.ndarray, CleanDiag]:
    """
    Clean a (B, C, T) whitened chunk. Returns (cleaned float32 array, diagnostics).
    Tiles smaller than (tile_freq_chans, tile_time_samp) at the edges are left
    untouched (passed through).
    """
    B, C, T = whitened.shape
    Fw = params.tile_freq_chans
    Tw = params.tile_time_samp
    nFt = C // Fw
    nTt = T // Tw
    C_use = nFt * Fw
    T_use = nTt * Tw

    diag = CleanDiag()
    cleaned = whitened.copy()

    if nFt == 0 or nTt == 0:
        return cleaned, diag

    # data: (B, C_use, T_use)
    data = whitened[:, :C_use, :T_use]
    # reshape -> (B, nFt, Fw, nTt, Tw) -> (nFt, nTt, B, Fw*Tw)
    data = (data
            .reshape(B, nFt, Fw, nTt, Tw)
            .transpose(1, 3, 0, 2, 4)
            .reshape(nFt, nTt, B, Fw * Tw))
    N = Fw * Tw

    # Per-tile cross-beam covariance: shape (nFt, nTt, B, B)
    cov = np.einsum("FTbn,FTcn->FTbc", data, data, optimize=True) / float(N)

    # Batched eigendecomp. eigvalsh returns ascending order.
    eigvals, eigvecs = np.linalg.eigh(cov)
    eigvals = eigvals[..., ::-1]                  # (nFt, nTt, B), descending
    eigvecs = eigvecs[..., ::-1]                  # (nFt, nTt, B, B), cols = eigvecs

    # Marchenko-Pastur upper edge for a B x N whitened sample covariance
    q = float(B) / float(N)
    mp_edge = (1.0 + math.sqrt(q)) ** 2
    threshold = mp_edge * params.safety

    significant = eigvals > threshold             # (nFt, nTt, B)

    # Participation ratio per eigenvector u_k:  P = 1 / sum_b u_b^4
    # eigvecs[..., :, k] is u_k; reduce over the B axis (-2)
    u_pow4_sum = (eigvecs ** 4).sum(axis=-2)      # (nFt, nTt, B)
    pr = 1.0 / np.maximum(u_pow4_sum, 1e-30)
    enough_support = pr >= float(params.min_support)

    is_rfi = significant & enough_support
    # Cap rank: if more than max_rank significant + supported eigenvecs in a
    # tile, keep only the top max_rank (by eigenvalue, which is sorted desc).
    if params.max_rank > 0 and params.max_rank < B:
        cum = np.cumsum(is_rfi.astype(np.int32), axis=-1)
        is_rfi = is_rfi & (cum <= params.max_rank)

    diag.n_tiles = int(nFt * nTt)
    diag.n_rfi_eigvals = int(is_rfi.sum())
    diag.n_tiles_with_rfi = int(is_rfi.any(axis=-1).sum())
    diag.max_eigval = float(eigvals[..., 0].max()) if eigvals.size else 0.0
    diag.mp_edge = mp_edge
    diag.safety_threshold = threshold

    if diag.n_rfi_eigvals == 0:
        # nothing to subtract: pass data through
        return cleaned, diag

    # Project: x_clean = x - U_I U_I^T x
    #        = U_keep U_keep^T x, where U_keep zeroes out RFI columns of U.
    # We avoid building U_keep explicitly; instead zero the RFI projections.
    keep_mask = (~is_rfi).astype(np.float32)       # (nFt, nTt, B)

    # v = U^T x : shape (nFt, nTt, B_eig, N)
    v = np.einsum("FTbk,FTbn->FTkn", eigvecs, data, optimize=True)
    v *= keep_mask[..., None]
    cleaned_tiles = np.einsum("FTbk,FTkn->FTbn", eigvecs, v, optimize=True)

    # Reshape back and write into output
    cleaned_tiles = (cleaned_tiles
                     .reshape(nFt, nTt, B, Fw, Tw)
                     .transpose(2, 0, 3, 1, 4)
                     .reshape(B, C_use, T_use))
    cleaned[:, :C_use, :T_use] = cleaned_tiles
    return cleaned, diag


# ==========================================================================
# Multi-beam I/O + pipeline
# ==========================================================================

@dataclass
class BeamFile:
    path: Path
    header: SigprocHeader
    fp_in: object = None              # input file handle (open)
    fp_out: object = None             # output file handle (open)
    data_bytes_total: int = 0         # data bytes available in input
    nsamples_total: int = 0           # samples available in input
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
                    f"{b.path.name}: {key}={v} differs from reference "
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
    # tstart within tolerance (days)
    ref_tstart = float(ref.get("tstart", 0.0))
    tol_days = tstart_tol_sec / 86400.0
    for b in beams[1:]:
        ts = float(b.header.get("tstart", 0.0))
        if abs(ts - ref_tstart) > tol_days:
            raise ValueError(
                f"{b.path.name}: tstart {ts} differs from {ref_tstart} by "
                f">{tstart_tol_sec}s tolerance")
    # nifs must be 1 (we don't handle polarized data here)
    if (ref.get("nifs") or 1) != 1:
        raise ValueError("multibeam_clean: only nifs=1 is supported")


def open_outputs(beams: List[BeamFile], outdir: Path, suffix: str,
                 out_nbits: Optional[int]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for b in beams:
        stem = b.path.stem
        b.out_path = outdir / f"{stem}{suffix}.fil"
        b.fp_out = open(b.out_path, "wb")
        # Write header. If user changed output nbits, update the field.
        hdr_to_write = SigprocHeader(fields=list(b.header.fields))
        if out_nbits is not None and out_nbits != b.header.get("nbits"):
            hdr_to_write.set("nbits", int(out_nbits))
        # Update rawdatafile to point at the cleaned file (cosmetic).
        if hdr_to_write.get("rawdatafile") is not None:
            hdr_to_write.set("rawdatafile", b.out_path.name)
        write_sigproc_header(b.fp_out, hdr_to_write)


def read_chunk_packed(beam: BeamFile, n_packed_bytes: int) -> bytes:
    return beam.fp_in.read(n_packed_bytes)


# ==========================================================================
# Main pipeline
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


def process(beams: List[BeamFile], params: PipelineParams,
            quiet: bool = False) -> Dict:
    """Stream all beams, clean, write outputs. Returns summary dict."""
    nchans = beams[0].header.get("nchans")
    nbits  = beams[0].header.get("nbits")
    tsamp  = float(beams[0].header.get("tsamp"))
    B = len(beams)
    out_nbits = params.out_nbits if params.out_nbits is not None else nbits

    chunk_samps = int(round(params.chunk_sec / tsamp))
    tile_time_samps = max(1, int(round(params.tile_time_ms * 1e-3 / tsamp)))

    # Force chunk_samps to be a multiple of tile_time_samps so all tiles inside
    # a chunk are full-size (simpler accounting).
    if chunk_samps % tile_time_samps:
        chunk_samps = (chunk_samps // tile_time_samps) * tile_time_samps
    if chunk_samps == 0:
        raise ValueError("chunk_sec too short for tile_time_ms; pick larger chunk")

    # The number of samples we can actually process is the min across beams.
    nsamps_avail = min(b.nsamples_total for b in beams)
    # For 2-bit packing, the per-beam sample count must align with byte
    # boundaries; chunk_samps is already a multiple of tile_time_samps so as
    # long as (chunk_samps * nchans * nbits) is a multiple of 8 we're fine.
    bytes_per_chunk = int((chunk_samps * nchans * nbits) // 8)
    if (chunk_samps * nchans * nbits) % 8:
        raise ValueError(
            "chunk size does not produce an integer number of bytes per beam; "
            "increase --tile-time-ms or --chunk-sec to align")

    n_chunks = nsamps_avail // chunk_samps

    alpha = float(params.chunk_sec) / float(params.whitening_tau_sec)
    alpha = max(min(alpha, 1.0), 1e-3)

    clean_params = CleanParams(
        tile_time_samp = tile_time_samps,
        tile_freq_chans = max(1, params.tile_freq_chans),
        min_support = params.min_support,
        safety = params.safety,
        max_rank = params.max_rank,
    )

    if not quiet:
        tobs = nsamps_avail * tsamp
        print(f"[multibeam_clean] beams={B}  nchans={nchans}  nbits={nbits}  "
              f"tsamp={tsamp*1e6:.3f} us")
        print(f"[multibeam_clean] usable samples per beam = {nsamps_avail}  "
              f"({tobs:.2f} s = {tobs/3600:.3f} hr)")
        print(f"[multibeam_clean] chunk = {chunk_samps} samples ({chunk_samps*tsamp:.4f} s);  "
              f"{n_chunks} chunks total")
        print(f"[multibeam_clean] tile_time = {tile_time_samps} samples "
              f"({tile_time_samps*tsamp*1e3:.3f} ms);  tile_freq = "
              f"{clean_params.tile_freq_chans} chans")
        print(f"[multibeam_clean] whitening tau = {params.whitening_tau_sec:.3f} s (alpha = {alpha:.4f})")
        print(f"[multibeam_clean] min_support = {params.min_support}  "
              f"safety = {params.safety}  max_rank = {params.max_rank}  "
              f"out_nbits = {out_nbits}")
        print(f"[multibeam_clean] cores = {params.cores}")

    # Whitening state per beam
    whit = [WhiteningState(np.zeros(nchans, np.float32),
                           np.ones (nchans, np.float32),
                           initialized=False) for _ in range(B)]

    diag_fp = open(params.diag_path, "w") if params.diag_path else None

    # Reuse a parallel thread pool for IO across beams
    io_pool = ThreadPoolExecutor(max_workers=min(B, max(2, params.cores)))

    t0 = time.monotonic()
    last_report = t0
    total_rfi = 0
    total_rfi_tiles = 0
    total_tiles = 0

    try:
        for ck in range(n_chunks):
            # 1) Read packed bytes for every beam in parallel
            def _read(b):
                return read_chunk_packed(b, bytes_per_chunk)
            packed_list = list(io_pool.map(_read, beams))
            if any(len(p) != bytes_per_chunk for p in packed_list):
                if not quiet:
                    print(f"[multibeam_clean] short read at chunk {ck}, stopping", file=sys.stderr)
                break

            # 2) Unpack to float32 (B, C, T) in parallel across beams
            def _unpack(arg):
                idx, packed = arg
                arr = unpack_bits_to_float32(
                    np.frombuffer(packed, dtype=np.uint8),
                    nbits, chunk_samps * nchans
                )
                # SIGPROC layout is sample-major over chans: spectrum after
                # spectrum, each spectrum has nchans values. So reshape:
                arr = arr.reshape(chunk_samps, nchans).T  # -> (C, T)
                return idx, arr
            unpacked = [None] * B
            for idx, arr in io_pool.map(_unpack, list(enumerate(packed_list))):
                unpacked[idx] = arr
            chunk = np.stack(unpacked, axis=0)  # (B, C, T) float32

            # 3) Update whitening state and whiten in place
            for bi in range(B):
                update_whitening(whit[bi], chunk[bi], alpha)
            mean_bf  = np.stack([w.mean  for w in whit], axis=0)  # (B, C)
            sigma_bf = np.stack([w.sigma for w in whit], axis=0)  # (B, C)
            sigma_bf = np.maximum(sigma_bf, 1e-3).astype(np.float32)
            whitened = (chunk - mean_bf[..., None]) / sigma_bf[..., None]

            # 4) Clean
            cleaned, diag = clean_chunk(whitened, clean_params)
            total_rfi      += diag.n_rfi_eigvals
            total_rfi_tiles += diag.n_tiles_with_rfi
            total_tiles    += diag.n_tiles

            # 5) De-whiten and re-quantize back to the input dynamic range
            out_floats = cleaned * sigma_bf[..., None] + mean_bf[..., None]

            # 6) Per-beam pack + write (in parallel)
            #    SIGPROC layout: time-major over chans (one spectrum after
            #    another), so transpose (C, T) -> (T, C) and flatten.
            def _pack_write(arg):
                idx, b = arg
                samples = out_floats[idx].T.reshape(-1)  # (T*C,)
                if out_nbits in (8, 4, 2, 1):
                    samples = np.rint(samples)
                packed = pack_float32_to_bits(samples, out_nbits)
                b.fp_out.write(packed.tobytes())
            list(io_pool.map(_pack_write, [(i, b) for i, b in enumerate(beams)]))

            # 7) Diagnostics + progress
            if diag_fp is not None:
                diag_fp.write(json.dumps({
                    "chunk": ck,
                    "t_sec": ck * chunk_samps * tsamp,
                    "n_tiles": diag.n_tiles,
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
                print(f"  chunk {ck+1:6d}/{n_chunks}  "
                      f"{100*done:5.1f}%  elapsed {elapsed/60:6.1f} min  "
                      f"ETA {eta/60:6.1f} min  "
                      f"tiles_with_rfi {tiles_rfi_pct:5.1f}%  "
                      f"max_lambda={diag.max_eigval:6.2f}", flush=True)
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
# CLI
# ==========================================================================

def estimate_runtime_sec(nbeams: int, nchans: int, nbits: int,
                         tobs_sec: float, cores: int) -> float:
    """Very rough wall-clock estimate; useful for --dry-run."""
    # IO read+write floor: bytes / 2 GB/s
    bytes_one_beam = nchans * nbits / 8.0 * tobs_sec * (1.0 / 50e-6) \
        if tobs_sec > 0 else 0
    total_bytes = 2 * nbeams * bytes_one_beam  # read + write
    io_sec = total_bytes / (2e9)
    # Compute estimate: B^2 * N_tiles * 1e-6 / cores
    return max(io_sec, io_sec * 0.5 + 60.0)  # crude lower bound


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="+",
                   help="SIGPROC filterbank files, one per beam (>=1)")
    p.add_argument("--outdir", type=Path, default=Path("cleaned"),
                   help="output directory (default: ./cleaned)")
    p.add_argument("--suffix", default="_clean",
                   help="suffix appended to each input stem (default: _clean)")
    p.add_argument("--out-nbits", type=int, default=None,
                   help="output bit depth (default: same as input)")
    p.add_argument("--cores", "-j", type=int, default=1,
                   help="number of CPU cores for BLAS + parallel IO (default: 1)")
    p.add_argument("--chunk-sec", type=float, default=1.0,
                   help="processing chunk length in seconds (default: 1.0)")
    p.add_argument("--whitening-tau-sec", type=float, default=5.0,
                   help="EMA timescale for whitening stats (default: 5.0 s)")
    p.add_argument("--tile-time-ms", type=float, default=100.0,
                   help="time tile size in ms for eigen-clean (default: 100)")
    p.add_argument("--tile-freq-chans", type=int, default=4,
                   help="freq tile size in channels (default: 4)")
    p.add_argument("--min-support", type=float, default=20.0,
                   help="participation-ratio threshold for eigenvector spatial "
                        "support; only subtract eigvecs whose effective beam "
                        "count >= this (default: 20)")
    p.add_argument("--safety", type=float, default=1.10,
                   help="multiplicative safety factor on Marchenko-Pastur edge "
                        "(default: 1.10)")
    p.add_argument("--max-rank", type=int, default=5,
                   help="cap on # eigenvectors subtracted per tile (default: 5)")
    p.add_argument("--tstart-tol-sec", type=float, default=1e-3,
                   help="tolerance on tstart agreement across beams (default: 1 ms)")
    p.add_argument("--diag", type=Path, default=None,
                   help="optional path for per-chunk JSON-lines diagnostics")
    p.add_argument("--dry-run", action="store_true",
                   help="parse headers, validate compatibility, print plan and exit")
    p.add_argument("--quiet", action="store_true",
                   help="suppress per-chunk progress")
    args = p.parse_args(argv)

    if args.cores < 1:
        print("multibeam_clean: --cores must be >= 1", file=sys.stderr)
        return 2

    input_paths = [Path(s) for s in args.inputs]
    for ip in input_paths:
        if not ip.exists():
            print(f"multibeam_clean: input not found: {ip}", file=sys.stderr)
            return 2

    # Open + read headers
    beams = [open_input(ip) for ip in input_paths]
    validate_compatible(beams, args.tstart_tol_sec)

    nchans = beams[0].header.get("nchans")
    nbits  = beams[0].header.get("nbits")
    tsamp  = float(beams[0].header.get("tsamp"))
    out_nbits = args.out_nbits if args.out_nbits is not None else nbits

    # Build pipeline params
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
    )

    if args.dry_run:
        tobs = min(b.nsamples_total for b in beams) * tsamp
        rough = estimate_runtime_sec(len(beams), nchans, nbits, tobs, args.cores)
        print(f"[dry-run] beams              : {len(beams)}")
        for b in beams:
            print(f"           {b.path.name:40s}  "
                  f"{b.nsamples_total:>12d} samp  "
                  f"({b.nsamples_total*tsamp:8.2f} s)")
        print(f"[dry-run] nchans / nbits      : {nchans} / {nbits} "
              f"(out_nbits={out_nbits})")
        print(f"[dry-run] tsamp               : {tsamp*1e6:.3f} us")
        print(f"[dry-run] usable tobs         : {tobs:.2f} s ({tobs/3600:.3f} hr)")
        print(f"[dry-run] outdir              : {args.outdir}")
        print(f"[dry-run] cores               : {args.cores}")
        print(f"[dry-run] crude runtime floor : {rough/60:.1f} min")
        for b in beams: b.fp_in.close()
        return 0

    # Open outputs
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
