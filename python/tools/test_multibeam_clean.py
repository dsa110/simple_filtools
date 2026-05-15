#!/usr/bin/env python3
"""Synthetic smoke test for multibeam_clean.py.

Builds B=24 fake 2-bit filterbanks (small) where the first M=20 beams
share a common-mode burst, then runs the cleaner and checks that the
common-mode is suppressed in those beams while clean beams are
preserved.

Run from anywhere:
    python3 simple_filtools/python/tools/test_multibeam_clean.py
"""
from __future__ import annotations

import os
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]            # simple_filtools/python
TOOLS = ROOT / "tools"
sys.path.insert(0, str(ROOT))
import multibeam_clean as MBC


def _write_str(buf, s):
    b = s.encode("ascii")
    buf += struct.pack("<i", len(b)); buf += b


def _build_hdr(nchans, nbits, fch1, foff, tstart, tsamp, name):
    buf = bytearray()
    _write_str(buf, "HEADER_START")
    _write_str(buf, "rawdatafile"); _write_str(buf, name)
    _write_str(buf, "source_name"); _write_str(buf, "TEST")
    for tag, val in [("machine_id", 1), ("telescope_id", 1), ("data_type", 1),
                     ("nchans", nchans), ("nbits", nbits), ("nifs", 1),
                     ("nbeams", 1), ("ibeam", 1)]:
        _write_str(buf, tag); buf += struct.pack("<i", val)
    for tag, val in [("fch1", fch1), ("foff", foff),
                     ("tstart", tstart), ("tsamp", tsamp)]:
        _write_str(buf, tag); buf += struct.pack("<d", val)
    _write_str(buf, "HEADER_END")
    return bytes(buf)


def make_beam(path, nchans, nsamps, tsamp, fch1, foff, tstart, samples_2bit, name):
    hdr = _build_hdr(nchans, 2, fch1, foff, tstart, tsamp, name)
    with open(path, "wb") as f:
        f.write(hdr)
        # samples_2bit shape: (nsamps, nchans), values in {0..3}
        flat = samples_2bit.reshape(-1).astype(np.uint8)
        packed = (flat[0::4] | (flat[1::4] << 2) | (flat[2::4] << 4) | (flat[3::4] << 6)).astype(np.uint8)
        f.write(packed.tobytes())


def main():
    rng = np.random.default_rng(123)
    nchans = 64
    nsamps = 8000            # 8000 spectra
    tsamp  = 500e-6          # 500 us
    fch1, foff = 1500.0, -1.0
    tstart = 60000.0

    B_total = 24
    B_rfi   = 20             # first 20 beams share common-mode
    rfi_t   = np.array([1200, 1201, 1202, 3000, 5400, 5401])  # planted burst times

    workdir = Path(tempfile.mkdtemp(prefix="mbc_test_"))
    print(f"[test] workdir = {workdir}")
    in_paths = []
    for bi in range(B_total):
        # iid uniform {0..3} noise
        samp = rng.integers(0, 4, size=(nsamps, nchans), dtype=np.uint8)
        if bi < B_rfi:
            # add a bright spectrum at the burst times: force all channels to 3
            samp[rfi_t, :] = 3
        p = workdir / f"beam{bi:02d}.fil"
        make_beam(p, nchans, nsamps, tsamp, fch1, foff, tstart, samp,
                  name=f"beam{bi:02d}.fil")
        in_paths.append(p)

    outdir = workdir / "cleaned"
    argv = [str(p) for p in in_paths] + [
        "--outdir", str(outdir),
        "--cores", "2",
        "--chunk-sec", "1.0",
        "--tile-time-ms", "10.0",      # tile_time = 20 samples at 500us
        "--tile-freq-chans", "4",
        "--min-support", "15",
        "--safety", "1.05",
        "--max-rank", "3",
        "--whitening-tau-sec", "1.0",
        "--diag", str(workdir / "diag.jsonl"),
        "--quiet",
    ]
    rc = MBC.main(argv)
    assert rc == 0, f"multibeam_clean returned {rc}"

    # Re-read the cleaned beams; measure suppression at planted burst times.
    def read_beam(p):
        with open(p, "rb") as f:
            hdr = MBC.read_sigproc_header(f)
            data = f.read()
        nb = hdr.get("nbits")
        nc = hdr.get("nchans")
        n_total = nsamps * nc
        flat = MBC.unpack_bits_to_float32(np.frombuffer(data, dtype=np.uint8),
                                          nb, n_total)
        return flat.reshape(nsamps, nc)

    # zero-DM time series (mean over channels) -- planted bursts will spike
    # in the input but should be suppressed in the output for the RFI beams.
    input_zdm  = np.zeros((B_total, nsamps))
    output_zdm = np.zeros((B_total, nsamps))
    for bi, ip in enumerate(in_paths):
        op = outdir / f"{ip.stem}_clean.fil"
        assert op.exists(), f"missing output: {op}"
        x = read_beam(ip)
        y = read_beam(op)
        input_zdm[bi]  = x.mean(axis=1)
        output_zdm[bi] = y.mean(axis=1)

    # Median input zdm at the burst samples for the RFI beams should be
    # near 3 (we set all channels to 3); the median over RFI beams in the
    # CLEANED output should be near the noise mean (~1.5).
    rfi_in  = np.median(input_zdm[:B_rfi][:, rfi_t])
    rfi_out = np.median(output_zdm[:B_rfi][:, rfi_t])
    # Non-RFI beams: cleaning should not boost or kill the noise floor.
    noise_in  = input_zdm[B_rfi:].mean()
    noise_out = output_zdm[B_rfi:].mean()

    print(f"[test] RFI-beam zdm at burst: input median = {rfi_in:.3f}, "
          f"output median = {rfi_out:.3f}  (expect ~3 -> ~1.5)")
    print(f"[test] non-RFI mean zdm:      input = {noise_in:.3f}, "
          f"output = {noise_out:.3f}  (expect both ~1.5)")

    # Pass criteria (loose to allow re-quantization noise):
    #  - planted spikes suppressed: rfi_out << rfi_in
    #  - non-RFI beams roughly unchanged in mean
    assert rfi_out < rfi_in - 0.5, \
        f"planted RFI not suppressed enough: {rfi_in:.2f} -> {rfi_out:.2f}"
    assert abs(noise_out - noise_in) < 0.2, \
        f"non-RFI beams shifted: {noise_in:.2f} -> {noise_out:.2f}"

    print("[test] PASS (end-to-end 2-bit roundtrip)")

    # --- Direct kernel test (no quantization) ---
    # Show the eigen-clean kernel exactly suppresses a planted rank-1 RFI
    # signal across >= min_support beams when fed whitened floats directly.
    B = 24; nF = 16; nT = 200; Fw = 4; Tw = 20
    nFt = nF // Fw; nTt = nT // Tw
    rng2 = np.random.default_rng(7)
    x = rng2.standard_normal((B, nF, nT)).astype(np.float32)
    rfi_mask_beams = np.array([1.0] * 20 + [0.0] * 4, dtype=np.float32)
    burst_t = np.arange(40, 60)
    x[:, :, burst_t] += (3.0 * rfi_mask_beams)[:, None, None]
    cp = MBC.CleanParams(Tw, Fw, min_support=15, safety=1.05, max_rank=3)
    cleaned, diag = MBC.clean_chunk(x, cp)
    print(f"[kernel] diag: tiles={diag.n_tiles} rfi_eigvals={diag.n_rfi_eigvals} "
          f"max_eig={diag.max_eigval:.2f} mp_edge={diag.mp_edge:.2f}")
    rfi_resid = float(np.abs(cleaned[:20, :, burst_t]).mean())
    nonrfi_resid = float(np.abs(cleaned[20:, :, burst_t] - x[20:, :, burst_t]).mean())
    print(f"[kernel] |cleaned| at burst on RFI beams       = {rfi_resid:.3f} "
          f"(was {float(np.abs(x[:20, :, burst_t]).mean()):.3f}; expect ~thermal)")
    print(f"[kernel] |cleaned - input| on non-RFI beams    = {nonrfi_resid:.3f} "
          f"(expect ~0)")
    assert rfi_resid < 1.5, f"kernel did not suppress planted RFI: {rfi_resid}"
    # Small leakage onto non-RFI beams is expected: the leading eigenvector is
    # estimated from finite samples and bleeds a little into orthogonal beams.
    assert nonrfi_resid < 0.1, f"kernel disturbed non-RFI beams too much: {nonrfi_resid}"
    print("[kernel] PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
