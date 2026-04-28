"""End-to-end validation: run rfidiag on a synthetic fil, then independently
recompute the same statistics in numpy and assert they agree.

Integer accumulators (sum, sumsq, hist) must be bit-exact.
Float arrays (cmean, crms, zerodm) must agree to <= 1e-5 relative.
"""

from __future__ import annotations

import os
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "python"))
from diagio import read_diag  # noqa: E402

TMP = REPO / "tests" / "_tmp"
RFIDIAG = REPO / "rfidiag"
GEN = REPO / "python" / "tools" / "make_test_fil.py"


def parse_fil_header(path: Path) -> tuple[dict, int]:
    """Tiny sigproc header parser, just enough for this test."""
    with open(path, "rb") as f:
        raw = f.read()
    pos = 0

    def get_str():
        nonlocal pos
        (n,) = struct.unpack_from("<i", raw, pos); pos += 4
        s = raw[pos:pos+n].decode("ascii"); pos += n
        return s

    tag = get_str()
    assert tag == "HEADER_START", tag
    hdr: dict = {"nifs": 1}
    expecting = None
    while True:
        tag = get_str()
        if tag == "HEADER_END":
            break
        if expecting:
            hdr[expecting] = tag
            expecting = None
            continue
        if tag in ("rawdatafile", "source_name"):
            expecting = tag
            continue
        if tag in ("machine_id", "telescope_id", "data_type", "nchans",
                   "nbits", "nifs", "nbeams", "ibeam"):
            (v,) = struct.unpack_from("<i", raw, pos); pos += 4
            hdr[tag] = int(v)
        elif tag in ("fch1", "foff", "tstart", "tsamp",
                     "src_raj", "src_dej", "az_start", "za_start", "refdm"):
            (v,) = struct.unpack_from("<d", raw, pos); pos += 8
            hdr[tag] = float(v)
        else:
            raise RuntimeError(f"unknown header tag: {tag}")
    return hdr, pos


def numpy_diag(fil_path: Path, chunk_sec: float):
    hdr, hbytes = parse_fil_header(fil_path)
    nchans = hdr["nchans"]
    nbits = hdr["nbits"]
    tsamp = hdr["tsamp"]
    assert nbits == 2, "this checker only handles nbits=2"

    raw = np.fromfile(fil_path, dtype=np.uint8, offset=hbytes)
    # Each byte holds 4 samples in sigproc 2-bit order: bits 0-1 = s0, ..., 6-7 = s3.
    s0 =  raw        & 0x3
    s1 = (raw >> 2)  & 0x3
    s2 = (raw >> 4)  & 0x3
    s3 = (raw >> 6)  & 0x3
    samples = np.stack([s0, s1, s2, s3], axis=-1).reshape(-1).astype(np.uint8)
    nsamp = samples.size // nchans
    samples = samples[: nsamp * nchans].reshape(nsamp, nchans)

    sumc = samples.astype(np.uint64).sum(axis=0)
    sumq = (samples.astype(np.uint64) ** 2).sum(axis=0)

    nhist = 1 << nbits
    hist = np.zeros((nchans, nhist), dtype=np.uint64)
    for v in range(nhist):
        hist[:, v] = (samples == v).sum(axis=0)

    chunk_n = int(round(chunk_sec / tsamp))
    # match the C-side byte-boundary rounding for nbits<8
    spb = 8 // nbits
    g = np.gcd(nchans, spb)
    m = spb // g
    if chunk_n % m != 0:
        chunk_n = ((chunk_n // m) + 1) * m
    nchunks = (nsamp + chunk_n - 1) // chunk_n

    cmean = np.zeros((nchunks, nchans), dtype=np.float32)
    crms  = np.zeros((nchunks, nchans), dtype=np.float32)
    for k in range(nchunks):
        s = samples[k*chunk_n : (k+1)*chunk_n].astype(np.float64)
        m = s.mean(axis=0)
        v = (s ** 2).mean(axis=0) - m * m
        v = np.clip(v, 0.0, None)
        cmean[k] = m.astype(np.float32)
        crms [k] = np.sqrt(v).astype(np.float32)

    zerodm = samples.astype(np.float32).sum(axis=1)
    return dict(
        nchans=nchans, nbits=nbits, nsamp=nsamp, chunk_n=chunk_n, nchunks=nchunks,
        sumc=sumc, sumq=sumq, hist=hist,
        cmean=cmean, crms=crms, zerodm=zerodm,
    )


def run(cmd: list[str]):
    print("$", " ".join(str(c) for c in cmd))
    r = subprocess.run(cmd, check=True, capture_output=True, text=True)
    if r.stdout: print(r.stdout)
    if r.stderr: print(r.stderr)


def check_one(fil: Path, chunk_sec: float = 0.5):
    diag_path = fil.with_suffix(".diag")
    run([str(RFIDIAG), str(fil), "-o", str(diag_path), "-c", str(chunk_sec), "-q"])

    d = read_diag(diag_path)
    ref = numpy_diag(fil, chunk_sec=chunk_sec)

    # Header-derived sizes
    assert d.nchans == ref["nchans"], (d.nchans, ref["nchans"])
    assert d.nbits == ref["nbits"]
    assert d.nsamples == ref["nsamp"], (d.nsamples, ref["nsamp"])
    assert d.nchunks == ref["nchunks"], (d.nchunks, ref["nchunks"])

    # Bit-exact integer comparisons
    assert np.array_equal(d.sumc, ref["sumc"]), "SUMC mismatch"
    assert np.array_equal(d.sumq, ref["sumq"]), "SUMQ mismatch"
    assert np.array_equal(d.hist, ref["hist"]), "HIST mismatch"

    # Float comparisons
    def close(a, b, name, atol=1e-5, rtol=1e-5):
        if not np.allclose(a, b, atol=atol, rtol=rtol):
            diff = np.abs(a - b)
            raise AssertionError(
                f"{name} mismatch: max |abs diff|={diff.max():.3g}, "
                f"max |rel diff|={(diff / (np.abs(b)+1e-30)).max():.3g}"
            )

    close(d.cmean, ref["cmean"], "CMEA")
    close(d.crms,  ref["crms"],  "CRMS")
    if d.zerodm is not None:
        close(d.zerodm, ref["zerodm"], "Z0DM")

    print(f"  OK: {fil.name}  nsamp={d.nsamples} nchunks={d.nchunks}")


def main():
    if not RFIDIAG.exists():
        print(f"rfidiag not built at {RFIDIAG}; run `make` first.", file=sys.stderr)
        sys.exit(2)

    TMP.mkdir(parents=True, exist_ok=True)
    cases = [
        ("uniform",        ["-m", "uniform",        "--nchans", "64",  "--nsamp", "20000"]),
        ("saturated_chan", ["-m", "saturated_chan", "--nchans", "64",  "--nsamp", "20000"]),
        ("transient",      ["-m", "transient",      "--nchans", "64",  "--nsamp", "20000"]),
        ("gaussian",       ["-m", "gaussian",       "--nchans", "256", "--nsamp", "200000"]),
    ]
    for name, args in cases:
        fil = TMP / f"check_{name}.fil"
        if fil.exists(): fil.unlink()
        run([sys.executable, str(GEN), str(fil), *args, "--seed", "1"])
        check_one(fil, chunk_sec=0.25)

    print("\nAll numpy checks passed.")


if __name__ == "__main__":
    main()
