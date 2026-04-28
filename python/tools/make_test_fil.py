"""Generate small synthetic SIGPROC filterbank files for testing rfidiag.

The header format matches what dsa110-sigproc/src/read_header.c parses:
HEADER_START / HEADER_END text-tag binary header followed by raw packed
data. We emit a minimal but complete header (telescope_id, machine_id,
data_type=1, nchans, nbits=2, nifs=1, fch1, foff, tstart, tsamp).

Modes (-m):
  uniform          -- iid uniform {0,1,2,3} per (sample, channel)
  saturated_chan   -- like uniform but one channel is forced to 3
  transient        -- like uniform plus one bright zero-DM time sample
  gaussian         -- quantize a Gaussian (per-channel) to {0,1,2,3} via
                      thresholds chosen to give ~uniform marginal hist.
                      Useful as a sanity-check for the sigma exceedance plot
                      since the zero-DM CLT-sums many channels.

The output files are tiny by design (default ~5 s, 256 chans).
"""

from __future__ import annotations

import argparse
import os
import struct
import sys
from pathlib import Path

import numpy as np


def _send_str(buf: bytearray, s: str) -> None:
    b = s.encode("ascii")
    if not (1 <= len(b) <= 80):
        raise ValueError(f"sigproc string too long: {s!r}")
    buf += struct.pack("<i", len(b))
    buf += b


def _send_int(buf: bytearray, name: str, value: int) -> None:
    _send_str(buf, name)
    buf += struct.pack("<i", int(value))


def _send_double(buf: bytearray, name: str, value: float) -> None:
    _send_str(buf, name)
    buf += struct.pack("<d", float(value))


def build_header(*, telescope_id: int, machine_id: int, data_type: int,
                 nchans: int, nbits: int, nifs: int,
                 fch1: float, foff: float,
                 tstart: float, tsamp: float,
                 source_name: str = "TEST",
                 rawdatafile: str = "test.fil") -> bytes:
    buf = bytearray()
    _send_str(buf, "HEADER_START")
    _send_str(buf, "rawdatafile"); _send_str(buf, rawdatafile)
    _send_str(buf, "source_name"); _send_str(buf, source_name)
    _send_int(buf, "machine_id", machine_id)
    _send_int(buf, "telescope_id", telescope_id)
    _send_int(buf, "data_type", data_type)
    _send_int(buf, "nchans", nchans)
    _send_int(buf, "nbits", nbits)
    _send_int(buf, "nifs", nifs)
    _send_int(buf, "nbeams", 1)
    _send_int(buf, "ibeam", 1)
    _send_double(buf, "fch1", fch1)
    _send_double(buf, "foff", foff)
    _send_double(buf, "tstart", tstart)
    _send_double(buf, "tsamp", tsamp)
    _send_str(buf, "HEADER_END")
    return bytes(buf)


def pack_2bit(samples: np.ndarray) -> np.ndarray:
    """Pack uint8 samples (values 0..3) into bytes, sigproc convention.

    Sample order matches dsa110-sigproc/src/pack_unpack.c:
        byte = s0 | (s1<<2) | (s2<<4) | (s3<<6)
    The flat order of `samples` is preserved; length must be a multiple of 4.
    """
    a = np.ascontiguousarray(samples, dtype=np.uint8)
    if a.size % 4 != 0:
        raise ValueError("pack_2bit: input length must be a multiple of 4")
    a = a.reshape(-1, 4)
    if (a > 3).any():
        raise ValueError("pack_2bit: values must be in {0,1,2,3}")
    out = (a[:, 0] | (a[:, 1] << 2) | (a[:, 2] << 4) | (a[:, 3] << 6)).astype(np.uint8)
    return out


def make_samples(mode: str, *, nsamp: int, nchans: int, rng: np.random.Generator) -> np.ndarray:
    if (nsamp * nchans) % 4 != 0:
        raise ValueError("nsamp*nchans must be a multiple of 4 for 2-bit packing")

    if mode == "uniform":
        return rng.integers(0, 4, size=(nsamp, nchans), dtype=np.uint8)

    if mode == "saturated_chan":
        x = rng.integers(0, 4, size=(nsamp, nchans), dtype=np.uint8)
        x[:, nchans // 3] = 3
        return x

    if mode == "transient":
        x = rng.integers(0, 4, size=(nsamp, nchans), dtype=np.uint8)
        # one bright sample: force every channel to 3 at t=nsamp//2.
        x[nsamp // 2, :] = 3
        return x

    if mode == "gaussian":
        # Quantize a unit-variance Gaussian per (sample, chan) to {0,1,2,3}
        # using thresholds picked so that the four output bins are equiprobable
        # on a Gaussian. Values are the inverse-normal CDF at p=0.25, 0.5, 0.75.
        edges = np.array([-0.6744897501960817, 0.0, 0.6744897501960817])
        g = rng.standard_normal(size=(nsamp, nchans))
        x = np.digitize(g, edges).astype(np.uint8)
        return x

    raise ValueError(f"unknown mode: {mode}")


def write_fil(path: Path, samples: np.ndarray, *,
              fch1: float, foff: float, tsamp: float, tstart: float,
              source_name: str = "TEST") -> None:
    nsamp, nchans = samples.shape
    hdr = build_header(
        telescope_id=0, machine_id=0, data_type=1,
        nchans=nchans, nbits=2, nifs=1,
        fch1=fch1, foff=foff,
        tstart=tstart, tsamp=tsamp,
        source_name=source_name,
        rawdatafile=path.name,
    )
    packed = pack_2bit(samples.reshape(-1))
    with open(path, "wb") as f:
        f.write(hdr)
        f.write(packed.tobytes())


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", help="path to output .fil file")
    p.add_argument("-m", "--mode", default="uniform",
                   choices=["uniform", "saturated_chan", "transient", "gaussian"])
    p.add_argument("--nchans", type=int, default=256)
    p.add_argument("--nsamp",  type=int, default=20000,
                   help="number of time samples (must give nsamp*nchans %% 4 == 0)")
    p.add_argument("--tsamp",  type=float, default=50e-6, help="sample interval, seconds")
    p.add_argument("--fch1",   type=float, default=1850.0, help="MHz, top of band")
    p.add_argument("--foff",   type=float, default=-(1850.0 - 1300.0) / 2048.0,
                   help="channel width MHz (negative for descending)")
    p.add_argument("--tstart", type=float, default=60000.0, help="MJD")
    p.add_argument("--seed",   type=int, default=0)
    args = p.parse_args(argv)

    rng = np.random.default_rng(args.seed)
    samples = make_samples(args.mode, nsamp=args.nsamp, nchans=args.nchans, rng=rng)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    write_fil(out_path, samples,
              fch1=args.fch1, foff=args.foff, tsamp=args.tsamp,
              tstart=args.tstart, source_name=f"SYNTH_{args.mode}")
    nbytes = out_path.stat().st_size
    print(f"wrote {out_path}  ({nbytes/1024:.1f} KiB, {args.nsamp} samp x {args.nchans} chan, mode={args.mode})")


if __name__ == "__main__":
    main()
