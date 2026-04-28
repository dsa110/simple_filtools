"""Reader for rfidiag .diag binary files.

The on-disk format (see src/diag.h):

    [DiagHeader]                                  # fixed-size struct
    { tag[4], len_u64, payload[len] } * N         # tagged sections

Sections currently emitted:
    "SUMC"  uint64[nchans]               per-channel sample sum
    "SUMQ"  uint64[nchans]               per-channel sample sumsq
    "HIST"  uint64[nchans * nhist]       per-channel value histograms (nhist=2^nbits)
    "CMEA"  float32[nchunks * nchans]    per-chunk per-channel mean
    "CRMS"  float32[nchunks * nchans]    per-chunk per-channel RMS
    "Z0DM"  float32[nsamples]            full-rate zero-DM time series
"""

from __future__ import annotations

import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

import numpy as np


# Must match #pragma pack(1) struct DiagHeader in src/diag.h.
# Layout: 8s magic, I version, I nchans, I nbits, I nchunks, Q nsamples,
#         Q nsamples_per_chunk, d tsamp, d fch1, d foff, d tstart,
#         d chunk_sec, d t_offset, I nhist, I reserved
_HEADER_FMT = "<8sIIIIQQddddddII"
_HEADER_SIZE = struct.calcsize(_HEADER_FMT)


@dataclass
class DiagFile:
    # Header fields
    magic: bytes
    version: int
    nchans: int
    nbits: int
    nchunks: int
    nsamples: int
    nsamples_per_chunk: int
    tsamp: float
    fch1: float
    foff: float
    tstart: float
    chunk_sec: float
    t_offset: float
    nhist: int

    # Section payloads (numpy arrays, names match section tags)
    sections: Dict[str, np.ndarray]

    # Convenience accessors
    @property
    def nfreq(self) -> int:
        return self.nchans

    @property
    def freqs_mhz(self) -> np.ndarray:
        return self.fch1 + self.foff * np.arange(self.nchans, dtype=np.float64)

    @property
    def chunk_times_s(self) -> np.ndarray:
        return self.t_offset + self.chunk_sec * (np.arange(self.nchunks) + 0.5)

    @property
    def sample_times_s(self) -> np.ndarray:
        return self.t_offset + self.tsamp * np.arange(self.nsamples, dtype=np.float64)

    @property
    def sumc(self) -> np.ndarray:
        return self.sections["SUMC"]

    @property
    def sumq(self) -> np.ndarray:
        return self.sections["SUMQ"]

    @property
    def hist(self) -> np.ndarray:
        return self.sections["HIST"].reshape(self.nchans, self.nhist)

    @property
    def cmean(self) -> np.ndarray:
        return self.sections["CMEA"].reshape(self.nchunks, self.nchans)

    @property
    def crms(self) -> np.ndarray:
        return self.sections["CRMS"].reshape(self.nchunks, self.nchans)

    @property
    def zerodm(self) -> Optional[np.ndarray]:
        return self.sections.get("Z0DM")


def _read_exact(fp, n: int) -> bytes:
    b = fp.read(n)
    if len(b) != n:
        raise IOError(f"short read: wanted {n}, got {len(b)}")
    return b


def read_diag(path: str | Path) -> DiagFile:
    path = Path(path)
    with open(path, "rb") as fp:
        raw = _read_exact(fp, _HEADER_SIZE)
        (magic, version, nchans, nbits, nchunks, nsamples,
         nsamples_per_chunk, tsamp, fch1, foff, tstart,
         chunk_sec, t_offset, nhist, _reserved) = struct.unpack(_HEADER_FMT, raw)

        if not magic.startswith(b"FILDIAG"):
            raise ValueError(f"{path}: not a rfidiag .diag file (magic={magic!r})")
        if version != 1:
            raise ValueError(f"{path}: unsupported diag version {version}")

        sections: Dict[str, np.ndarray] = {}
        while True:
            head = fp.read(12)
            if not head:
                break
            if len(head) != 12:
                raise IOError(f"{path}: truncated section header")
            tag = head[:4].rstrip(b"\x00").decode("ascii", errors="replace")
            (length,) = struct.unpack("<Q", head[4:12])
            payload = _read_exact(fp, length)
            if tag in ("SUMC", "SUMQ", "HIST"):
                arr = np.frombuffer(payload, dtype=np.uint64).copy()
            elif tag in ("CMEA", "CRMS", "Z0DM"):
                arr = np.frombuffer(payload, dtype=np.float32).copy()
            else:
                # Unknown tag -> store as raw bytes so the file format stays
                # forward-compatible.
                arr = np.frombuffer(payload, dtype=np.uint8).copy()
            sections[tag] = arr

    return DiagFile(
        magic=magic,
        version=version,
        nchans=nchans,
        nbits=nbits,
        nchunks=nchunks,
        nsamples=nsamples,
        nsamples_per_chunk=nsamples_per_chunk,
        tsamp=tsamp,
        fch1=fch1,
        foff=foff,
        tstart=tstart,
        chunk_sec=chunk_sec,
        t_offset=t_offset,
        nhist=nhist,
        sections=sections,
    )
