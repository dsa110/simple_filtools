#!/usr/bin/env python3
"""
Compare three SIGPROC filterbank files (raw, cleaned, filtool-style flagged).

If --flagged is not provided, this script generates a filtool-style flagged file
from the raw input using a lightweight approximation of:
    rfi_flag "kadaneF 4 8 zdot"

Outputs:
  - zero_dm_overlay.png              (overplotted time series)
  - avg_spectrum_overlay.png         (overplotted time-averaged bandpass)
  - waterfall_4x16.png               (3-panel waterfall, freq x4 / time x16)
  - zero_dm_powerspectra_overlay.png (overplotted PSD)
  - summary.json                     (basic run metadata)
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from multibeam_clean import (
    SigprocHeader,
    open_input,
    write_sigproc_header,
    unpack_bits_to_float32,
    pack_float32_to_bits,
)


@dataclass
class Stats:
    tsamp: float
    freqs_mhz: np.ndarray
    nchans: int
    sum_spec: np.ndarray
    n_samp: int = 0
    zdm_chunks: list = None
    wf_chunks: list = None

    def __post_init__(self):
        self.zdm_chunks = []
        self.wf_chunks = []

    def add_chunk(self, x: np.ndarray, max_wf_time_bins: int):
        # x: (T, C), float32
        self.sum_spec += x.sum(axis=0, dtype=np.float64)
        self.n_samp += x.shape[0]
        self.zdm_chunks.append(x.mean(axis=1, dtype=np.float64))

        c4 = (x.shape[1] // 4) * 4
        t16 = (x.shape[0] // 16) * 16
        if c4 == 0 or t16 == 0:
            return
        wf = x[:t16, :c4].reshape(t16, c4 // 4, 4).mean(axis=2)
        wf = wf.reshape(t16 // 16, 16, c4 // 4).mean(axis=1)  # (T/16, C/4)
        self.wf_chunks.append(wf.astype(np.float16, copy=False))

        # Soft memory guard: if we grow too large in time bins, coarsen by x2.
        n_bins = sum(w.shape[0] for w in self.wf_chunks)
        if n_bins > max_wf_time_bins * 2:
            merged = np.vstack(self.wf_chunks).astype(np.float32, copy=False)
            n2 = (merged.shape[0] // 2) * 2
            merged = merged[:n2].reshape(-1, 2, merged.shape[1]).mean(axis=1)
            self.wf_chunks = [merged.astype(np.float16)]

    def finalize(self, max_wf_time_bins: int) -> Dict[str, np.ndarray]:
        zdm = np.concatenate(self.zdm_chunks) if self.zdm_chunks else np.zeros(0)
        spec = self.sum_spec / max(self.n_samp, 1)
        wf = np.vstack(self.wf_chunks).astype(np.float32) if self.wf_chunks else np.zeros((0, self.nchans // 4))
        if wf.shape[0] > max_wf_time_bins:
            f = int(np.ceil(wf.shape[0] / max_wf_time_bins))
            n = (wf.shape[0] // f) * f
            wf = wf[:n].reshape(-1, f, wf.shape[1]).mean(axis=1)
            wf_dt = self.tsamp * 16 * f
        else:
            wf_dt = self.tsamp * 16
        return {"zdm": zdm, "spec": spec, "wf": wf, "wf_dt": wf_dt}


def _robust_z(x: np.ndarray) -> np.ndarray:
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    s = 1.4826 * mad if mad > 0 else (x.std() + 1e-6)
    return (x - med) / max(s, 1e-6)


def _kadane_positive_segments(a: np.ndarray):
    """Return positive-sum max-subarray segments iteratively."""
    x = a.copy()
    segs = []
    for _ in range(64):
        best = 0.0
        cur = 0.0
        s = 0
        bs = -1
        be = -1
        for i, v in enumerate(x):
            if cur <= 0:
                cur = v
                s = i
            else:
                cur += v
            if cur > best:
                best = cur
                bs = s
                be = i
        if best <= 0 or bs < 0:
            break
        segs.append((bs, be + 1, best))
        x[bs:be + 1] = 0.0
    return segs


def filtool_style_kadanef_zdot(x: np.ndarray, kadane_widths=(4, 8), snr_thre=7.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Approximation of filtool's "kadaneF 4 8 zdot":
      1) kadaneF-like contiguous frequency masking from channel variance excess
      2) zdot-like linear trend removal in zero-DM time series
    Returns (x_flagged, channel_mask_bool).
    """
    y = x.astype(np.float32, copy=True)
    t = np.arange(y.shape[0], dtype=np.float32)

    # --- kadaneF-like channel score
    ch_mu = y.mean(axis=0)
    ch_sd = y.std(axis=0)
    ch_sd[ch_sd <= 1e-3] = 1e-3
    z = (y - ch_mu[None, :]) / ch_sd[None, :]
    score = np.mean(z * z, axis=0) - 1.0
    score = _robust_z(score)

    mask = np.zeros(y.shape[1], dtype=bool)
    for w in kadane_widths:
        w = max(1, int(w))
        ker = np.ones(w, dtype=np.float32) / float(w)
        smooth = np.convolve(score, ker, mode="same")
        arr = smooth - float(snr_thre)
        for a, b, _ in _kadane_positive_segments(arr):
            if (b - a) >= w:
                mask[a:b] = True

    if mask.any() and (~mask).any():
        fill = np.median(y[:, ~mask], axis=1).astype(np.float32, copy=False)
        y[:, mask] = fill[:, None]

    # --- zdot-like zero-DM linear detrend
    zdm = y.mean(axis=1)
    t0 = t - t.mean()
    denom = float((t0 * t0).sum()) + 1e-9
    slope = float((t0 * (zdm - zdm.mean())).sum()) / denom
    trend = slope * t0 + zdm.mean()
    y -= trend[:, None]
    y += float(np.median(x))  # keep quantization range anchored

    return y, mask


def process_file(path: Path, chunk_samps: int, max_wf_time_bins: int, apply_flagger=False,
                 flagged_out: Optional[Path] = None,
                 kadane_widths=(4, 8), snr_thre=7.0):
    bf = open_input(path)
    hdr = bf.header
    nchans = int(hdr.get("nchans"))
    nbits = int(hdr.get("nbits"))
    tsamp = float(hdr.get("tsamp"))
    freqs = float(hdr.get("fch1")) + np.arange(nchans, dtype=np.float64) * float(hdr.get("foff"))

    out_fp = None
    if apply_flagger and flagged_out is not None:
        flagged_out.parent.mkdir(parents=True, exist_ok=True)
        out_fp = open(flagged_out, "wb")
        out_hdr = SigprocHeader(fields=list(hdr.fields))
        if out_hdr.get("rawdatafile") is not None:
            out_hdr.set("rawdatafile", flagged_out.name)
        write_sigproc_header(out_fp, out_hdr)

    st = Stats(tsamp=tsamp, freqs_mhz=freqs, nchans=nchans,
               sum_spec=np.zeros(nchans, dtype=np.float64))
    nflag = 0
    nseen = 0

    bps = (nbits * nchans) // 8
    remain = bf.nsamples_total
    while remain > 0:
        t = min(chunk_samps, remain)
        flat = t * nchans
        nbytes = (flat * nbits + 7) // 8
        packed = np.frombuffer(bf.fp_in.read(nbytes), dtype=np.uint8)
        x = unpack_bits_to_float32(packed, nbits, flat).reshape(t, nchans)

        if apply_flagger:
            y, m = filtool_style_kadanef_zdot(x, kadane_widths=kadane_widths, snr_thre=snr_thre)
            nflag += int(m.sum())
            nseen += m.size
            if out_fp is not None:
                packed_out = pack_float32_to_bits(y.reshape(-1), nbits)
                out_fp.write(packed_out.tobytes())
            st.add_chunk(y, max_wf_time_bins=max_wf_time_bins)
        else:
            st.add_chunk(x, max_wf_time_bins=max_wf_time_bins)

        remain -= t

    bf.fp_in.close()
    if out_fp is not None:
        out_fp.close()

    out = st.finalize(max_wf_time_bins=max_wf_time_bins)
    out["tsamp"] = tsamp
    out["freqs"] = freqs
    out["flag_frac"] = (nflag / max(nseen, 1)) if apply_flagger else 0.0
    return out


def welch_like(zdm: np.ndarray, tsamp: float, nperseg: int = 65536) -> Tuple[np.ndarray, np.ndarray]:
    if zdm.size < 16:
        return np.zeros(0), np.zeros(0)
    z = zdm.astype(np.float64) - float(np.mean(zdm))
    nper = min(nperseg, z.size)
    step = nper // 2
    win = np.hanning(nper)
    norm = (win * win).sum()
    acc = None
    nseg = 0
    for i in range(0, z.size - nper + 1, max(step, 1)):
        seg = z[i:i + nper] * win
        p = np.abs(np.fft.rfft(seg)) ** 2 / max(norm, 1.0)
        acc = p if acc is None else (acc + p)
        nseg += 1
    if nseg == 0:
        p = np.abs(np.fft.rfft(z * np.hanning(z.size))) ** 2
        f = np.fft.rfftfreq(z.size, d=tsamp)
        return f, p
    p = acc / nseg
    f = np.fft.rfftfreq(nper, d=tsamp)
    return f, p


def _plot(args, raw, clean, flag):
    args.outdir.mkdir(parents=True, exist_ok=True)
    labels = ["raw", "cleaned", "flagged(kadaneF 4 8 zdot-style)"]
    data = [raw, clean, flag]
    colors = ["C0", "C2", "C3"]

    # 1) zero-DM timeseries
    fig, ax = plt.subplots(figsize=(12, 4))
    for lab, d, c in zip(labels, data, colors):
        z = d["zdm"]
        if z.size == 0:
            continue
        stride = max(1, z.size // 200000)
        t = np.arange(0, z.size, stride) * d["tsamp"]
        ax.plot(t, z[::stride], lw=0.7, alpha=0.85, color=c, label=lab)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("zero-DM intensity")
    ax.set_title("Zero-DM time series")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.savefig(args.outdir / "zero_dm_overlay.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # 2) time-averaged spectrum
    fig, ax = plt.subplots(figsize=(12, 4))
    f = raw["freqs"]
    for lab, d, c in zip(labels, data, colors):
        ax.plot(f, d["spec"], lw=0.9, alpha=0.9, color=c, label=lab)
    ax.set_xlabel("frequency (MHz)")
    ax.set_ylabel("mean intensity")
    ax.set_title("Time-averaged spectrum")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.savefig(args.outdir / "avg_spectrum_overlay.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # 3) waterfall binned by 4x freq and 16x time
    fig, axs = plt.subplots(3, 1, figsize=(12, 8), sharex=True, sharey=True)
    f4 = f[: (len(f) // 4) * 4].reshape(-1, 4).mean(axis=1)
    for ax, lab, d in zip(axs, labels, data):
        wf = d["wf"]
        if wf.size == 0:
            ax.set_title(f"{lab}: no waterfall data")
            continue
        tmax = wf.shape[0] * d["wf_dt"]
        v = np.log10(np.maximum(wf, 1e-3))
        im = ax.imshow(v.T, aspect="auto", origin="lower",
                       extent=[0, tmax, f4.min(), f4.max()], cmap="magma")
        ax.set_title(lab)
        ax.set_ylabel("freq (MHz)")
        fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label="log10(I)")
    axs[-1].set_xlabel("time (s)")
    fig.suptitle("Waterfall (freq x4, time x16; plus adaptive extra time decimation if needed)")
    fig.savefig(args.outdir / "waterfall_4x16.png", dpi=140, bbox_inches="tight")
    plt.close(fig)

    # 4) zero-DM power spectra
    fig, ax = plt.subplots(figsize=(10, 4))
    for lab, d, c in zip(labels, data, colors):
        ff, pp = welch_like(d["zdm"], d["tsamp"])
        if ff.size == 0:
            continue
        ax.loglog(ff[1:], np.maximum(pp[1:], 1e-20), lw=0.9, alpha=0.9, color=c, label=lab)
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("PSD (arb.)")
    ax.set_title("Zero-DM power spectra")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    fig.savefig(args.outdir / "zero_dm_powerspectra_overlay.png", dpi=140, bbox_inches="tight")
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("raw_fil", type=Path)
    p.add_argument("cleaned_fil", type=Path)
    p.add_argument("--flagged", type=Path, default=None,
                   help="Existing flagged .fil to compare; if absent, generate from raw.")
    p.add_argument("--flagged-out", type=Path, default=None,
                   help="Where to write generated flagged file (default: <outdir>/<rawstem>_kadaneF48zdot.fil)")
    p.add_argument("--outdir", type=Path, default=Path("triplet_compare"))
    p.add_argument("--chunk-samps", type=int, default=32768)
    p.add_argument("--kadane-widths", default="4,8", help="Comma list, e.g. 4,8")
    p.add_argument("--kadane-snr", type=float, default=7.0)
    p.add_argument("--max-waterfall-time-bins", type=int, default=4096)
    args = p.parse_args(argv)

    if args.flagged is None:
        flagged_out = args.flagged_out or (args.outdir / f"{args.raw_fil.stem}_kadaneF48zdot.fil")
    else:
        flagged_out = args.flagged

    widths = tuple(int(w.strip()) for w in args.kadane_widths.split(",") if w.strip())
    if len(widths) == 0:
        raise ValueError("No --kadane-widths provided")

    print(f"[1/3] processing raw: {args.raw_fil}")
    raw = process_file(args.raw_fil, args.chunk_samps, args.max_waterfall_time_bins,
                       apply_flagger=False)
    print(f"[2/3] processing cleaned: {args.cleaned_fil}")
    clean = process_file(args.cleaned_fil, args.chunk_samps, args.max_waterfall_time_bins,
                         apply_flagger=False)

    if args.flagged is None:
        print(f"[3/3] generating + processing flagged: {flagged_out}")
        flag = process_file(args.raw_fil, args.chunk_samps, args.max_waterfall_time_bins,
                            apply_flagger=True, flagged_out=flagged_out,
                            kadane_widths=widths, snr_thre=args.kadane_snr)
    else:
        print(f"[3/3] processing existing flagged: {flagged_out}")
        flag = process_file(flagged_out, args.chunk_samps, args.max_waterfall_time_bins,
                            apply_flagger=False)

    # compatibility check
    if len(raw["spec"]) != len(clean["spec"]) or len(raw["spec"]) != len(flag["spec"]):
        raise ValueError("raw/cleaned/flagged have incompatible channel counts")

    _plot(args, raw, clean, flag)

    meta = {
        "raw_fil": str(args.raw_fil),
        "cleaned_fil": str(args.cleaned_fil),
        "flagged_fil": str(flagged_out),
        "generated_flagged": args.flagged is None,
        "kadane_widths": list(widths),
        "kadane_snr": args.kadane_snr,
        "raw_nsamp": int(raw["zdm"].size),
        "clean_nsamp": int(clean["zdm"].size),
        "flag_nsamp": int(flag["zdm"].size),
        "flag_channel_fraction_est": float(flag.get("flag_frac", 0.0)),
    }
    args.outdir.mkdir(parents=True, exist_ok=True)
    with open(args.outdir / "summary.json", "w") as f:
        json.dump(meta, f, indent=2)
    print(f"[done] wrote plots + summary to {args.outdir}")


if __name__ == "__main__":
    main()

