"""Generate PNG diagnostic plots from a rfidiag .diag file.

Usage:
    python plot_diagnostics.py <out.diag> [--outdir plots/]
                               [--fdec FDEC] [--tdec TDEC]
                               [--n-spectra 5] [--spectra-sec 20]

Produces:
    chan_summary.png            mean/min/max/RMS spectrum + 2-bit clipping fractions
    waterfall.png               downsampled time-frequency waterfall of mean power
    bandpass_stability.png      per-channel RMS-vs-time map (chunk_rms / chan-mean)
    zerodm_05s.png              total power vs time at 0.5 s resolution
    zerodm_powerspectra.png     5x 20-s power-spectra (log-log) of zero-DM
    zerodm_sigma_exceedance.png exceedance fraction vs N-sigma compared to Gaussian
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path
from typing import Optional

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

REPO = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO))
from diagio import DiagFile, read_diag  # noqa: E402


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def _safe_log(x: np.ndarray, floor: float | None = None) -> np.ndarray:
    """log of positive values; floor non-positive entries to avoid -inf."""
    out = np.array(x, dtype=np.float64, copy=True)
    if floor is None:
        pos = out[out > 0]
        floor = pos.min() * 0.1 if pos.size else 1e-30
    out[out <= 0] = floor
    return np.log10(out)


def _decimate_axis(arr: np.ndarray, axis: int, factor: int, how: str = "mean") -> np.ndarray:
    if factor <= 1:
        return arr
    n = arr.shape[axis]
    nused = (n // factor) * factor
    sl = [slice(None)] * arr.ndim
    sl[axis] = slice(0, nused)
    a = arr[tuple(sl)]
    new_shape = list(a.shape)
    new_shape[axis] = nused // factor
    new_shape.insert(axis + 1, factor)
    a = a.reshape(new_shape)
    if how == "mean":
        return a.mean(axis=axis + 1)
    if how == "max":
        return a.max(axis=axis + 1)
    raise ValueError(how)


def _per_chan_stats(d: DiagFile):
    n = d.nsamples
    sumc = d.sumc.astype(np.float64)
    sumq = d.sumq.astype(np.float64)
    mean = sumc / n
    var = np.maximum(sumq / n - mean * mean, 0.0)
    rms = np.sqrt(var)
    cmean = d.cmean
    minspec = cmean.min(axis=0)
    maxspec = cmean.max(axis=0)
    return mean, rms, minspec, maxspec


def _apply_default_style():
    plt.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 130,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
    })


# --------------------------------------------------------------------------
# Plots
# --------------------------------------------------------------------------
def plot_chan_summary(d: DiagFile, out_path: Path):
    f = d.freqs_mhz
    mean, rms, minspec, maxspec = _per_chan_stats(d)
    n = d.nsamples
    if d.nbits == 2:
        clip0 = d.hist[:, 0] / n
        clip3 = d.hist[:, d.nhist - 1] / n
    else:
        clip0 = d.hist[:, 0] / n
        clip3 = d.hist[:, d.nhist - 1] / n

    fig, axes = plt.subplots(4, 1, figsize=(10, 9), sharex=True)
    ax = axes[0]
    ax.plot(f, mean, lw=0.7, color="C0")
    ax.set_ylabel("mean")
    ax.set_title(f"Per-channel diagnostics (n={n} samples, {n*d.tsamp:.2f} s)")
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.fill_between(f, minspec, maxspec, alpha=0.3, color="C2", label="min..max over chunks")
    ax.plot(f, mean, lw=0.7, color="C0", label="mean")
    ax.set_ylabel("min / max")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)

    ax = axes[2]
    ax.plot(f, rms, lw=0.7, color="C1")
    ax.set_ylabel("direct RMS")
    ax.grid(alpha=0.3)

    ax = axes[3]
    label_lo = "v=0 (low clip)"
    label_hi = f"v={d.nhist - 1} (high clip)"
    ax.semilogy(f, np.clip(clip0, 1e-9, None), lw=0.7, label=label_lo, color="C3")
    ax.semilogy(f, np.clip(clip3, 1e-9, None), lw=0.7, label=label_hi, color="C4")
    ax.axhline(1.0 / d.nhist, color="grey", ls="--", lw=0.6,
               label=f"uniform expectation = 1/{d.nhist}")
    ax.set_ylabel("clipping fraction")
    ax.set_xlabel("frequency (MHz)")
    ax.set_ylim(1e-6, 1.0)
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3, which="both")

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_waterfall(d: DiagFile, out_path: Path, fdec: int, tdec: int):
    """Time-frequency waterfall of chunk_mean, downsampled for tractability."""
    cmean = d.cmean
    wf = _decimate_axis(cmean, axis=1, factor=fdec, how="mean")
    wf = _decimate_axis(wf,    axis=0, factor=tdec, how="mean")
    nT, nF = wf.shape
    f = d.freqs_mhz
    f_centers = f[: (len(f) // fdec) * fdec].reshape(-1, fdec).mean(axis=1)
    t_centers = d.chunk_times_s[: (d.nchunks // tdec) * tdec].reshape(-1, tdec).mean(axis=1)

    fig, ax = plt.subplots(figsize=(10, 6))
    extent = (f_centers[0], f_centers[-1], t_centers[-1], t_centers[0])
    pos = wf.copy()
    pos[pos <= 0] = pos[pos > 0].min() * 0.1 if (pos > 0).any() else 1e-3
    im = ax.imshow(pos, aspect="auto", extent=extent,
                   norm=LogNorm(vmin=pos.min(), vmax=pos.max()),
                   cmap="viridis", interpolation="nearest")
    ax.set_xlabel("frequency (MHz)")
    ax.set_ylabel("time since data start (s)")
    ax.set_title(f"Waterfall (chunk_mean, fdec={fdec} -> {nF} chans, "
                 f"tdec={tdec} -> {nT} bins)")
    fig.colorbar(im, ax=ax, label="mean power per channel (log)")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_bandpass_stability(d: DiagFile, out_path: Path, fdec: int, tdec: int):
    crms = d.crms
    # Normalize each channel by its time-median RMS so we see relative drift.
    norm = np.median(crms, axis=0)
    norm = np.where(norm > 0, norm, 1.0)
    rel = crms / norm[None, :]
    rel = _decimate_axis(rel, axis=1, factor=fdec, how="mean")
    rel = _decimate_axis(rel, axis=0, factor=tdec, how="mean")
    nT, nF = rel.shape
    f = d.freqs_mhz
    f_centers = f[: (len(f) // fdec) * fdec].reshape(-1, fdec).mean(axis=1)
    t_centers = d.chunk_times_s[: (d.nchunks // tdec) * tdec].reshape(-1, tdec).mean(axis=1)

    # Symmetric color limits around 1.0 in log space, robust to outliers.
    finite = rel[np.isfinite(rel) & (rel > 0)]
    if finite.size:
        lo, hi = np.percentile(finite, [2, 98])
    else:
        lo, hi = 0.5, 2.0

    fig, ax = plt.subplots(figsize=(10, 6))
    extent = (f_centers[0], f_centers[-1], t_centers[-1], t_centers[0])
    im = ax.imshow(rel, aspect="auto", extent=extent,
                   vmin=lo, vmax=hi, cmap="RdBu_r", interpolation="nearest")
    ax.set_xlabel("frequency (MHz)")
    ax.set_ylabel("time since data start (s)")
    ax.set_title(f"Bandpass stability: chunk_rms / median(channel)  "
                 f"(fdec={fdec}, tdec={tdec})")
    fig.colorbar(im, ax=ax, label="ratio (1.0 = stable)")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def _zerodm_05s(d: DiagFile) -> tuple[np.ndarray, np.ndarray]:
    if d.zerodm is None:
        # Fall back to chunk-mean summed across channels
        z = d.cmean.sum(axis=1) * d.nsamples_per_chunk
        t = d.chunk_times_s
        return t, z
    z = d.zerodm
    block_n = max(1, int(round(0.5 / d.tsamp)))
    nblocks = z.size // block_n
    z = z[: nblocks * block_n].reshape(nblocks, block_n).mean(axis=1)
    t = d.t_offset + 0.5 * (np.arange(nblocks) + 0.5)
    return t, z


def plot_zerodm_05s(d: DiagFile, out_path: Path):
    t, z = _zerodm_05s(d)
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(t, z, lw=0.7)
    ax.set_xlabel("time since data start (s)")
    ax.set_ylabel("mean total power (sum over chans)")
    ax.set_title("Zero-DM time series, 0.5 s blocks")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_zerodm_powerspectra(d: DiagFile, out_path: Path,
                             n_spectra: int, slice_sec: float):
    if d.zerodm is None:
        # Without full-rate zero-DM we cannot resolve fast power spectra.
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5,
                "no full-rate zero-DM in this .diag\n(rerun without -Z)",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        fig.savefig(out_path)
        plt.close(fig)
        return

    z = d.zerodm.astype(np.float64)
    z = z - z.mean()
    nslice = int(round(slice_sec / d.tsamp))
    if nslice < 8:
        raise ValueError("slice too short for FFT")
    if nslice * n_spectra > z.size:
        # shrink slice or n
        n_spectra = max(1, z.size // nslice)
    starts = np.linspace(0, z.size - nslice, n_spectra, dtype=np.int64)

    freqs = np.fft.rfftfreq(nslice, d=d.tsamp)
    # drop DC bin for log-log
    keep = slice(1, None)

    fig, ax = plt.subplots(figsize=(10, 5))
    for i, s0 in enumerate(starts):
        seg = z[s0 : s0 + nslice]
        seg = seg - seg.mean()
        # boxcar window is fine; user requested simple FFT
        spec = np.abs(np.fft.rfft(seg)) ** 2 / nslice
        ax.loglog(freqs[keep], spec[keep], lw=0.6,
                  label=f"t = {(d.t_offset + s0*d.tsamp):.1f} s")
    ax.set_xlabel("frequency (Hz)")
    ax.set_ylabel("|FFT|^2 / N  (arbitrary units)")
    ax.set_title(f"Zero-DM power spectra of {n_spectra} x {slice_sec:.1f}-s slices")
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="upper right", ncol=2)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_sigma_exceedance(d: DiagFile, out_path: Path):
    if d.zerodm is None:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5,
                "no full-rate zero-DM in this .diag\n(rerun without -Z)",
                ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        fig.savefig(out_path)
        plt.close(fig)
        return

    z = d.zerodm.astype(np.float64)
    mu = z.mean()
    sigma = z.std()
    if sigma <= 0:
        sigma = 1.0
    s = (z - mu) / sigma
    n = s.size

    ks = np.arange(3.0, 10.5, 0.5)
    pos = np.array([(s > k).sum() / n for k in ks])
    neg = np.array([(s < -k).sum() / n for k in ks])
    # |s| > k => 2-sided
    both = np.array([(np.abs(s) > k).sum() / n for k in ks])

    # Gaussian expectation: P(|X|>k) = erfc(k/sqrt(2)),  P(X>k) = 0.5*erfc(k/sqrt(2))
    from math import erfc, sqrt
    gauss_two = np.array([erfc(k / sqrt(2)) for k in ks])
    gauss_one = 0.5 * gauss_two

    floor = 0.5 / n  # smallest measurable rate

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(ks, np.clip(pos,  floor, None), "o-", color="C3", label="positive tail")
    ax.semilogy(ks, np.clip(neg,  floor, None), "s-", color="C0", label="negative tail")
    ax.semilogy(ks, np.clip(both, floor, None), "^-", color="C2", label="|x| tail")
    ax.semilogy(ks, gauss_one, "--", color="grey", label="Gaussian one-sided")
    ax.semilogy(ks, gauss_two, ":",  color="grey", label="Gaussian two-sided")
    ax.set_xlabel("threshold k (sigmas, rms-based)")
    ax.set_ylabel("fraction of zero-DM samples")
    ax.set_title(f"Zero-DM sigma exceedance vs Gaussian "
                 f"(mu={mu:.3g}, rms={sigma:.3g}, n={n})")
    ax.legend(loc="lower left")
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------
def main(argv: Optional[list[str]] = None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("diag", help="path to a .diag file produced by rfidiag")
    p.add_argument("--outdir", default=None,
                   help="output directory (default: alongside .diag)")
    p.add_argument("--fdec", type=int, default=None,
                   help="freq decimation for waterfall/stability plots "
                        "(default: chosen so output has <=512 channels)")
    p.add_argument("--tdec", type=int, default=1,
                   help="time decimation for waterfall/stability plots (default: 1)")
    p.add_argument("--n-spectra", type=int, default=5,
                   help="number of zero-DM power-spectrum slices")
    p.add_argument("--spectra-sec", type=float, default=20.0,
                   help="duration of each zero-DM power-spectrum slice (s)")
    args = p.parse_args(argv)

    diag_path = Path(args.diag)
    outdir = Path(args.outdir) if args.outdir else diag_path.parent / (diag_path.stem + "_plots")
    outdir.mkdir(parents=True, exist_ok=True)

    d = read_diag(diag_path)

    if args.fdec is None:
        args.fdec = max(1, math.ceil(d.nchans / 512))

    _apply_default_style()
    plot_chan_summary       (d, outdir / "chan_summary.png")
    plot_waterfall          (d, outdir / "waterfall.png",
                             fdec=args.fdec, tdec=args.tdec)
    plot_bandpass_stability (d, outdir / "bandpass_stability.png",
                             fdec=args.fdec, tdec=args.tdec)
    plot_zerodm_05s         (d, outdir / "zerodm_05s.png")
    plot_zerodm_powerspectra(d, outdir / "zerodm_powerspectra.png",
                             n_spectra=args.n_spectra, slice_sec=args.spectra_sec)
    plot_sigma_exceedance   (d, outdir / "zerodm_sigma_exceedance.png")

    print(f"Wrote 6 plots to {outdir}")


if __name__ == "__main__":
    main()
