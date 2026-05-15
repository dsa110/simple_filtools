"""Compare per-channel RFI flags across multiple .diag files (one per beam).

For each input .diag, derive a per-channel "bad" mask using a robust
z-score (median absolute deviation) on both the per-channel mean and
direct RMS. Channels whose mean OR rms is more than `--nsigma` MADs
from the across-channel median are flagged.

Aggregating across beams:
    n_beams_flagged[c] = number of beams that flag channel c
    common[c] = (n_beams_flagged[c] == nbeams)
    unique[c] = (n_beams_flagged[c] == 1)
    shared[c] = (1 < n_beams_flagged[c] < nbeams)
    clean[c]  = (n_beams_flagged[c] == 0)

If the .diag files contain zero-DM time series (Z0DM section, present
unless `rfidiag -Z` was used), this script also performs a parallel
analysis in the time domain: each beam's zero-DM is binned to
`--zdm-bin-sec` (default 1 ms), MAD-normalised, then compared across
beams via pairwise Pearson correlation and burst-coincidence at
`--zdm-burst-sigma` (default 6).

Usage
-----
    python compare_beam_rfi.py beam0.diag beam1.diag ... [options]

Spatial-filter feasibility (Kocz et al. 2010, AJ 140, 2086)
-----------------------------------------------------------
Kocz et al. demonstrate RFI removal on the Parkes Multibeam by SVD of the
*voltage* cross-spectrum covariance matrix C(f) = <S_i(f) S_j(f)*>; the
key requirement is that C(f) be dominated by a small number of strong
eigenvalues (the RFI subspace, paper Eq. 11). With detected total-power
filterbanks we cannot form that complex covariance directly, but the
real-valued cross-beam intensity covariance is a faithful proxy for
predicting whether voltage-domain Kocz would work on this RFI
environment: if RFI couples coherently across beams, a single eigenvalue
will dominate at each contaminated channel; if not, voltage-domain SVD
won't help either. This script computes:

  - whitening / gain spread across beams           (paper Section 3)
  - per-channel cross-beam intensity SVD            (analog of paper Fig. 2)
  - rank-k residual ("spatial-filter suppression") on the chunk-mean cube
  - PCA of the cross-beam zero-DM time series

and prints a verdict on whether the technique is likely to succeed if you
ever capture voltages.

Outputs (in --outdir, default ./compare_rfi/):
    flag_matrix.png            beams x channels heatmap of flag/no-flag
    n_beams_per_chan.png       n_beams_flagged vs frequency
    mean_rms_overlay.png       per-beam mean & rms spectra, common/unique bands
    summary_counts.png         #channels in each commonality bucket
    per_chan.csv               one row per channel
    zdm_overlay.png            per-beam zero-DM time series, MAD-normalised
    zdm_correlation.png        pairwise Pearson-correlation matrix of zero-DM
    zdm_sigma_exceedance.png   exceedance fraction vs k, per beam
    zdm_burst_coincidence.png  z-score and #beams-above-threshold vs time
    spatial_eigvals.png        per-channel cross-beam SVD spectrum (Kocz Fig. 2)
    spatial_eigratio.png       lambda_1 / lambda_med vs frequency (gain proxy)
    spatial_suppression.png    residual variance after rank-k subtraction
    spatial_zdm_pca.png        zero-DM PC time series + rank-k residual
    spatial_whitening.png      per-beam noise-floor spread (whitening check)

Options
-------
    --outdir PATH         where to write plots+csv  (default: ./compare_rfi/)
    --names a,b,c         explicit beam labels      (default: file stems)
    --nsigma N            MAD-z thresh for chan flag (default: 5)
    --no-mean-flag        don't use mean z-score
    --no-rms-flag         don't use rms z-score
    --top-k K             print K worst items / beam (default: 10)
    --zdm-bin-sec SEC     time bin for zero-DM cross-beam stats (default: 1e-3)
    --zdm-burst-sigma N   coincidence threshold (default: 6.0)
    --no-zdm              skip zero-DM analysis even if data is present
    --no-spatial          skip spatial-filter feasibility analysis
    --spatial-rank K      rank-k subtraction depth for suppression curve
                          (default: max(1, nbeams // 4))

Assumes all .diag files share nchans, fch1, foff, nbits. Raises if not.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO))
from diagio import DiagFile, read_diag  # noqa: E402


# --------------------------------------------------------------------------
# Per-beam analysis
# --------------------------------------------------------------------------
@dataclass
class BeamStats:
    name: str
    path: Path
    diag: DiagFile
    mean: np.ndarray   # (nchans,) per-channel mean intensity
    rms: np.ndarray    # (nchans,) per-channel direct RMS
    z_mean: np.ndarray # (nchans,) MAD-normalised z-score of mean
    z_rms: np.ndarray  # (nchans,) MAD-normalised z-score of rms
    flagged: np.ndarray  # (nchans,) bool


def _zmad(x: np.ndarray) -> np.ndarray:
    """Robust z-score: (x - median(x)) / (1.4826 * MAD(x))."""
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    if mad == 0.0:
        # degenerate: all values equal, no spread to measure
        return np.zeros_like(x)
    return (x - med) / (1.4826 * mad)


def analyse_beam(path: Path, name: str, nsigma: float,
                 use_mean: bool, use_rms: bool) -> BeamStats:
    d = read_diag(path)
    n = float(d.nsamples)
    if n <= 0:
        raise ValueError(f"{path}: empty diag (nsamples=0)")
    sumc = d.sumc.astype(np.float64)
    sumq = d.sumq.astype(np.float64)
    mean = sumc / n
    var = np.maximum(sumq / n - mean * mean, 0.0)
    rms = np.sqrt(var)
    z_mean = _zmad(mean)
    z_rms = _zmad(rms)
    flag = np.zeros(d.nchans, dtype=bool)
    if use_mean:
        flag |= np.abs(z_mean) > nsigma
    if use_rms:
        flag |= np.abs(z_rms) > nsigma
    if not (use_mean or use_rms):
        raise ValueError("must enable at least one of --use-mean / --use-rms")
    return BeamStats(name=name, path=path, diag=d,
                     mean=mean, rms=rms,
                     z_mean=z_mean, z_rms=z_rms,
                     flagged=flag)


def _check_compatible(beams: Sequence[BeamStats]) -> None:
    ref = beams[0].diag
    for b in beams[1:]:
        d = b.diag
        if (d.nchans, d.nbits) != (ref.nchans, ref.nbits):
            raise ValueError(
                f"{b.path}: incompatible (nchans={d.nchans} nbits={d.nbits}) "
                f"vs reference {beams[0].path} (nchans={ref.nchans} nbits={ref.nbits})")
        if not (np.isclose(d.fch1, ref.fch1) and np.isclose(d.foff, ref.foff)):
            raise ValueError(
                f"{b.path}: frequency grid (fch1={d.fch1}, foff={d.foff}) "
                f"does not match reference (fch1={ref.fch1}, foff={ref.foff})")


# --------------------------------------------------------------------------
# Helpers for printing / range listing
# --------------------------------------------------------------------------
def _runs(mask: np.ndarray) -> List[tuple]:
    """Return a list of (start, stop) for each contiguous True-run in `mask`,
    where stop is exclusive (numpy slice convention)."""
    if not mask.any():
        return []
    diff = np.diff(mask.astype(np.int8))
    starts = list(np.where(diff == 1)[0] + 1)
    stops  = list(np.where(diff == -1)[0] + 1)
    if mask[0]:
        starts.insert(0, 0)
    if mask[-1]:
        stops.append(len(mask))
    return list(zip(starts, stops))


def _format_runs(mask: np.ndarray, freqs: np.ndarray, top_k: int) -> List[str]:
    runs = _runs(mask)
    # sort descending by run length
    runs.sort(key=lambda ab: ab[1] - ab[0], reverse=True)
    out = []
    for (a, b) in runs[:top_k]:
        # freqs[a]..freqs[b-1]; foff may be negative, so use min/max for clarity
        f_lo = min(freqs[a], freqs[b - 1])
        f_hi = max(freqs[a], freqs[b - 1])
        out.append(f"  chans [{a:>5d}..{b-1:<5d}]  ({b-a:>4d} ch)  "
                   f"f={f_lo:9.3f}..{f_hi:9.3f} MHz")
    return out


# --------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------
def _save(fig, path: Path, dpi: int = 130) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}")


def plot_flag_matrix(beams: List[BeamStats], outpath: Path) -> None:
    nbeams = len(beams)
    nchans = beams[0].diag.nchans
    flag_matrix = np.stack([b.flagged for b in beams])  # (nbeams, nchans)
    freqs = beams[0].diag.freqs_mhz

    fig, ax = plt.subplots(figsize=(min(2 + 0.005 * nchans, 18),
                                    1.2 + 0.4 * nbeams))
    extent = [0, nchans, nbeams, 0]
    ax.imshow(flag_matrix, aspect="auto", interpolation="nearest",
              cmap="Greys", vmin=0, vmax=1, extent=extent)
    ax.set_yticks(np.arange(nbeams) + 0.5)
    ax.set_yticklabels([b.name for b in beams])
    ax.set_xlabel("channel index")
    ax.set_title(f"per-channel flags ({int(flag_matrix.sum())} flag-pixels "
                 f"across {nbeams} beams x {nchans} channels)")

    # secondary axis: frequency
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    n_ticks = 6
    tick_chans = np.linspace(0, nchans - 1, n_ticks).astype(int)
    ax2.set_xticks(tick_chans)
    ax2.set_xticklabels([f"{freqs[i]:.0f}" for i in tick_chans])
    ax2.set_xlabel("frequency (MHz)")

    _save(fig, outpath)


def plot_n_beams_per_chan(beams: List[BeamStats], outpath: Path) -> None:
    nbeams = len(beams)
    flag_matrix = np.stack([b.flagged for b in beams])
    n_flagged = flag_matrix.sum(axis=0)
    freqs = beams[0].diag.freqs_mhz

    fig, ax = plt.subplots(figsize=(11, 3.2))
    ax.fill_between(freqs, 0, n_flagged, step="mid", color="0.4", alpha=0.7)

    # overlay common (=nbeams) and unique (=1) channels
    common = (n_flagged == nbeams)
    unique = (n_flagged == 1)
    if common.any():
        ax.scatter(freqs[common], n_flagged[common], s=12,
                   color="C3", label=f"common ({int(common.sum())} ch)")
    if unique.any():
        ax.scatter(freqs[unique], n_flagged[unique], s=12,
                   color="C1", label=f"unique-to-1-beam ({int(unique.sum())} ch)")
    ax.set_ylim(-0.4, nbeams + 0.4)
    ax.set_yticks(np.arange(nbeams + 1))
    ax.set_xlabel("frequency (MHz)")
    ax.set_ylabel("# beams flagging")
    ax.set_title("# beams flagging each channel")
    if common.any() or unique.any():
        ax.legend(loc="upper right", fontsize=9, framealpha=0.85)
    ax.grid(alpha=0.3)
    _save(fig, outpath)


def plot_mean_rms_overlay(beams: List[BeamStats], outpath: Path) -> None:
    nbeams = len(beams)
    freqs = beams[0].diag.freqs_mhz
    flag_matrix = np.stack([b.flagged for b in beams])
    n_flagged = flag_matrix.sum(axis=0)
    common = (n_flagged == nbeams)
    unique = (n_flagged == 1)

    fig, axs = plt.subplots(2, 1, figsize=(11, 5.2), sharex=True)

    cmap = plt.get_cmap("tab10" if nbeams <= 10 else "tab20")
    for i, b in enumerate(beams):
        c = cmap(i % cmap.N)
        axs[0].plot(freqs, b.mean, lw=0.8, alpha=0.85, color=c, label=b.name)
        axs[1].plot(freqs, b.rms,  lw=0.8, alpha=0.85, color=c)

    # vertical bands for common (red) and unique (orange)
    for ax in axs:
        for (a, e) in _runs(common):
            ax.axvspan(freqs[a], freqs[e - 1], color="C3", alpha=0.15, lw=0)
        for (a, e) in _runs(unique):
            ax.axvspan(freqs[a], freqs[e - 1], color="C1", alpha=0.15, lw=0)

    axs[0].set_ylabel("per-channel mean")
    axs[1].set_ylabel("per-channel rms")
    axs[1].set_xlabel("frequency (MHz)")
    axs[0].set_title("per-beam spectra; red=flagged in all beams, "
                     "orange=flagged in exactly one beam")
    axs[0].legend(loc="upper right", fontsize=8, ncol=min(nbeams, 4),
                  framealpha=0.85)
    for ax in axs:
        ax.grid(alpha=0.3)
    _save(fig, outpath)


def plot_summary_counts(beams: List[BeamStats], outpath: Path) -> None:
    nbeams = len(beams)
    flag_matrix = np.stack([b.flagged for b in beams])
    n_flagged = flag_matrix.sum(axis=0)
    counts = np.bincount(n_flagged, minlength=nbeams + 1)
    nchans = beams[0].diag.nchans

    fig, ax = plt.subplots(figsize=(7, 3.5))
    bars = ax.bar(np.arange(nbeams + 1), counts, color="0.55", edgecolor="k")
    if nbeams >= 1:
        bars[0].set_color("C2")
        bars[-1].set_color("C3")
        if nbeams >= 2:
            bars[1].set_color("C1")
    ax.set_xticks(np.arange(nbeams + 1))
    ax.set_xlabel("# beams flagging")
    ax.set_ylabel("# channels")
    ax.set_title(f"channel commonality ({nchans} channels, {nbeams} beams)\n"
                 f"green=clean, orange=unique, red=common")
    for x, h in enumerate(counts):
        if h > 0:
            ax.text(x, h, f"{h}", ha="center", va="bottom", fontsize=9)
    ax.grid(alpha=0.3, axis="y")
    _save(fig, outpath)


# --------------------------------------------------------------------------
# CSV
# --------------------------------------------------------------------------
def write_csv(beams: List[BeamStats], outpath: Path) -> None:
    nchans = beams[0].diag.nchans
    freqs = beams[0].diag.freqs_mhz
    flag_matrix = np.stack([b.flagged for b in beams])
    n_flagged = flag_matrix.sum(axis=0)

    outpath.parent.mkdir(parents=True, exist_ok=True)
    cols = ["chan", "freq_MHz", "n_beams_flagged", "beams_flagged"]
    for b in beams:
        cols.append(f"mean_{b.name}")
        cols.append(f"rms_{b.name}")
        cols.append(f"flag_{b.name}")
    with open(outpath, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for c in range(nchans):
            flagged_names = ";".join(b.name for i, b in enumerate(beams)
                                     if flag_matrix[i, c])
            row = [c, f"{freqs[c]:.6f}", int(n_flagged[c]), flagged_names]
            for b in beams:
                row.append(f"{b.mean[c]:.6g}")
                row.append(f"{b.rms[c]:.6g}")
                row.append(int(bool(b.flagged[c])))
            w.writerow(row)
    print(f"  wrote {outpath}")


# --------------------------------------------------------------------------
# Zero-DM cross-beam analysis
# --------------------------------------------------------------------------
@dataclass
class ZdmStats:
    name: str
    raw: np.ndarray       # decimated zero-DM time series
    znorm: np.ndarray     # (raw - median) / (1.4826 * MAD)
    t: np.ndarray         # bin centre times (s, relative to file start)
    median: float
    mad: float            # median absolute deviation (raw units, NOT *1.4826)
    bin_sec: float        # actual bin size used (factor * tsamp)


def _decimate_mean(x: np.ndarray, factor: int) -> np.ndarray:
    """Block-mean decimate `x` by `factor`. Drops the trailing partial block."""
    if factor <= 1:
        return np.asarray(x, dtype=np.float64)
    n = (len(x) // factor) * factor
    if n == 0:
        return np.asarray(x[:0], dtype=np.float64)
    return x[:n].astype(np.float64).reshape(-1, factor).mean(axis=1)


def analyse_beam_zdm(beam: BeamStats, bin_sec: float):
    """Bin-decimate zero-DM, compute MAD-normalised z. Returns None if no Z0DM."""
    z = beam.diag.zerodm
    if z is None or len(z) == 0:
        return None
    factor = max(1, int(round(bin_sec / beam.diag.tsamp)))
    raw = _decimate_mean(z, factor)
    if raw.size == 0:
        return None
    med = float(np.median(raw))
    mad = float(np.median(np.abs(raw - med)))
    sigma = mad * 1.4826
    if sigma == 0.0:
        # constant series; avoid div-by-zero
        znorm = np.zeros_like(raw)
    else:
        znorm = (raw - med) / sigma
    actual_bin = factor * beam.diag.tsamp
    t = beam.diag.t_offset + (np.arange(len(raw)) + 0.5) * actual_bin
    return ZdmStats(name=beam.name, raw=raw, znorm=znorm, t=t,
                    median=med, mad=mad, bin_sec=actual_bin)


def _truncate_to_min(zdms: List["ZdmStats"]) -> int:
    """Return the common length to use across beams (min of all)."""
    return min(z.raw.size for z in zdms)


def correlation_matrix(zdms: List["ZdmStats"]) -> np.ndarray:
    n = _truncate_to_min(zdms)
    stack = np.stack([z.znorm[:n] for z in zdms])  # (nbeams, nT)
    if n < 2:
        return np.eye(stack.shape[0])
    return np.corrcoef(stack)


def find_coincident_bursts(zdms: List["ZdmStats"], sigma: float):
    """Group time bins where any beam exceeds +sigma into contiguous events.

    Returns a list of dicts with:
        t0, t1     -- start/stop time (s) of the run
        nbins      -- run length in bins
        n_beams    -- # beams crossing the threshold somewhere in the run
        beam_max   -- (nbeams,) peak z within the run, per beam
        beams_hit  -- list of beam names that crossed threshold
    """
    if not zdms:
        return []
    n = _truncate_to_min(zdms)
    stack = np.stack([z.znorm[:n] for z in zdms])  # (nbeams, nT)
    above = stack > sigma                          # (nbeams, nT)
    any_above = above.any(axis=0)
    runs = _runs(any_above)
    t = zdms[0].t[:n]
    bin_dt = zdms[0].bin_sec
    events = []
    for (a, b) in runs:
        seg = stack[:, a:b]
        beam_max = seg.max(axis=1)
        beams_hit = [z.name for k, z in enumerate(zdms) if (above[k, a:b]).any()]
        events.append({
            "t0": float(t[a] - 0.5 * bin_dt),
            "t1": float(t[b - 1] + 0.5 * bin_dt),
            "nbins": int(b - a),
            "n_beams": int(len(beams_hit)),
            "beam_max": beam_max,
            "beams_hit": beams_hit,
        })
    return events


# --------------------------------------------------------------------------
# Zero-DM plots
# --------------------------------------------------------------------------
def plot_zdm_overlay(zdms: List["ZdmStats"], outpath: Path,
                     stack_offset_sigma: float = 8.0) -> None:
    fig, ax = plt.subplots(figsize=(11, max(2.5, 0.7 * len(zdms) + 1.5)))
    cmap = plt.get_cmap("tab10" if len(zdms) <= 10 else "tab20")
    for i, z in enumerate(zdms):
        c = cmap(i % cmap.N)
        ax.plot(z.t, z.znorm + i * stack_offset_sigma,
                lw=0.6, alpha=0.9, color=c)
        ax.text(z.t[0], i * stack_offset_sigma + 0.5, f" {z.name}",
                color=c, fontsize=9, va="bottom")
    ax.set_xlabel("time since file start (s)")
    ax.set_ylabel("MAD-z (each beam offset by "
                  f"{stack_offset_sigma:g} sigma)")
    ax.set_title(f"zero-DM time series across beams "
                 f"(bin = {zdms[0].bin_sec*1e3:.3g} ms)")
    ax.grid(alpha=0.3)
    _save(fig, outpath)


def plot_zdm_correlation(zdms: List["ZdmStats"], outpath: Path) -> None:
    corr = correlation_matrix(zdms)
    nbeams = corr.shape[0]
    names = [z.name for z in zdms]

    fig, ax = plt.subplots(figsize=(0.7 + 0.5 * nbeams, 0.7 + 0.5 * nbeams))
    im = ax.imshow(corr, vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(nbeams)); ax.set_xticklabels(names, rotation=45, ha="right")
    ax.set_yticks(range(nbeams)); ax.set_yticklabels(names)
    for i in range(nbeams):
        for j in range(nbeams):
            ax.text(j, i, f"{corr[i, j]:+.2f}", ha="center", va="center",
                    fontsize=8,
                    color="white" if abs(corr[i, j]) > 0.6 else "black")
    ax.set_title(f"zero-DM Pearson correlation\n(bin = {zdms[0].bin_sec*1e3:.3g} ms)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    _save(fig, outpath)


def plot_zdm_sigma_exceedance(zdms: List["ZdmStats"], outpath: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ks = np.linspace(0.0, 8.0, 81)
    cmap = plt.get_cmap("tab10" if len(zdms) <= 10 else "tab20")
    for i, z in enumerate(zdms):
        n = z.znorm.size
        if n == 0:
            continue
        # one-sided exceedance fraction (positive z); zero-DM is positive-definite
        # signal of interest (CLT-summed channels), so positive tail is what matters.
        fracs = np.array([(z.znorm > k).sum() / n for k in ks])
        # clamp zeros so log scale doesn't crash
        floor = 0.5 / max(n, 1)
        fracs = np.where(fracs > 0, fracs, floor)
        ax.semilogy(ks, fracs, lw=1.0, alpha=0.9,
                    color=cmap(i % cmap.N), label=z.name)
    # Gaussian reference
    from math import erfc, sqrt
    gauss = np.array([0.5 * erfc(k / sqrt(2.0)) for k in ks])
    ax.semilogy(ks, np.maximum(gauss, 1e-12), "k--", lw=1.0,
                label="Gaussian (one-sided)")
    ax.set_xlabel("k  (MAD-sigma)")
    ax.set_ylabel("fraction of samples with z > k")
    ax.set_title("zero-DM positive-tail exceedance, per beam")
    ax.set_ylim(bottom=floor / 2 if zdms else 1e-9)
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="upper right", fontsize=8, ncol=2)
    _save(fig, outpath)


def plot_zdm_burst_coincidence(zdms: List["ZdmStats"], outpath: Path,
                               sigma: float) -> None:
    n = _truncate_to_min(zdms)
    stack = np.stack([z.znorm[:n] for z in zdms])  # (nbeams, nT)
    above = stack > sigma
    n_above = above.sum(axis=0)
    t = zdms[0].t[:n]

    fig, axs = plt.subplots(2, 1, figsize=(11, 5.5),
                            gridspec_kw=dict(height_ratios=[3, 1]),
                            sharex=True)
    cmap = plt.get_cmap("tab10" if len(zdms) <= 10 else "tab20")
    for i, z in enumerate(zdms):
        axs[0].plot(t, z.znorm[:n], lw=0.5, alpha=0.85,
                    color=cmap(i % cmap.N), label=z.name)
    axs[0].axhline(sigma, color="k", ls="--", lw=0.7, alpha=0.7,
                   label=f"{sigma:g}-sigma")
    axs[0].set_ylabel("MAD-z (zero-DM)")
    axs[0].legend(loc="upper right", fontsize=8, ncol=min(len(zdms), 4))
    axs[0].grid(alpha=0.3)
    axs[0].set_title(f"zero-DM bursts, threshold = {sigma:g} sigma  "
                     f"(bin = {zdms[0].bin_sec*1e3:.3g} ms)")

    nbeams = len(zdms)
    axs[1].fill_between(t, 0, n_above, step="mid", color="0.4", alpha=0.7)
    axs[1].set_ylim(-0.4, nbeams + 0.4)
    axs[1].set_yticks(np.arange(nbeams + 1))
    axs[1].set_xlabel("time since file start (s)")
    axs[1].set_ylabel("# beams\nabove thresh")
    axs[1].grid(alpha=0.3)
    _save(fig, outpath)


# --------------------------------------------------------------------------
# Spatial-filter feasibility (Kocz et al. 2010 in total-power)
# --------------------------------------------------------------------------
#
# The Kocz method does SVD of the complex voltage covariance C(f) =
# <S_i(f) S_j(f)*> per spectral channel, then nulls / projects out the
# top q eigenvalues (the RFI subspace). We don't have voltages, only
# detected total-power intensities. The total-power analog is to form,
# for each channel c, the (nbeams, nchunks) matrix M_c of mean-subtracted
# beam intensities, and take the SVD of (M_c M_c^T) / nchunks. This is
# the cross-beam intensity covariance: if the RFI is rank-1 across beams
# (one common interferer hitting every beam with an arbitrary per-beam
# gain) the eigenvalue spectrum will have lambda_1 >> lambda_2..p, just
# like in the voltage case. If the eigenvalues are flat, the RFI either
# isn't actually correlated across beams or it's lost in noise -- and
# voltage-domain Kocz won't help in either case.
#
# Whitening (paper Eq. 11): the SVD partitioning into "interference" and
# "noise" subspaces only holds if the system noise is spatially white,
# i.e. <|N_i|^2> = sigma_N^2 for all beams. We check this by looking at
# the spread of per-channel medians and per-beam noise levels.
# --------------------------------------------------------------------------
@dataclass
class SpatialStats:
    eigvals_per_chan: np.ndarray   # (nchans, nbeams), descending sorted
    eig_ratio: np.ndarray          # (nchans,) lambda_1 / median(lambda_2..p)
    mean_levels: np.ndarray        # (nbeams,) median of mean spectrum per beam
    noise_levels: np.ndarray       # (nbeams,) median CRMS over chunks+chans
    rank_resid_per_chan: np.ndarray  # (nchans, nbeams+1), frac of variance left
    zdm_eigvals: np.ndarray        # (nbeams,) eigenvalues of zero-DM cov matrix
    zdm_pcs: np.ndarray            # (nbeams, nT) principal-component time series
    zdm_t: np.ndarray              # (nT,)
    zdm_explained: np.ndarray      # (nbeams,) cumulative variance fraction


def _stack_chunk_means(beams: List[BeamStats]) -> np.ndarray:
    """Return (nbeams, nchunks, nchans) cube of CMEA, truncated to common
    nchunks across beams."""
    n = min(b.diag.nchunks for b in beams)
    arr = np.stack([b.diag.cmean[:n].astype(np.float64) for b in beams])
    return arr  # (nbeams, nchunks, nchans)


def _per_beam_noise_level(beams: List[BeamStats]) -> np.ndarray:
    """Robust per-beam noise estimate: median of CRMS across chunks and
    channels. For 2-bit uniform data this is ~1.118; deviations indicate
    gain/normalisation differences across beams."""
    return np.array([float(np.median(b.diag.crms)) for b in beams])


def analyse_spatial_feasibility(beams: List[BeamStats],
                                whiten: bool = True) -> SpatialStats:
    """Per-channel SVD of cross-beam chunk-mean intensities, plus PCA
    of the cross-beam zero-DM time series.

    For each channel c:
        M = (CMEA[:, :, c] - mean_over_chunks_per_beam)            (p x N)
        if whiten: M /= rms_over_chunks_per_beam                   (per-beam)
        eigenvalues of (M @ M.T) / nchunks                         (p,)
        rank_k_residual[k] = (sum of eigvals[k:]) / sum(eigvals)
    """
    nbeams = len(beams)
    cube = _stack_chunk_means(beams)            # (nbeams, nchunks, nchans)
    nch_eff = cube.shape[1]
    nchans = cube.shape[2]

    # center each (beam, channel) time series
    M = cube - cube.mean(axis=1, keepdims=True)
    if whiten:
        s = M.std(axis=1, keepdims=True)
        s[s == 0.0] = 1.0
        M = M / s

    # cross-beam covariance matrix per channel: shape (nchans, nbeams, nbeams)
    # einsum is fast and avoids a python loop over channels.
    C = np.einsum("itc,jtc->cij", M, M) / max(nch_eff, 1)

    # eigvalsh is for Hermitian/symmetric matrices; returns ascending order
    eigs = np.linalg.eigvalsh(C)                # (nchans, nbeams)
    eigs = eigs[:, ::-1]                        # descending
    eigs = np.maximum(eigs, 0.0)

    # rank-k residual variance fraction
    cum = np.cumsum(eigs, axis=1)               # (nchans, nbeams)
    total = cum[:, -1:].copy()
    total[total == 0] = 1.0
    rank_resid = np.concatenate([np.ones((nchans, 1)),
                                 1.0 - cum / total], axis=1)  # (nchans, nbeams+1)
    rank_resid = np.clip(rank_resid, 0.0, 1.0)

    # diagnostic ratio: top eigenvalue over median of the rest. A clean
    # channel has all eigenvalues ~1 (after whitening) so ratio ~ 1; a
    # rank-1 RFI channel has lambda_1 >> lambda_2..p so ratio >> 1.
    if nbeams >= 2:
        rest_med = np.median(eigs[:, 1:], axis=1)
        rest_med[rest_med == 0] = np.nan
        eig_ratio = eigs[:, 0] / rest_med
    else:
        eig_ratio = np.full(nchans, np.nan)

    mean_levels = np.array([float(np.median(b.mean)) for b in beams])
    noise_levels = _per_beam_noise_level(beams)

    # ---- zero-DM PCA (single covariance over whole observation) -----------
    # We re-derive it from each beam's full-rate Z0DM if available, decimated
    # to a coarser bin to keep memory sane. If Z0DM is missing for any beam
    # we fall back to the chunk-mean sum across channels.
    zdms_raw = []
    for b in beams:
        if b.diag.zerodm is not None:
            zdms_raw.append(b.diag.zerodm.astype(np.float64))
        else:
            zdms_raw.append(b.diag.cmean.sum(axis=1).astype(np.float64))
    n_t = min(z.size for z in zdms_raw)
    Z = np.stack([z[:n_t] for z in zdms_raw])    # (nbeams, n_t)
    # decimate to keep computation tractable
    target_n = 5000
    factor = max(1, n_t // target_n)
    if factor > 1:
        m = (n_t // factor) * factor
        Z = Z[:, :m].reshape(nbeams, -1, factor).mean(axis=2)
    Z = Z - Z.mean(axis=1, keepdims=True)
    s = Z.std(axis=1, keepdims=True)
    s[s == 0] = 1.0
    Zw = Z / s
    Czdm = (Zw @ Zw.T) / max(Zw.shape[1], 1)
    w, V = np.linalg.eigh(Czdm)                  # ascending
    w = w[::-1]
    V = V[:, ::-1]
    w = np.maximum(w, 0.0)
    pcs = V.T @ Zw                                # (nbeams, nT_dec)
    explained = np.cumsum(w) / max(w.sum(), 1.0)

    # build a time axis for the decimated zero-DM (best-effort)
    tsamp = beams[0].diag.tsamp
    if any(b.diag.zerodm is not None for b in beams):
        bin_sec = factor * tsamp
        t_axis = beams[0].diag.t_offset + (np.arange(Zw.shape[1]) + 0.5) * bin_sec
    else:
        chunk_sec = beams[0].diag.chunk_sec
        bin_sec = factor * chunk_sec
        t_axis = beams[0].diag.t_offset + (np.arange(Zw.shape[1]) + 0.5) * bin_sec

    return SpatialStats(
        eigvals_per_chan=eigs, eig_ratio=eig_ratio,
        mean_levels=mean_levels, noise_levels=noise_levels,
        rank_resid_per_chan=rank_resid,
        zdm_eigvals=w, zdm_pcs=pcs, zdm_t=t_axis, zdm_explained=explained,
    )


# --- spatial-filter plots --------------------------------------------------
def plot_spatial_eigvals(beams, st: SpatialStats, outpath: Path) -> None:
    nchans, nbeams = st.eigvals_per_chan.shape
    freqs = beams[0].diag.freqs_mhz
    eigs = st.eigvals_per_chan.copy()
    eigs[eigs <= 0] = np.nan

    fig, ax = plt.subplots(figsize=(11, 3.0 + 0.18 * nbeams))
    im = ax.imshow(np.log10(eigs.T), aspect="auto", origin="lower",
                   extent=[0, nchans, 0.5, nbeams + 0.5],
                   interpolation="nearest", cmap="magma")
    ax.set_yticks(np.arange(1, nbeams + 1))
    ax.set_ylabel("eigenvalue index (1=largest)")
    ax.set_xlabel("channel index")
    ax.set_title("per-channel cross-beam SVD spectrum (Kocz Fig. 2 analog)\n"
                 "log10(eigvalue);  rank-1 RFI -> bright top row, dim below")
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("log10 eigenvalue")
    # secondary frequency axis
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    n_ticks = 6
    tick_chans = np.linspace(0, nchans - 1, n_ticks).astype(int)
    ax2.set_xticks(tick_chans)
    ax2.set_xticklabels([f"{freqs[i]:.0f}" for i in tick_chans])
    ax2.set_xlabel("frequency (MHz)")
    _save(fig, outpath)


def plot_spatial_eigratio(beams, st: SpatialStats, outpath: Path) -> None:
    freqs = beams[0].diag.freqs_mhz
    fig, ax = plt.subplots(figsize=(11, 3.2))
    ratio = st.eig_ratio.copy()
    ax.semilogy(freqs, np.where(np.isfinite(ratio), ratio, np.nan),
                lw=0.8, color="C3")
    ax.axhline(1.0, color="k", ls=":", lw=0.7, label="flat (no spatial structure)")
    ax.axhline(5.0, color="0.4", ls="--", lw=0.7,
               label="spatial filtering useful (>5x)")
    ax.set_xlabel("frequency (MHz)")
    ax.set_ylabel(r"$\lambda_1 / \mathrm{median}(\lambda_{2..p})$")
    ax.set_title("per-channel cross-beam dominance ratio\n"
                 "(predicts where Kocz-style voltage SVD would help)")
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.85)
    _save(fig, outpath)


def plot_spatial_suppression(beams, st: SpatialStats, outpath: Path,
                             rank: int) -> None:
    """Pick the worst-RFI, median, and cleanest channel; show what fraction
    of variance survives after rank-k subtraction, k = 0..nbeams. This is
    the total-power analog of Kocz Tables 2-3 (RFI-on vs RFI-off S/N)."""
    freqs = beams[0].diag.freqs_mhz
    ratio = st.eig_ratio
    finite = np.where(np.isfinite(ratio))[0]
    if finite.size == 0:
        return
    order = finite[np.argsort(ratio[finite])]
    cleanest = order[0]
    worst = order[-1]
    mid = order[len(order) // 2]
    picks = [("worst RFI", worst), ("median", mid), ("cleanest", cleanest)]

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    nbeams = st.eigvals_per_chan.shape[1]
    ks = np.arange(nbeams + 1)
    for label, c in picks:
        f = freqs[c]
        r = st.rank_resid_per_chan[c]
        ax.plot(ks, r, marker="o", lw=1.2,
                label=f"ch {c} ({f:.1f} MHz, ratio={ratio[c]:.2g})  [{label}]")
    ax.axvline(rank, color="0.4", ls=":", lw=1.0,
               label=f"chosen rank-k = {rank}")
    ax.set_xlabel("rank-k subtracted")
    ax.set_ylabel("residual variance fraction")
    ax.set_title("rank-k spatial-filter residual\n"
                 "(low at small k => RFI is low-rank => filter would work)")
    ax.set_yscale("log")
    ax.grid(alpha=0.3, which="both")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.85)
    _save(fig, outpath)


def plot_spatial_zdm_pca(beams, st: SpatialStats, outpath: Path) -> None:
    nbeams = st.zdm_eigvals.size
    fig, axs = plt.subplots(2, 1, figsize=(11, 5.4),
                            gridspec_kw=dict(height_ratios=[1, 2]))

    # top: bar chart of normalised eigenvalues
    eigs = st.zdm_eigvals / max(st.zdm_eigvals.sum(), 1.0)
    bars = axs[0].bar(np.arange(1, nbeams + 1), eigs, color="0.55", edgecolor="k")
    bars[0].set_color("C3")
    axs[0].set_xlabel("eigenvalue index"); axs[0].set_ylabel("variance fraction")
    axs[0].set_xticks(np.arange(1, nbeams + 1))
    axs[0].set_title("zero-DM cross-beam PCA  "
                     f"(top PC explains {100*eigs[0]:.1f}% of variance)")
    axs[0].grid(alpha=0.3, axis="y")

    # bottom: top-3 PC time series (or fewer if nbeams < 3)
    npc = min(3, nbeams)
    cmap = plt.get_cmap("tab10")
    for k in range(npc):
        pc = st.zdm_pcs[k]
        # offset for clarity
        axs[1].plot(st.zdm_t,
                    pc / (np.std(pc) if np.std(pc) > 0 else 1.0)
                    + 6 * (npc - 1 - k),
                    lw=0.7, color=cmap(k),
                    label=f"PC{k+1} ({100*eigs[k]:.1f}% var)")
    axs[1].set_xlabel("time since file start (s)")
    axs[1].set_ylabel("normalised PC amplitude (offset)")
    axs[1].legend(loc="upper right", fontsize=9, framealpha=0.85)
    axs[1].grid(alpha=0.3)
    _save(fig, outpath)


def plot_spatial_whitening(beams, st: SpatialStats, outpath: Path) -> None:
    names = [b.name for b in beams]
    nbeams = len(beams)
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.4))
    axs[0].bar(np.arange(nbeams), st.mean_levels, color="C0")
    axs[0].set_xticks(np.arange(nbeams)); axs[0].set_xticklabels(names, rotation=45, ha="right")
    axs[0].set_ylabel("median(mean spectrum)")
    axs[0].set_title("per-beam mean intensity")
    axs[0].grid(alpha=0.3, axis="y")

    axs[1].bar(np.arange(nbeams), st.noise_levels, color="C2")
    axs[1].set_xticks(np.arange(nbeams)); axs[1].set_xticklabels(names, rotation=45, ha="right")
    axs[1].set_ylabel("median(CRMS)")
    axs[1].set_title("per-beam noise level "
                     f"(spread = {100*st.noise_levels.std()/max(st.noise_levels.mean(),1e-9):.1f}%)")
    axs[1].grid(alpha=0.3, axis="y")
    fig.suptitle("whitening diagnostic (Kocz Eq. 11 requires <|N_i|^2> ~ const)",
                 fontsize=10)
    _save(fig, outpath)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("diags", nargs="+", help="one or more .diag files")
    p.add_argument("--outdir", default="compare_rfi", type=Path,
                   help="output directory (default: ./compare_rfi)")
    p.add_argument("--names", default=None,
                   help="comma-separated beam labels; default = filename stems")
    p.add_argument("--nsigma", type=float, default=5.0,
                   help="MAD-z threshold for flagging (default: 5.0)")
    p.add_argument("--no-mean-flag", action="store_true",
                   help="do not flag on per-channel mean z-score")
    p.add_argument("--no-rms-flag", action="store_true",
                   help="do not flag on per-channel rms z-score")
    p.add_argument("--top-k", type=int, default=10,
                   help="print top-K longest bad bands per beam (default: 10)")
    p.add_argument("--csv", default=None, type=Path,
                   help="path for per-channel CSV (default: <outdir>/per_chan.csv)")
    p.add_argument("--zdm-bin-sec", type=float, default=1e-3,
                   help="time bin (sec) for zero-DM cross-beam stats (default: 1e-3)")
    p.add_argument("--zdm-burst-sigma", type=float, default=6.0,
                   help="MAD-z threshold for zero-DM burst coincidence (default: 6.0)")
    p.add_argument("--no-zdm", action="store_true",
                   help="skip zero-DM cross-beam analysis even if Z0DM is present")
    p.add_argument("--no-spatial", action="store_true",
                   help="skip Kocz-style spatial-filter feasibility analysis")
    p.add_argument("--spatial-rank", type=int, default=0,
                   help="rank-k subtraction depth for spatial suppression curve "
                        "(default: max(1, nbeams // 4))")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    paths = [Path(p) for p in args.diags]
    if len(paths) < 2:
        print("warning: comparing across <2 beams is not very interesting; "
              "the script will still produce per-beam stats.", file=sys.stderr)

    if args.names:
        names = [s.strip() for s in args.names.split(",")]
        if len(names) != len(paths):
            sys.exit(f"--names has {len(names)} entries but {len(paths)} files were given")
    else:
        names = [p.stem for p in paths]

    print(f"loading {len(paths)} .diag files ...")
    beams: List[BeamStats] = []
    for path, name in zip(paths, names):
        b = analyse_beam(path, name, nsigma=args.nsigma,
                         use_mean=not args.no_mean_flag,
                         use_rms=not args.no_rms_flag)
        d = b.diag
        print(f"  {name:>12s}  nchans={d.nchans:5d}  nsamples={d.nsamples:>10d}  "
              f"tobs={d.nsamples * d.tsamp:6.2f}s  "
              f"flagged={int(b.flagged.sum())}/{d.nchans} ch")
        beams.append(b)
    _check_compatible(beams)

    nbeams = len(beams)
    nchans = beams[0].diag.nchans
    freqs = beams[0].diag.freqs_mhz
    flag_matrix = np.stack([b.flagged for b in beams])
    n_flagged = flag_matrix.sum(axis=0)

    # ---- print summary -----------------------------------------------------
    print(f"\nchannel-commonality summary  (nbeams={nbeams}, nchans={nchans}):")
    counts = np.bincount(n_flagged, minlength=nbeams + 1)
    print(f"  clean (0 beams)            : {counts[0]:5d} ch  "
          f"({100.0 * counts[0] / nchans:5.1f}%)")
    if nbeams >= 2:
        print(f"  unique to one beam         : {counts[1]:5d} ch  "
              f"({100.0 * counts[1] / nchans:5.1f}%)")
        shared = int(counts[2:nbeams].sum()) if nbeams >= 3 else 0
        if nbeams >= 3:
            print(f"  shared by 2..{nbeams - 1} beams       : {shared:5d} ch  "
                  f"({100.0 * shared / nchans:5.1f}%)")
        print(f"  common (all {nbeams} beams)        : {counts[nbeams]:5d} ch  "
              f"({100.0 * counts[nbeams] / nchans:5.1f}%)")

    # ---- per-beam: top bad bands ------------------------------------------
    print(f"\nper-beam bad bands (top-{args.top_k} longest contiguous):")
    for b in beams:
        n = int(b.flagged.sum())
        print(f"\n  [{b.name}]  flagged {n} / {b.diag.nchans} channels")
        if n == 0:
            continue
        lines = _format_runs(b.flagged, freqs, args.top_k)
        for s in lines:
            print(s)

    # ---- per-beam: uniquely-flagged channels ------------------------------
    if nbeams >= 2:
        print(f"\nchannels flagged ONLY by one beam:")
        unique_mask = (n_flagged == 1)
        for i, b in enumerate(beams):
            mine = b.flagged & unique_mask
            n = int(mine.sum())
            print(f"\n  [{b.name}]  unique to this beam: {n} ch")
            if n == 0:
                continue
            for s in _format_runs(mine, freqs, args.top_k):
                print(s)

    # ---- common channels ---------------------------------------------------
    if nbeams >= 2:
        common = (n_flagged == nbeams)
        n = int(common.sum())
        print(f"\nchannels flagged by ALL {nbeams} beams: {n} ch")
        for s in _format_runs(common, freqs, args.top_k):
            print(s)

    # ---- plots + csv -------------------------------------------------------
    print(f"\nwriting plots to {args.outdir}/ ...")
    plot_flag_matrix       (beams, args.outdir / "flag_matrix.png")
    plot_n_beams_per_chan  (beams, args.outdir / "n_beams_per_chan.png")
    plot_mean_rms_overlay  (beams, args.outdir / "mean_rms_overlay.png")
    plot_summary_counts    (beams, args.outdir / "summary_counts.png")
    csv_path = args.csv if args.csv else (args.outdir / "per_chan.csv")
    write_csv(beams, csv_path)

    # ---- zero-DM cross-beam analysis --------------------------------------
    if args.no_zdm:
        print("\n--no-zdm: skipping zero-DM cross-beam analysis.")
        print("done.")
        return
    zdms: List[ZdmStats] = []
    missing = []
    for b in beams:
        zs = analyse_beam_zdm(b, args.zdm_bin_sec)
        if zs is None:
            missing.append(b.name)
        else:
            zdms.append(zs)
    if missing:
        print(f"\nzero-DM missing for {len(missing)} beam(s): {','.join(missing)} "
              f"(was rfidiag run with -Z?)")
    if len(zdms) < 1:
        print("no beams have Z0DM; skipping zero-DM analysis.")
        print("done.")
        return

    # alignment sanity check
    n_common = _truncate_to_min(zdms)
    n_per_beam = [z.raw.size for z in zdms]
    if len(set(n_per_beam)) > 1:
        print(f"warning: zero-DM bin counts differ across beams "
              f"({n_per_beam}); truncating to {n_common}.")
    tstarts = [b.diag.tstart for b in beams if b.name in {z.name for z in zdms}]
    if max(tstarts) - min(tstarts) > 1e-9:
        print(f"warning: tstart differs across beams "
              f"(range {(max(tstarts)-min(tstarts))*86400.0:.6f} s); "
              f"time alignment is approximate.")

    print(f"\nzero-DM cross-beam summary "
          f"(bin = {zdms[0].bin_sec*1e3:.3g} ms, "
          f"{n_common} bins, {n_common * zdms[0].bin_sec:.3f} s):")
    sigma = args.zdm_burst_sigma
    print(f"  {'beam':>14s}  {'median':>10s}  {'mad':>8s}  "
          f"{'max_z':>7s}  {'#> '+f'{sigma:g}σ':>9s}")
    for z in zdms:
        n_burst = int((z.znorm > sigma).sum())
        max_z = float(z.znorm.max()) if z.znorm.size else 0.0
        print(f"  {z.name:>14s}  {z.median:10.3f}  {z.mad:8.3f}  "
              f"{max_z:7.2f}  {n_burst:9d}")

    if len(zdms) >= 2:
        corr = correlation_matrix(zdms)
        print("\npairwise zero-DM Pearson correlation:")
        head = "         " + "".join(f"{z.name:>10s}" for z in zdms)
        print(head)
        for i, zi in enumerate(zdms):
            row = f"  {zi.name:>7s}" + "".join(
                f"  {corr[i, j]:+.4f}" for j in range(len(zdms)))
            print(row)

    # burst events
    events = find_coincident_bursts(zdms, sigma)
    if events:
        # sort by max z across any beam, descending
        events.sort(key=lambda e: float(e["beam_max"].max()), reverse=True)
        print(f"\nburst events with z > {sigma:g} (any beam), "
              f"top {min(args.top_k, len(events))} of {len(events)} "
              f"(sorted by peak z):")
        names = [z.name for z in zdms]
        head = (f"  {'t (s)':>10s}  {'dur(ms)':>7s}  {'#beams':>6s}  "
                + "  ".join(f"{n:>7s}" for n in names))
        print(head)
        for e in events[:args.top_k]:
            t_mid = 0.5 * (e["t0"] + e["t1"])
            dur_ms = (e["t1"] - e["t0"]) * 1e3
            cells = "  ".join(
                (f"{e['beam_max'][k]:7.2f}*"
                 if names[k] in e["beams_hit"]
                 else f"{e['beam_max'][k]:7.2f} ")
                for k in range(len(names)))
            print(f"  {t_mid:10.4f}  {dur_ms:7.2f}  {e['n_beams']:6d}  {cells}")
        n_in_k = np.zeros(len(zdms) + 1, dtype=int)
        for e in events:
            n_in_k[e["n_beams"]] += 1
        print("  events seen in N beams:  " +
              "  ".join(f"N={k}: {n_in_k[k]}" for k in range(1, len(zdms) + 1)))
    else:
        print(f"\nno zero-DM samples cross {sigma:g}-sigma in any beam.")

    plot_zdm_overlay          (zdms, args.outdir / "zdm_overlay.png")
    plot_zdm_correlation      (zdms, args.outdir / "zdm_correlation.png")
    plot_zdm_sigma_exceedance (zdms, args.outdir / "zdm_sigma_exceedance.png")
    plot_zdm_burst_coincidence(zdms, args.outdir / "zdm_burst_coincidence.png",
                               sigma=sigma)

    # ---- spatial-filter feasibility (Kocz et al. 2010) --------------------
    if args.no_spatial:
        print("\n--no-spatial: skipping spatial-filter feasibility analysis.")
        print("done.")
        return
    if nbeams < 2:
        print("\nspatial-filter analysis skipped: need at least 2 beams.")
        print("done.")
        return

    print(f"\nspatial-filter feasibility (Kocz et al. 2010 total-power analog):")
    st = analyse_spatial_feasibility(beams, whiten=True)

    # whitening diagnostic: gain spread across beams
    nm = st.noise_levels
    spread_pct = 100.0 * nm.std() / max(nm.mean(), 1e-9)
    mn = st.mean_levels
    mn_spread_pct = 100.0 * mn.std() / max(mn.mean(), 1e-9)
    print(f"  whitening (Kocz Eq. 11 requires <|N_i|^2> approx const across beams):")
    print(f"    per-beam noise level (median CRMS) : "
          f"min={nm.min():.4f}  max={nm.max():.4f}  spread={spread_pct:5.2f}%")
    print(f"    per-beam mean level (median CMEA)  : "
          f"min={mn.min():.4f}  max={mn.max():.4f}  spread={mn_spread_pct:5.2f}%")
    if spread_pct > 20.0:
        print(f"    >> WARNING: noise-level spread {spread_pct:.1f}% is large; "
              f"normalise gains before voltage-domain Kocz, else the noise "
              f"subspace will leak into the interferer eigenvalues.")

    # per-channel SVD verdict
    ratio = st.eig_ratio
    finite = np.isfinite(ratio)
    n_kocz_useful = int(((ratio > 5.0) & finite).sum())
    n_strong       = int(((ratio > 20.0) & finite).sum())
    n_clean        = int(((ratio < 1.5) & finite).sum())
    median_ratio   = float(np.nanmedian(ratio[finite])) if finite.any() else float("nan")
    max_ratio      = float(np.nanmax(ratio[finite])) if finite.any() else float("nan")
    nchans = beams[0].diag.nchans
    print(f"  per-channel cross-beam SVD (lambda_1 / median(lambda_2..p)):")
    print(f"    median across channels : {median_ratio:7.2f}")
    print(f"    max    across channels : {max_ratio:7.2f}")
    print(f"    #chans with ratio > 5  : {n_kocz_useful:5d} / {nchans} "
          f"({100.0*n_kocz_useful/nchans:5.1f}%)   "
          f"<- expect substantial Kocz gain in these")
    print(f"    #chans with ratio > 20 : {n_strong:5d} / {nchans} "
          f"({100.0*n_strong/nchans:5.1f}%)   "
          f"<- rank-1 RFI dominated, very strong gain expected")
    print(f"    #chans with ratio < 1.5: {n_clean:5d} / {nchans} "
          f"({100.0*n_clean/nchans:5.1f}%)   "
          f"<- clean; do NOT spatial-filter (Kocz Fig. 3 caveat)")

    # rank-k residual at the worst channel
    rank = args.spatial_rank if args.spatial_rank > 0 else max(1, nbeams // 4)
    if finite.any():
        worst_c = int(np.nanargmax(ratio))
        f_worst = beams[0].diag.freqs_mhz[worst_c]
        resid = st.rank_resid_per_chan[worst_c]
        # dB suppression at chosen rank (vs no subtraction)
        sup_db = 10.0 * np.log10(max(resid[rank], 1e-12) / max(resid[0], 1e-12))
        print(f"  rank-{rank} suppression at worst channel "
              f"(ch {worst_c}, {f_worst:.1f} MHz): "
              f"{sup_db:+6.1f} dB residual variance "
              f"({100.0*resid[rank]:5.2f}% left after rank-{rank})")

    # zero-DM PCA verdict
    zdm_top_var = 100.0 * st.zdm_eigvals[0] / max(st.zdm_eigvals.sum(), 1.0)
    print(f"  zero-DM PCA: top PC explains {zdm_top_var:5.1f}% of cross-beam "
          f"variance; rank-1 explains {100*st.zdm_explained[0]:5.1f}%, "
          f"rank-2 {100*st.zdm_explained[min(1, nbeams-1)]:5.1f}%, "
          f"rank-{rank} {100*st.zdm_explained[min(rank-1, nbeams-1)]:5.1f}%.")

    # short verdict
    print("\nverdict:")
    if n_kocz_useful >= max(1, nchans // 50) and spread_pct < 30.0:
        print(f"  Kocz-style voltage-domain spatial filtering is LIKELY TO HELP")
        print(f"  on {n_kocz_useful} channels ({100.0*n_kocz_useful/nchans:.1f}%) of the band.")
        print(f"  the worst-affected channels show rank-1 or rank-2 RFI in")
        print(f"  total-power proxy, which is the regime Kocz Fig. 4/Table 2 cover")
        print(f"  (their pre/post S/N ratios suggest 5-7x gain).")
    elif n_kocz_useful == 0:
        print(f"  Kocz-style spatial filtering is UNLIKELY TO HELP: no channels")
        print(f"  show rank-1 RFI in the cross-beam intensity matrix. RFI is")
        print(f"  either too weak, uncorrelated across beams, or already gone.")
    else:
        print(f"  marginal: {n_kocz_useful} channels show useful spatial structure")
        print(f"  but whitening spread ({spread_pct:.1f}%) is large; tighten")
        print(f"  gain calibration before relying on voltage-domain Kocz.")
    print("  caveat: this is the real-valued total-power proxy; actual voltage")
    print("    SVD may resolve interferers that overlap in intensity (because of")
    print("    phase) and may lose channels where AGC/saturation distorts gains")
    print("    (paper Section 5.1, Fig. 7).")

    # plots
    plot_spatial_eigvals     (beams, st, args.outdir / "spatial_eigvals.png")
    plot_spatial_eigratio    (beams, st, args.outdir / "spatial_eigratio.png")
    plot_spatial_suppression (beams, st, args.outdir / "spatial_suppression.png",
                              rank=rank)
    plot_spatial_zdm_pca     (beams, st, args.outdir / "spatial_zdm_pca.png")
    plot_spatial_whitening   (beams, st, args.outdir / "spatial_whitening.png")
    print("done.")


if __name__ == "__main__":
    main()
