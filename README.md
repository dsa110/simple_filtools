# simple_filtools

Standalone, dependency-free RFI diagnostics for SIGPROC-format filterbank data.
Designed around 2-bit, 2048-channel, ~50 us, ~40-min DSA-2000 search-mode
files (~25 GB packed), but works on any SIGPROC `nbits in {1,2,4,8}` file
with `nifs == 1`.

The C tool (`rfidiag`) streams the file in one pass, writes a single
self-describing binary diagnostics file (`*.diag`), and a companion Python
script (`python/plot_diagnostics.py`) turns that into a set of PNG plots.

The SIGPROC header reader is ported from
[`dsa110-sigproc`](../dsa110-sigproc) (`read_header.c`, `pack_unpack.c`,
`strings_equal.c`) so this repo has no external C dependencies beyond
`libc + libm`.

---

## Build

```bash
make            # builds ./rfidiag
make test       # builds tests, runs unpack unit test + numpy cross-check
make clean
```

Tests need Python with `numpy` (matplotlib is required only for plotting).

## Run

```bash
# 1) C streaming pass: produces obs.diag
./rfidiag obs.fil -o obs.diag

# 2) Plot to PNGs
python3 python/plot_diagnostics.py obs.diag --outdir obs_plots/
```

`./rfidiag --help` for full options. The most useful are:

| flag           | meaning                                         | default |
| -------------- | ----------------------------------------------- | ------- |
| `-s start_sec` | skip the first `start_sec` seconds              | 0       |
| `-d dur_sec`   | process only `dur_sec` seconds                  | all     |
| `-c chunk_sec` | chunk size for chunk-averaged stats / waterfall | 0.5     |
| `-Z`           | do **not** store full-rate zero-DM in `.diag`   | off     |
| `-q`           | quiet (no progress)                             | off     |

`plot_diagnostics.py --help` for plot-side options (`--fdec`, `--tdec`,
`--n-spectra`, `--spectra-sec`).

## Diagnostics produced

The `.diag` file holds (see `src/diag.h`):

| section | content                                    | dtype     | shape                 |
| ------- | ------------------------------------------ | --------- | --------------------- |
| `SUMC`  | per-channel sum                            | uint64    | `(nchans,)`           |
| `SUMQ`  | per-channel sum-of-squares                 | uint64    | `(nchans,)`           |
| `HIST`  | per-channel value histogram                | uint64    | `(nchans, 2^nbits)`   |
| `CMEA`  | per-chunk per-channel mean                 | float32   | `(nchunks, nchans)`   |
| `CRMS`  | per-chunk per-channel direct RMS           | float32   | `(nchunks, nchans)`   |
| `Z0DM`  | full-rate zero-DM time series (optional)   | float32   | `(nsamples,)`         |

`plot_diagnostics.py` produces 6 PNGs:

1. **`chan_summary.png`** — per-channel mean spectrum, min/max envelope
   over chunks, direct RMS, and 2-bit clipping fractions (counts at
   `v=0` and `v=2^nbits-1`).
2. **`waterfall.png`** — log-stretched time-frequency waterfall of
   `chunk_mean`, downsampled by `--fdec` and `--tdec`.
3. **`bandpass_stability.png`** — per-channel `chunk_rms / median(channel)`
   heatmap; a stable channel is uniform 1.0, drifts/RFI show up as
   stripes.
4. **`zerodm_05s.png`** — total power vs time at 0.5 s resolution.
5. **`zerodm_powerspectra.png`** — log-log overlay of `--n-spectra`
   evenly spaced FFTs of `--spectra-sec`-second slices of the zero-DM
   time series.
6. **`zerodm_sigma_exceedance.png`** — fraction of zero-DM samples above
   `k * rms` (positive, negative, two-sided), compared to the
   one-/two-sided Gaussian expectation `0.5 * erfc(k/sqrt(2))`. The
   intensity / chi-squared nature of the data shows up as deviations
   from Gaussian; transient events show up as flat tails far above the
   Gaussian curve.

## Notes on conventions

- Sample values are treated as **raw integer intensities** (`{0,1,2,3}`
  for 2-bit data). The data are positive-definite / chi-squared, so we
  use the **mean** and **direct RMS**, not median/MAD.
- Sigma exceedance is reported in units of the direct RMS, not MAD.
- The `Z0DM` section is large (`4 * nsamples` bytes; ~192 MB for 40 min
  at 50 us). Use `-Z` to drop it; the power-spectra and sigma-exceedance
  plots will then be skipped.

## Memory and throughput

At default settings on a 40-min/2-bit/2048-chan file:

| buffer                                                | size     |
| ----------------------------------------------------- | -------- |
| per-chunk packed + unpacked I/O buffers (`-c 0.5`)    | ~25 MB   |
| per-channel accumulators (`SUMC`, `SUMQ`, `HIST`)     | < 1 MB   |
| per-chunk arrays (`CMEA`, `CRMS`)                     | ~38 MB   |
| full-rate zero-DM (`Z0DM`)                            | ~192 MB  |
| **total RAM**                                         | ~260 MB  |

Throughput is dominated by sequential I/O from disk; on a fast SSD
expect a 25 GB file to be processed in a couple of minutes.

## Testing strategy

Layered, no real-data dependency:

- `tests/test_unpack.c` — exhaustive (256 single bytes) + 1 MiB random
  test of the 2-bit unpacker against the per-byte reference adapted
  from `dsa110-sigproc`'s `pack_unpack.c`. Plus `nbits=1,2,4,8`
  dispatch checks.
- `python/tools/make_test_fil.py` — generates short SIGPROC files with
  known content (`uniform`, `saturated_chan`, `transient`, `gaussian`).
- `tests/check_against_numpy.py` — runs `rfidiag` on each synthetic
  case and asserts the output matches an independent numpy
  computation **bit-exactly** for the integer accumulators
  (`SUMC`, `SUMQ`, `HIST`) and to better than 1e-5 relative for the
  float arrays (`CMEA`, `CRMS`, `Z0DM`).

Run all of the above with `make test`.

## Layout

```
simple_filtools/
├── Makefile
├── README.md
├── src/
│   ├── filhdr.{h,c}           # ported SIGPROC header reader (no globals)
│   ├── unpack.{h,c}           # n-bit unpack, fast 2-bit path
│   ├── diag.{h,c}             # tagged-section binary output writer
│   └── rfidiag.c              # streaming main
├── python/
│   ├── diagio.py              # .diag reader (numpy)
│   ├── plot_diagnostics.py    # 6 PNGs from a .diag
│   ├── requirements.txt
│   └── tools/
│       └── make_test_fil.py   # synthetic SIGPROC fil generator
└── tests/
    ├── test_unpack.c
    └── check_against_numpy.py
```
