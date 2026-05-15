# simple_filtools

Standalone, dependency-free utilities for SIGPROC-format filterbank data.
Designed around 2-bit, 2048-channel, ~50 us, ~40-min DSA-2000 search-mode
files (~25 GB packed), but works on any SIGPROC `nbits in {1,2,4,8}` file
with `nifs == 1`.

Three C binaries:

| binary     | purpose                                                       |
| ---------- | ------------------------------------------------------------- |
| `rfidiag`  | streaming RFI diagnostics, writes a `*.diag` binary           |
| `chop_fil` | cut a `[start, start+dur]` slice into a new `.fil`            |
| `header`   | print SIGPROC header fields (mirrors sigproc's `header` tool) |

Two companion Python tools:

| script                          | purpose                                                              |
| ------------------------------- | -------------------------------------------------------------------- |
| `python/plot_diagnostics.py`    | turn a `.diag` into a set of PNG plots                               |
| `python/compare_beam_rfi.py`    | compare RFI signatures across multiple beams' `.diag` files          |
| `python/multibeam_clean.py`     | clean cross-beam-correlated RFI from N filterbank files in parallel  |

The SIGPROC header reader is ported from
[`dsa110-sigproc`](../dsa110-sigproc) (`read_header.c`, `pack_unpack.c`,
`strings_equal.c`) so this repo has no external C dependencies beyond
`libc + libm`.

---

## Build

```bash
make            # builds rfidiag, chop_fil, header
make test       # builds tests, runs unpack unit test + numpy cross-check
make clean
```

Tests need Python with `numpy` (matplotlib is required only for plotting).

## Run

### `rfidiag` — RFI diagnostics

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

### `chop_fil` — cut a time slice

```bash
./chop_fil obs.fil -s <start_sec> -d <dur_sec> [-o out.fil]
```

Writes a new SIGPROC filterbank file containing samples
`[round(start_sec/tsamp), +round(dur_sec/tsamp))` from `obs.fil`. The
header is copied verbatim from the input (preserving any non-standard
tags) with a single edit: `tstart` is incremented by
`start_samp * tsamp / 86400` days. Default output is `<input>.cut.fil`.

If `dur_sec` would run past the end of the file, the cut is silently
truncated and a warning is printed to stderr. If `start_sec` is past the
end of the file, `chop_fil` exits with an error.

### `header` — print header fields

```bash
./header obs.fil                     # full human-readable summary
./header obs.fil -tstart -tsamp      # one value per line, for shell use
```

Standalone reimplementation of dsa110-sigproc's `header` tool on top of
the simple_filtools header reader. Supported flags (`./header --help`):

```
-source_name -telescope -machine -datatype -data_type
-frame -barycentric -pulsarcentric
-headersize -datasize -nsamples -tobs
-fch1 -foff -bandwidth -fmid -nchans -nbits -nifs -nbeam -ibeam
-tstart -tsamp -mjd -utstart -date
-src_raj -src_dej -ra_deg -dec_deg -az_start -za_start
-refdm (alias -dm) -frequencies
```

Differences from sigproc's tool: WAPP/PSPM/BPP raw formats are not
supported; `-frequencies` always prints `fch1 + i*foff` rather than
reading a `FREQUENCY_START`/`FREQUENCY_END` channel table; `-obsdb`,
`-scan_number`, and the `signed` 8-bit flag are not implemented.

### `multibeam_clean.py` — remove RFI correlated across beam groups

```bash
python3 python/multibeam_clean.py \
    beam00.fil beam01.fil ... beam71.fil \
    --outdir cleaned/ \
    --cores 20 \
    --min-support 20 \
    --chunk-sec 1.0 \
    --tile-time-ms 100 \
    --tile-freq-chans 4 \
    --diag run.jsonl
```

Streams an arbitrary set of SIGPROC filterbank files (one per beam,
any duration), reads all data-layout parameters (`nchans`, `nbits`,
`fch1`, `foff`, `tsamp`, `tstart`, ...) from the headers, and writes
one cleaned `.fil` per input into `--outdir`. Within each time chunk:

1. Per-(beam, channel) whitening using a robust running scale (MAD,
   EMA-smoothed across chunks with timescale `--whitening-tau-sec`).
2. For each `(tile_time_ms, tile_freq_chans)` tile, form the `BxB`
   cross-beam covariance and diagonalize it.
3. Treat eigenpair `(lambda_k, u_k)` as RFI iff
   * `lambda_k > MP_edge(B, N) * --safety`  (Marchenko-Pastur null),
   * `participation_ratio(u_k) >= --min-support`  (>= N beams), and
   * fewer than `--max-rank` components have already been removed.
4. Subtract that subspace from the data and de-whiten.
5. Re-quantize back to the input `nbits` and append to each output.

Pass `--dry-run` to validate headers and print the processing plan
without doing any work. Pass `--diag run.jsonl` for per-chunk
JSON-lines diagnostics (`max_eigval`, `n_tiles_with_rfi`, ...).

All knobs (with sensible defaults; see `--help`):

| flag                     | meaning                                                   |
| ------------------------ | --------------------------------------------------------- |
| `--cores N`              | BLAS + parallel-IO threads (sets `OMP_NUM_THREADS` etc.)  |
| `--chunk-sec`            | streaming chunk length                                    |
| `--whitening-tau-sec`    | EMA timescale for running mean / scale                    |
| `--tile-time-ms`         | time tile size for eigen-clean                            |
| `--tile-freq-chans`      | freq tile size in channels                                |
| `--min-support`          | participation-ratio threshold ("only flag RFI that hits >= N beams") |
| `--safety`               | multiplicative safety on MP edge                          |
| `--max-rank`             | cap on # eigenvectors subtracted per tile                 |
| `--out-nbits`            | output bit depth (default: input)                         |
| `--tstart-tol-sec`       | tolerance on tstart alignment across beams                |
| `--diag`                 | per-chunk JSON-lines log                                  |
| `--dry-run`              | parse + validate + print plan, no IO                      |

Smoke test:

```bash
python3 python/tools/test_multibeam_clean.py
```

Plants a common-mode burst on 20 of 24 synthetic beams and verifies
both the end-to-end 2-bit roundtrip and the eigen-clean kernel
suppress the planted signal while leaving uninvolved beams untouched.

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
│   ├── rfidiag.c              # streaming RFI diagnostics main
│   ├── chop_fil.c             # cut a time slice into a new .fil
│   └── header.c               # print SIGPROC header fields
├── python/
│   ├── diagio.py              # .diag reader (numpy)
│   ├── plot_diagnostics.py    # 6 PNGs from a .diag
│   ├── compare_beam_rfi.py    # multi-beam .diag comparison + spatial-filter feasibility
│   ├── multibeam_clean.py     # cross-beam-correlated RFI cleaner (eigen subtraction)
│   ├── requirements.txt
│   └── tools/
│       ├── make_test_fil.py            # synthetic SIGPROC fil generator
│       └── test_multibeam_clean.py     # synthetic smoke test for multibeam_clean.py
└── tests/
    ├── test_unpack.c
    └── check_against_numpy.py
```
