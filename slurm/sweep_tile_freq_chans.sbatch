#!/bin/bash
#SBATCH --job-name=mbclean_sweep_f
#SBATCH --output=mbclean_sweep_f_%j.out
#SBATCH --error=mbclean_sweep_f_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G

# Sweep multibeam_clean over a range of --tile-freq-chans values.
#
# Default use:
#   sbatch sweep_tile_freq_chans.sbatch
#
# Override examples:
#   TILE_FREQ_LIST="1 2 4 8 16 32" sbatch sweep_tile_freq_chans.sbatch
#   OUT_BASE=/resnick/groups/radiolab/vikram/scratch/clean_sweep sbatch sweep_tile_freq_chans.sbatch
#   CHUNK_SEC=0.4 MIN_SUPPORT=10 SAFETY=1.03 MAX_RANK=8 sbatch sweep_tile_freq_chans.sbatch

set -euo pipefail

CODE_DIR=${CODE_DIR:-"/home/vikram/installs/simple_filtools"}
CONDA_ENV=${CONDA_ENV:-"science"}
MAMBA_EXE=${MAMBA_EXE:-"/resnick/home/vikram/.local/bin/micromamba"}
MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-"/resnick/home/vikram/micromamba"}

INPUT_GLOB=${INPUT_GLOB:-"/resnick/groups/radiolab/fverdi/lowlat_survey/tile_17/test_filterbanks/*.fil"}
OUT_BASE=${OUT_BASE:-"/resnick/groups/radiolab/vikram/scratch/clean_sweep_${SLURM_JOB_ID:-local}"}

CORES=${CORES:-${SLURM_CPUS_PER_TASK:-20}}
CHUNK_SEC=${CHUNK_SEC:-0.4}
TILE_TIME_MS=${TILE_TIME_MS:-25}
TILE_FREQ_LIST=${TILE_FREQ_LIST:-"1 2 4 8 16 32"}
FREQ_BATCH_TILES=${FREQ_BATCH_TILES:-32}
MIN_SUPPORT=${MIN_SUPPORT:-10}
SAFETY=${SAFETY:-1.03}
MAX_RANK=${MAX_RANK:-8}
ROBUST_STATS=${ROBUST_STATS:-none}
NO_PREFILTER=${NO_PREFILTER:-0}
OUT_NBITS=${OUT_NBITS:-}

# ---------------------------------------------------------------------------
# Activate environment (micromamba first, conda fallback).
# ---------------------------------------------------------------------------
if [[ -x "$MAMBA_EXE" ]]; then
    __mamba_setup="$("$MAMBA_EXE" shell hook --shell bash --root-prefix "$MAMBA_ROOT_PREFIX" 2>/dev/null)"
    if [[ $? -eq 0 ]]; then
        eval "$__mamba_setup"
    else
        alias micromamba="$MAMBA_EXE"
    fi
    unset __mamba_setup
    micromamba activate "$CONDA_ENV"
elif command -v conda >/dev/null 2>&1 || [[ -n "${CONDA_EXE:-}" ]]; then
    if [[ -z "${CONDA_EXE:-}" ]]; then
        for candidate in \
            "$HOME/miniconda3/etc/profile.d/conda.sh" \
            "$HOME/anaconda3/etc/profile.d/conda.sh" \
            "/opt/conda/etc/profile.d/conda.sh"; do
            if [[ -f "$candidate" ]]; then
                # shellcheck disable=SC1090
                source "$candidate"
                break
            fi
        done
    else
        # shellcheck disable=SC1091
        source "$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh"
    fi
    conda activate "$CONDA_ENV"
else
    echo "ERROR: neither micromamba nor conda was found on this node." >&2
    exit 2
fi

export OMP_NUM_THREADS=$CORES
export MKL_NUM_THREADS=$CORES
export OPENBLAS_NUM_THREADS=$CORES
export NUMEXPR_NUM_THREADS=$CORES

shopt -s nullglob
INPUTS=( $INPUT_GLOB )
if (( ${#INPUTS[@]} == 0 )); then
    echo "ERROR: no inputs matched INPUT_GLOB=$INPUT_GLOB" >&2
    exit 2
fi

mkdir -p "$OUT_BASE"

echo "=========================================================================="
echo "host             : $(hostname)"
echo "job id           : ${SLURM_JOB_ID:-<no slurm>}"
echo "cpus granted     : $CORES"
echo "python           : $(which python)"
echo "input files      : ${#INPUTS[@]}"
echo "input glob       : $INPUT_GLOB"
echo "output base      : $OUT_BASE"
echo "tile_time_ms     : $TILE_TIME_MS"
echo "tile_freq_list   : $TILE_FREQ_LIST"
echo "chunk_sec        : $CHUNK_SEC"
echo "min_support      : $MIN_SUPPORT"
echo "safety           : $SAFETY"
echo "max_rank         : $MAX_RANK"
echo "robust_stats     : $ROBUST_STATS"
echo "prefilter        : $([[ \"$NO_PREFILTER\" == \"1\" ]] && echo off || echo on)"
echo "=========================================================================="

PREFILTER_ARG=()
if [[ "$NO_PREFILTER" == "1" ]]; then
    PREFILTER_ARG=(--no-prefilter)
fi

OUTNBITS_ARG=()
if [[ -n "$OUT_NBITS" ]]; then
    OUTNBITS_ARG=(--out-nbits "$OUT_NBITS")
fi

SUMMARY="$OUT_BASE/sweep_summary.tsv"
echo -e "tile_freq_chans\tstatus\telapsed_min\trt_factor\ttiles_with_rfi_pct\tdiag_json" > "$SUMMARY"

for FCH in $TILE_FREQ_LIST; do
    RUN_DIR="$OUT_BASE/fch_${FCH}"
    OUTDIR="$RUN_DIR/cleaned"
    DIAG_JSON="$RUN_DIR/mbclean_diag.jsonl"
    RUN_LOG="$RUN_DIR/run.log"
    mkdir -p "$RUN_DIR" "$OUTDIR"

    echo
    echo "----- sweep run: tile-freq-chans=$FCH -----"
    t0=$(date +%s)

    set +e
    srun --cpu-bind=cores python "$CODE_DIR/python/multibeam_clean.py" \
        "${INPUTS[@]}" \
        --outdir "$OUTDIR" \
        --cores "$CORES" \
        --chunk-sec "$CHUNK_SEC" \
        --tile-time-ms "$TILE_TIME_MS" \
        --tile-freq-chans "$FCH" \
        --freq-batch-tiles "$FREQ_BATCH_TILES" \
        --min-support "$MIN_SUPPORT" \
        --safety "$SAFETY" \
        --max-rank "$MAX_RANK" \
        --robust-stats "$ROBUST_STATS" \
        --diag "$DIAG_JSON" \
        "${PREFILTER_ARG[@]}" \
        "${OUTNBITS_ARG[@]}" \
        2>&1 | tee "$RUN_LOG"
    rc=${PIPESTATUS[0]}
    set -e

    t1=$(date +%s)
    elapsed_min=$(python - <<PY
print(f"{(${t1}-${t0})/60.0:.2f}")
PY
)

    rt="NA"
    pct="NA"
    if [[ -f "$RUN_LOG" ]]; then
        rt_line=$(grep -E "rt-factor" "$RUN_LOG" | tail -n 1 || true)
        pct_line=$(grep -E "tiles_with_rfi" "$RUN_LOG" | tail -n 1 || true)
        if [[ -n "$rt_line" ]]; then
            rt=$(echo "$rt_line" | sed -E 's/.*rt-factor ([0-9.]+)x.*/\1/' || echo "NA")
        fi
        if [[ -n "$pct_line" ]]; then
            pct=$(echo "$pct_line" | sed -E 's/.*\(([0-9.]+)%\).*/\1/' || echo "NA")
        fi
    fi

    if [[ $rc -eq 0 ]]; then
        status="ok"
    else
        status="fail($rc)"
    fi

    echo -e "${FCH}\t${status}\t${elapsed_min}\t${rt}\t${pct}\t${DIAG_JSON}" >> "$SUMMARY"
done

echo
echo "----- sweep complete -----"
echo "summary: $SUMMARY"
echo
cat "$SUMMARY"

