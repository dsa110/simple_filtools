#!/bin/bash
#SBATCH --job-name=mbclean_bench
#SBATCH --output=mbclean_bench_%j.out
#SBATCH --error=mbclean_bench_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G

# Short job that runs multibeam_clean.py in --benchmark mode (synthetic data,
# no IO) and prints the threading banner + per-stage timings. Use it to
# confirm BLAS / threadpoolctl / cgroup affinity are all happy before
# submitting a real job.

set -euo pipefail

CODE_DIR=${CODE_DIR:-"/home/vikram/installs/simple_filtools"}
CONDA_ENV=${CONDA_ENV:-"science"}
MAMBA_EXE=${MAMBA_EXE:-"/resnick/home/vikram/.local/bin/micromamba"}
MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-"/resnick/home/vikram/micromamba"}
CORES=${CORES:-${SLURM_CPUS_PER_TASK:-20}}
CHUNK_SEC=${CHUNK_SEC:-0.4}
TILE_TIME_MS=${TILE_TIME_MS:-100}
TILE_FREQ_CHANS=${TILE_FREQ_CHANS:-4}
FREQ_BATCH_TILES=${FREQ_BATCH_TILES:-8}
MIN_SUPPORT=${MIN_SUPPORT:-20}
SAFETY=${SAFETY:-1.10}
MAX_RANK=${MAX_RANK:-5}
ROBUST_STATS=${ROBUST_STATS:-none}
NO_PREFILTER=${NO_PREFILTER:-0}
# Bench needs at least one real input file (so it can read nchans/nbits/tsamp
# from the header). It WILL NOT touch the data; everything is synthetic.
INPUT_GLOB=${INPUT_GLOB:-"/resnick/groups/radiolab/fverdi/lowlat_survey/tile_17/test_filterbanks/*.fil"}

# Activate environment (micromamba first; conda fallback).
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
    echo "Set MAMBA_EXE/MAMBA_ROOT_PREFIX or ensure conda is installed." >&2
    exit 2
fi

export OMP_NUM_THREADS=$CORES
export MKL_NUM_THREADS=$CORES
export OPENBLAS_NUM_THREADS=$CORES

shopt -s nullglob
INPUTS=( $INPUT_GLOB )
if (( ${#INPUTS[@]} == 0 )); then
    echo "ERROR: no input files matched: $INPUT_GLOB" >&2
    exit 2
fi
echo "bench will read header from: ${INPUTS[0]} (+ ${#INPUTS[@]}-1 more for B)"
echo "cores=$CORES chunk_sec=$CHUNK_SEC tile_time_ms=$TILE_TIME_MS tile_freq=$TILE_FREQ_CHANS batch=$FREQ_BATCH_TILES robust=$ROBUST_STATS no_prefilter=$NO_PREFILTER"

PREFILTER_ARG=()
if [[ "$NO_PREFILTER" == "1" ]]; then
    PREFILTER_ARG=(--no-prefilter)
fi

srun --cpu-bind=cores python "$CODE_DIR/python/multibeam_clean.py" \
    --benchmark \
    --cores "$CORES" \
    --chunk-sec "$CHUNK_SEC" \
    --tile-time-ms "$TILE_TIME_MS" \
    --tile-freq-chans "$TILE_FREQ_CHANS" \
    --freq-batch-tiles "$FREQ_BATCH_TILES" \
    --min-support "$MIN_SUPPORT" \
    --safety "$SAFETY" \
    --max-rank "$MAX_RANK" \
    --robust-stats "$ROBUST_STATS" \
    "${PREFILTER_ARG[@]}" \
    "${INPUTS[@]}"
