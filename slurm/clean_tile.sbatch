#!/bin/bash
#SBATCH --job-name=mbclean
#SBATCH --output=mbclean_%j.out
#SBATCH --error=mbclean_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=128G
# Add --partition / --account / --qos as needed for your cluster, e.g.:
# #SBATCH --partition=expansion

# ---------------------------------------------------------------------------
# Multi-beam RFI clean for one tile of filterbank data.
#
# Usage:
#   sbatch clean_tile.sbatch                                   # uses defaults
#   sbatch clean_tile.sbatch /path/to/other/*.fil              # override inputs
#   sbatch --export=OUTDIR=/scratch/me/cleaned clean_tile.sbatch
#
# All knobs below can be overridden by exporting before sbatch, e.g.
#   CORES=32 CHUNK_SEC=0.5 sbatch clean_tile.sbatch
# ---------------------------------------------------------------------------

set -euo pipefail

INPUT_GLOB=${INPUT_GLOB:-"/resnick/groups/radiolab/fverdi/lowlat_survey/tile_17/test_filterbanks/*.fil"}
OUTDIR=${OUTDIR:-"/resnick/groups/radiolab/fverdi/lowlat_survey/tile_17/cleaned"}
DIAG_JSON=${DIAG_JSON:-"${OUTDIR}/mbclean_diag_${SLURM_JOB_ID:-local}.json"}

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
ROBUST_STATS=${ROBUST_STATS:-none}    # 'none' is fastest; use 'trimmed' for hard RFI
OUT_NBITS=${OUT_NBITS:-}              # leave empty to match input nbits

# ---------------------------------------------------------------------------
# Activate environment in non-interactive Slurm shells.
# Prefer micromamba (your setup), with conda fallback.
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
    echo "Set MAMBA_EXE/MAMBA_ROOT_PREFIX or ensure conda is installed." >&2
    exit 2
fi

# ---------------------------------------------------------------------------
# Pin BLAS / OpenMP thread pools to the number of cores Slurm allocated.
# multibeam_clean.py also exports these from --cores before importing numpy,
# but setting them here too is harmless and protects pre-import code paths.
# ---------------------------------------------------------------------------
export OMP_NUM_THREADS=$CORES
export MKL_NUM_THREADS=$CORES
export OPENBLAS_NUM_THREADS=$CORES
export NUMEXPR_NUM_THREADS=$CORES
export VECLIB_MAXIMUM_THREADS=$CORES

# Make outdir on a path Slurm can reach. Keep it adjacent to the input.
mkdir -p "$OUTDIR"

# ---------------------------------------------------------------------------
# Expand the glob NOW (inside this script, after conda is active), so the
# resulting list is what gets passed to python.
# ---------------------------------------------------------------------------
shopt -s nullglob
INPUTS=( $INPUT_GLOB )
if (( ${#INPUTS[@]} == 0 )); then
    echo "ERROR: no input files matched: $INPUT_GLOB" >&2
    exit 2
fi

echo "=========================================================================="
echo "host          : $(hostname)"
echo "job id        : ${SLURM_JOB_ID:-<no slurm>}"
echo "node          : ${SLURM_NODELIST:-<no slurm>}"
echo "cpus granted  : $CORES"
echo "code dir      : $CODE_DIR"
echo "conda env     : $CONDA_ENV"
echo "python        : $(which python)"
echo "n input files : ${#INPUTS[@]}"
echo "outdir        : $OUTDIR"
echo "diag json     : $DIAG_JSON"
echo "chunk-sec     : $CHUNK_SEC"
echo "min-support   : $MIN_SUPPORT  (RFI must touch >= this many beams)"
echo "robust-stats  : $ROBUST_STATS"
echo "=========================================================================="

# ---------------------------------------------------------------------------
# Run the cleaner.
# ---------------------------------------------------------------------------
OUTNBITS_ARG=()
if [[ -n "$OUT_NBITS" ]]; then
    OUTNBITS_ARG=(--out-nbits "$OUT_NBITS")
fi

srun --cpu-bind=cores python "$CODE_DIR/python/multibeam_clean.py" \
    "${INPUTS[@]}" \
    --outdir "$OUTDIR" \
    --cores "$CORES" \
    --chunk-sec "$CHUNK_SEC" \
    --tile-time-ms "$TILE_TIME_MS" \
    --tile-freq-chans "$TILE_FREQ_CHANS" \
    --freq-batch-tiles "$FREQ_BATCH_TILES" \
    --min-support "$MIN_SUPPORT" \
    --safety "$SAFETY" \
    --max-rank "$MAX_RANK" \
    --robust-stats "$ROBUST_STATS" \
    --diag "$DIAG_JSON" \
    "${OUTNBITS_ARG[@]}"

echo "----- mbclean done at $(date) -----"
