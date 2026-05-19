#!/bin/bash
#SBATCH --job-name=mbdiag_compare
#SBATCH --output=mbdiag_compare_%j.out
#SBATCH --error=mbdiag_compare_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# Build before/after diagnostics from dirty+clean beam sets.
#
# Default paths:
#   dirty: /resnick/groups/radiolab/fverdi/lowlat_survey/tile_17/test_filterbanks/*.fil
#   clean: /resnick/groups/radiolab/vikram/scratch/cleaned/*.fil
#
# Submit:
#   sbatch diag_compare_dirty_clean.sbatch
#
# Override knobs:
#   DIRTY_GLOB='/path/dirty/*.fil' CLEAN_GLOB='/path/clean/*.fil' \
#   OUT_ROOT='/path/out' NSIGMA=5 TOPK=20 sbatch diag_compare_dirty_clean.sbatch

set -euo pipefail

CODE_DIR=${CODE_DIR:-"/home/vikram/installs/simple_filtools"}
CONDA_ENV=${CONDA_ENV:-"science"}
MAMBA_EXE=${MAMBA_EXE:-"/resnick/home/vikram/.local/bin/micromamba"}
MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-"/resnick/home/vikram/micromamba"}

DIRTY_GLOB=${DIRTY_GLOB:-"/resnick/groups/radiolab/fverdi/lowlat_survey/tile_17/test_filterbanks/*.fil"}
CLEAN_GLOB=${CLEAN_GLOB:-"/resnick/groups/radiolab/vikram/scratch/cleaned/*.fil"}
OUT_ROOT=${OUT_ROOT:-"/resnick/groups/radiolab/vikram/scratch/diag_compare_${SLURM_JOB_ID:-local}"}

CHUNK_SEC=${CHUNK_SEC:-0.5}
NSIGMA=${NSIGMA:-5}
TOPK=${TOPK:-20}
ZDM_BIN_SEC=${ZDM_BIN_SEC:-1e-3}

CORES=${CORES:-${SLURM_CPUS_PER_TASK:-8}}

# ---------------------------------------------------------------------------
# Activate environment in non-interactive Slurm shells.
# Prefer micromamba, with conda fallback.
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

export OMP_NUM_THREADS=$CORES
export MKL_NUM_THREADS=$CORES
export OPENBLAS_NUM_THREADS=$CORES
export NUMEXPR_NUM_THREADS=$CORES

shopt -s nullglob
DIRTY_FILS=( $DIRTY_GLOB )
CLEAN_FILS=( $CLEAN_GLOB )
if (( ${#DIRTY_FILS[@]} == 0 )); then
    echo "ERROR: no dirty files matched DIRTY_GLOB=$DIRTY_GLOB" >&2
    exit 2
fi
if (( ${#CLEAN_FILS[@]} == 0 )); then
    echo "ERROR: no clean files matched CLEAN_GLOB=$CLEAN_GLOB" >&2
    exit 2
fi

DIRTY_DIAG_DIR="$OUT_ROOT/dirty_diag"
CLEAN_DIAG_DIR="$OUT_ROOT/clean_diag"
COMPARE_DIRTY="$OUT_ROOT/compare_dirty"
COMPARE_CLEAN="$OUT_ROOT/compare_clean"
mkdir -p "$DIRTY_DIAG_DIR" "$CLEAN_DIAG_DIR" "$COMPARE_DIRTY" "$COMPARE_CLEAN"

echo "=========================================================================="
echo "host          : $(hostname)"
echo "job id        : ${SLURM_JOB_ID:-<no slurm>}"
echo "node          : ${SLURM_NODELIST:-<no slurm>}"
echo "cpus granted  : $CORES"
echo "python        : $(which python)"
echo "code dir      : $CODE_DIR"
echo "conda env     : $CONDA_ENV"
echo "dirty files   : ${#DIRTY_FILS[@]} (glob: $DIRTY_GLOB)"
echo "clean files   : ${#CLEAN_FILS[@]} (glob: $CLEAN_GLOB)"
echo "out root      : $OUT_ROOT"
echo "rfidiag chunk : $CHUNK_SEC sec"
echo "compare nsigma: $NSIGMA ; topk: $TOPK ; zdm-bin: $ZDM_BIN_SEC sec"
echo "=========================================================================="

# ---------------------------------------------------------------------------
# 1) Build .diag for dirty beams
# ---------------------------------------------------------------------------
echo "--- building dirty .diag files ---"
for f in "${DIRTY_FILS[@]}"; do
    b=$(basename "$f" .fil)
    out="$DIRTY_DIAG_DIR/${b}.diag"
    echo "dirty: $b"
    srun --cpu-bind=cores --ntasks=1 --cpus-per-task=1 \
        "$CODE_DIR/rfidiag" "$f" -o "$out" -c "$CHUNK_SEC"
done

# ---------------------------------------------------------------------------
# 2) Build .diag for clean beams
# ---------------------------------------------------------------------------
echo "--- building clean .diag files ---"
for f in "${CLEAN_FILS[@]}"; do
    b=$(basename "$f" .fil)
    out="$CLEAN_DIAG_DIR/${b}.diag"
    echo "clean: $b"
    srun --cpu-bind=cores --ntasks=1 --cpus-per-task=1 \
        "$CODE_DIR/rfidiag" "$f" -o "$out" -c "$CHUNK_SEC"
done

# ---------------------------------------------------------------------------
# 3) Compare dirty beams against each other
# ---------------------------------------------------------------------------
echo "--- compare dirty set ---"
srun --cpu-bind=cores python "$CODE_DIR/python/compare_beam_rfi.py" \
    "$DIRTY_DIAG_DIR"/*.diag \
    --outdir "$COMPARE_DIRTY" \
    --nsigma "$NSIGMA" \
    --top-k "$TOPK" \
    --zdm-bin-sec "$ZDM_BIN_SEC"

# ---------------------------------------------------------------------------
# 4) Compare clean beams against each other
# ---------------------------------------------------------------------------
echo "--- compare clean set ---"
srun --cpu-bind=cores python "$CODE_DIR/python/compare_beam_rfi.py" \
    "$CLEAN_DIAG_DIR"/*.diag \
    --outdir "$COMPARE_CLEAN" \
    --nsigma "$NSIGMA" \
    --top-k "$TOPK" \
    --zdm-bin-sec "$ZDM_BIN_SEC"

# ---------------------------------------------------------------------------
# 5) Emit a compact before/after summary from per_chan.csv
# ---------------------------------------------------------------------------
echo "--- writing before/after summary ---"
export DIRTY_CSV="$COMPARE_DIRTY/per_chan.csv"
export CLEAN_CSV="$COMPARE_CLEAN/per_chan.csv"
export SUMMARY_TXT="$OUT_ROOT/before_after_summary.txt"
export SUMMARY_CSV="$OUT_ROOT/before_after_summary.csv"

python - <<'PY'
import csv
import os
import numpy as np

dirty_csv = os.environ["DIRTY_CSV"]
clean_csv = os.environ["CLEAN_CSV"]
summary_txt = os.environ["SUMMARY_TXT"]
summary_csv = os.environ["SUMMARY_CSV"]

def load_counts(path):
    vals = []
    with open(path) as f:
        r = csv.DictReader(f)
        for row in r:
            vals.append(int(row["n_beams_flagged"]))
    return np.asarray(vals, dtype=np.int64)

d = load_counts(dirty_csv)
c = load_counts(clean_csv)
n = min(d.size, c.size)
d = d[:n]
c = c[:n]

metrics = [
    ("n_channels_compared", n, n, 0),
    ("any_beam_flagged", int((d > 0).sum()), int((c > 0).sum()), int((c > 0).sum() - (d > 0).sum())),
    ("all_beams_flagged", int((d == d.max()).sum()), int((c == c.max()).sum()), int((c == c.max()).sum() - (d == d.max()).sum())),
    ("mean_n_beams_flagged", float(d.mean()), float(c.mean()), float(c.mean() - d.mean())),
    ("p95_n_beams_flagged", float(np.percentile(d, 95)), float(np.percentile(c, 95)), float(np.percentile(c, 95) - np.percentile(d, 95))),
]

with open(summary_txt, "w") as f:
    f.write("before/after summary (dirty -> clean)\n")
    f.write(f"dirty_csv={dirty_csv}\n")
    f.write(f"clean_csv={clean_csv}\n\n")
    for name, vd, vc, delta in metrics:
        f.write(f"{name:24s} dirty={vd} clean={vc} delta={delta}\n")

with open(summary_csv, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["metric", "dirty", "clean", "delta_clean_minus_dirty"])
    for row in metrics:
        w.writerow(row)

print(f"wrote {summary_txt}")
print(f"wrote {summary_csv}")
PY

echo "----- done -----"
echo "outputs:"
echo "  $OUT_ROOT"
echo "  $COMPARE_DIRTY"
echo "  $COMPARE_CLEAN"
echo "  $SUMMARY_TXT"
echo "  $SUMMARY_CSV"
