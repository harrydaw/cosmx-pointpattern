#!/bin/bash
# SLURM array job — extended L-R panel (v2) at n_sim=99.
# Companion to scripts/run_array.sh — v1 was the 438-job CellChatDB run at n_sim=199.
# This v2 runs the ~80-pair pre-registered extended panel from
#   notes/extended_panel_rationale.md
# at n_sim=99 (justified by A.6: envelope is indistinguishable at 99 vs 199
# under the production window; running at 99 is purely a speed choice).
#
# Pre-flight (locally or on HPC):
#   python scripts/build_v2_manifest.py
# This prints the exact sbatch command, e.g.:
#   sbatch --array=1-215%2 scripts/run_array_v2.sh
#
# Re-run only specific failing tasks:
#   sbatch --array=12,87,118 scripts/run_array_v2.sh

#SBATCH --job-name=cosmx_v2
#SBATCH --partition=msc_appbio
#SBATCH --array=1-243%2
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs_v2/%A/%a.out
#SBATCH --error=logs_v2/%A/%a.err

set -euo pipefail

module load python/3.11.6-gcc-13.2.0
cd /scratch/users/$USER/NoSeggs/cosmx-pointpattern
source .venv/bin/activate

mkdir -p "results/per_pair_v2" "logs_v2/$SLURM_ARRAY_JOB_ID"

# Read row N+1 (skip header) from v2 manifest
ROW=$(awk -v n=$SLURM_ARRAY_TASK_ID 'NR==n+1' results/hpc_job_manifest_v2.csv)
if [[ -z "$ROW" ]]; then
    echo "ERROR: no row $SLURM_ARRAY_TASK_ID in v2 manifest" >&2
    exit 2
fi

LIGAND=$(echo "$ROW"   | cut -d, -f1)
RECEPTOR=$(echo "$ROW" | cut -d, -f2)
STRIP=$(echo "$ROW"    | cut -d, -f3)

OUT="results/per_pair_v2/${LIGAND}_${RECEPTOR}_${STRIP}.parquet"

if [[ -f "$OUT" ]]; then
    echo "SKIP: $OUT already exists"
    exit 0
fi

echo "Task $SLURM_ARRAY_TASK_ID: $LIGAND x $RECEPTOR | $STRIP"

python scripts/batch_k_analysis.py \
    --data    data/processed/s1_all_strips_cleaned.parquet \
    --fns     notebooks/00_functions.ipynb \
    --r_vals  data/processed/r_vals.npy \
    --pair    "${LIGAND},${RECEPTOR}" \
    --strip   "$STRIP" \
    --window  concave --concave_ratio 0.1 \
    --n_sim   99 --seed 42 \
    --out     "$OUT"
