#!/bin/bash
# SLURM array job — Ripley's K analysis for one (ligand, receptor, strip) per task.
#
# Reads the (1-indexed) row number from $SLURM_ARRAY_TASK_ID out of
# results/hpc_job_manifest.csv (skipping the header). Idempotent: skips
# tasks whose output parquet already exists.
#
# Submit full panel (438 jobs, max 5 concurrent):
#     sbatch scripts/run_array.sh
#
# Submit a small test (first 5 jobs only):
#     sbatch --array=1-5%2 scripts/run_array.sh
#
# Re-run only specific failing tasks:
#     sbatch --array=12,87,331 scripts/run_array.sh

#SBATCH --job-name=cosmx_arr
#SBATCH --partition=msc_appbio
#SBATCH --array=1-438%5
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%A/%a.out
#SBATCH --error=logs/%A/%a.err

set -euo pipefail

module load python/3.11.6-gcc-13.2.0
cd /scratch/users/$USER/NoSeggs/cosmx-pointpattern
source .venv/bin/activate

mkdir -p "results/per_pair" "logs/$SLURM_ARRAY_JOB_ID"

# Read row N+1 (skip header) from manifest
ROW=$(awk -v n=$SLURM_ARRAY_TASK_ID 'NR==n+1' results/hpc_job_manifest.csv)
if [[ -z "$ROW" ]]; then
    echo "ERROR: no row $SLURM_ARRAY_TASK_ID in manifest" >&2
    exit 2
fi

LIGAND=$(echo "$ROW"   | cut -d, -f1)
RECEPTOR=$(echo "$ROW" | cut -d, -f2)
STRIP=$(echo "$ROW"    | cut -d, -f3)

OUT="results/per_pair/${LIGAND}_${RECEPTOR}_${STRIP}.parquet"

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
    --n_sim   199 --seed 42 \
    --out     "$OUT"
