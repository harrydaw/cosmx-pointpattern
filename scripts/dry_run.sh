#!/bin/bash
# Single-job HPC dry run — KRT8 x KRT18 on strip_2, n_sim=9.
# Cheap smoke test before submitting the full array.
#
# Submit with:  sbatch scripts/dry_run.sh
# Inspect with: tail -f logs/dryrun_<jobid>.out
#               sacct -j <jobid> --format=JobID,State,ExitCode,MaxRSS,Elapsed

#SBATCH --job-name=cosmx_dryrun
#SBATCH --partition=msc_appbio
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/dryrun_%j.out
#SBATCH --error=logs/dryrun_%j.err

set -euo pipefail

module load python/3.11.6-gcc-13.2.0
cd /scratch/users/$USER/NoSeggs/cosmx-pointpattern
source .venv/bin/activate

mkdir -p results/per_pair logs

python scripts/batch_k_analysis.py \
    --data    data/processed/s1_all_strips_noise_flagged.parquet \
    --fns     notebooks/00_functions.ipynb \
    --r_vals  data/processed/r_vals.npy \
    --pair    "KRT8,KRT18" \
    --strip   strip_2 \
    --window  concave --concave_ratio 0.1 \
    --n_sim   9 --seed 42 \
    --out     results/per_pair/KRT8_KRT18_strip_2.parquet
