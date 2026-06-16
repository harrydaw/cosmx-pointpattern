#!/bin/bash
# Reproduce the control-pair analyses with the EXACT production pipeline
# (concave hull ratio 0.1 + Shapely edge correction, n_sim=199, seed=42), on all 3 strips.
# Mirrors scripts/run_array.sh but for the named control pairs, with output kept
# separate from the canonical panel parquets.
#
# Run from the repo root with the venv active:
#     source .venv/bin/activate          # Windows: source .venv/Scripts/activate
#     bash scripts/run_controls.sh
#
# Controls:
#   Positive    KRT8 x KRT18    -> strongest co-localisation (epithelial keratins)
#   Negative 1  MALAT1 x KRT18  -> ubiquitous x epithelial (no gene-specific co-loc)
#   Negative 2  KRT8 x SCGB3A1  -> epithelial vs secretory airway (distinct compartments)
#
# Outputs: results/controls/<L>_<R>_<strip>.parquet  (+ .r_vals.npy)
# Verdict is read from l_obs vs [l_lo, l_hi]; see Final_Writeup/METHODS_LOG.md.

set -euo pipefail

mkdir -p results/controls

PAIRS=("KRT8,KRT18" "MALAT1,KRT18" "KRT8,SCGB3A1")

for PAIR in "${PAIRS[@]}"; do
    L=${PAIR%,*}
    R=${PAIR#*,}
    for S in strip_1 strip_2 strip_3; do
        OUT="results/controls/${L}_${R}_${S}.parquet"
        echo "=== ${L} x ${R} | ${S} ==="
        python scripts/batch_k_analysis.py \
            --data    data/processed/s1_all_strips_cleaned.parquet \
            --fns     notebooks/00_functions.ipynb \
            --r_vals  data/processed/r_vals.npy \
            --pair    "${L},${R}" \
            --strip   "${S}" \
            --window  concave --concave_ratio 0.1 \
            --n_sim   199 --seed 42 \
            --out     "${OUT}"
    done
done

echo "All control jobs complete -> results/controls/"
