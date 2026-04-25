#!/usr/bin/env python
"""
Batch Ripley's K-function analysis for all viable LR pairs.

Loads all function definitions from 00_functions.ipynb, then runs
run_pair_analysis() serially across every viable (ligand, receptor, strip)
combination found in the LR panel index.

Usage
-----
On HPC (after transferring the cosmx-pointpattern directory):

    python scripts/batch_k_analysis.py \\
        --data   data/processed/s1_all_strips_cleaned.parquet \\
        --lr     data/processed/lr_panel_index.json \\
        --fns    notebooks/00_functions.ipynb \\
        --out    results/batch_results.parquet

Required flags:
    --r_vals        path to .npy file containing the r_vals array
                    (written alongside the parquet by the notebook)

Optional flags:
    --window        concave | hull | rect | custom   (default: concave)
    --concave_ratio float                            (default: 0.1, only used when --window=concave)
    --n_sim         int                              (default: 99)
    --seed          int                              (default: 42)
    --strips        strip_1,strip_2,strip_3  comma-separated (default: all three)

Results parquet schema
----------------------
ligand, receptor, strip, n_a, n_b, l_obs (list), l_lo (list), l_hi (list),
window_type, window_area_px2, status

r_vals is serialised next to the output parquet as `<out>.r_vals.npy`
so downstream code can load it with np.load(path).

Requirements
------------
pandas, numpy, scipy, shapely, scikit-learn, pyarrow, nbformat
"""

import argparse
import json
import os
import sys
import time

import nbformat
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Load all function definitions from 00_functions.ipynb at runtime
# ---------------------------------------------------------------------------

def load_functions_from_notebook(nb_path: str) -> None:
    """
    Execute every code cell in a Jupyter notebook, injecting all definitions
    into the calling module's global namespace.

    This lets the batch script use functions defined in 00_functions.ipynb
    without a separate pip-installable package.
    """
    nb = nbformat.read(nb_path, as_version=4)
    g = globals()
    for cell in nb.cells:
        if cell.cell_type == "code":
            try:
                exec(compile(cell.source, nb_path, "exec"), g)
            except Exception as e:
                print(f"  [warn] cell exec failed ({e}): {cell.source[:60]!r}")


# ---------------------------------------------------------------------------
# Core batch logic
# ---------------------------------------------------------------------------

def build_job_list(lr_index: dict, genes_in_data: set, strips: list) -> list:
    """Return list of (ligand, receptor, strip) tuples for all viable jobs."""
    jobs = []
    for key, strip_data in lr_index.items():
        lig, rec = key.split("|")
        if lig not in genes_in_data or rec not in genes_in_data:
            continue
        for strip in strips:
            if strip_data.get(strip, {}).get("viable", False):
                jobs.append((lig, rec, strip))
    return jobs


def run_batch(
    clean_df: pd.DataFrame,
    jobs: list,
    r_vals: np.ndarray,
    window_type: str,
    concave_ratio: float,
    n_sim: int,
    seed: int,
    resolution: int,
) -> list:
    """
    Run run_pair_analysis() for every job serially.

    Returns a list of result dicts suitable for pd.DataFrame().
    """
    STRIPS = ["strip_1", "strip_2", "strip_3"]
    strip_dfs = {s: clean_df[clean_df["strip"] == s] for s in STRIPS}

    results = []
    t_start = time.time()

    for i, (lig, rec, strip) in enumerate(jobs):
        t0 = time.time()
        print(f"[{i+1}/{len(jobs)}] {lig} × {rec} | {strip} ", end="", flush=True)
        try:
            res = run_pair_analysis(          # noqa: F821  (injected at runtime)
                strip_dfs[strip],
                lig,
                rec,
                r_vals,
                window_type=window_type,
                concave_ratio=concave_ratio,
                n_sim=n_sim,
                seed=seed,
                resolution=resolution,
            )
            window_area = (
                res["window"].area
                if not isinstance(res["window"], tuple)
                else (
                    (res["window"][1] - res["window"][0])
                    * (res["window"][3] - res["window"][2])
                )
            )
            results.append(
                {
                    "ligand":          lig,
                    "receptor":        rec,
                    "strip":           strip,
                    "n_a":             res["n_a"],
                    "n_b":             res["n_b"],
                    "l_obs":           res["l_obs"].tolist(),
                    "l_lo":            res["l_lo"].tolist(),
                    "l_hi":            res["l_hi"].tolist(),
                    "window_type":     window_type,
                    "window_area_px2": round(window_area),
                    "status":          "ok",
                }
            )
            elapsed = time.time() - t0
            print(f"ok ({elapsed:.0f}s)")

        except Exception as e:
            print(f"ERROR: {e}")
            results.append(
                {
                    "ligand":    lig,
                    "receptor":  rec,
                    "strip":     strip,
                    "status":    f"error: {e}",
                }
            )

        # Running ETA
        done      = i + 1
        total_t   = time.time() - t_start
        rate      = done / total_t
        remaining = len(jobs) - done
        eta_s     = remaining / rate if rate > 0 else 0
        print(
            f"    Progress: {done}/{len(jobs)} | "
            f"elapsed {total_t/60:.1f} min | ETA {eta_s/60:.1f} min"
        )

    return results


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch Ripley's K-function analysis for all viable LR pairs."
    )
    parser.add_argument(
        "--data", required=True,
        help="Path to s1_all_strips_cleaned.parquet",
    )
    parser.add_argument(
        "--lr", required=True,
        help="Path to lr_panel_index.json",
    )
    parser.add_argument(
        "--fns", required=True,
        help="Path to 00_functions.ipynb",
    )
    parser.add_argument(
        "--out", required=True,
        help="Output parquet path for results",
    )
    parser.add_argument(
        "--window", default="concave",
        choices=["concave", "hull", "rect", "custom"],
        help="Window type (default: concave)",
    )
    parser.add_argument(
        "--concave_ratio", type=float, default=0.1,
        help="Concave hull ratio (0=tightest, 1=convex). Only used when --window=concave.",
    )
    parser.add_argument(
        "--r_vals", required=True,
        help="Path to .npy file containing r_vals array (written by the notebook).",
    )
    parser.add_argument("--n_sim",  type=int,   default=99)
    parser.add_argument("--seed",   type=int,   default=42)
    parser.add_argument("--resolution", type=int, default=64)
    parser.add_argument(
        "--strips",
        default="strip_1,strip_2,strip_3",
        help="Comma-separated list of strips to include (default: all three)",
    )
    args = parser.parse_args()

    strips = [s.strip() for s in args.strips.split(",")]

    # -- Load functions -------------------------------------------------------
    print("=" * 60)
    print("Loading functions from:", args.fns)
    load_functions_from_notebook(args.fns)
    print("Functions loaded.")

    # -- Load data ------------------------------------------------------------
    print("\nLoading transcript data from:", args.data)
    df = pd.read_parquet(args.data)
    clean = df[
        ~df["is_noise"] &
        ~df["is_small_cluster"] &
        ~df["manually_excluded"]
    ].copy()
    genes_in_data = set(clean["target"].unique())
    print(f"  Clean transcripts: {len(clean):,}  |  Unique genes: {len(genes_in_data)}")

    # -- Load LR index --------------------------------------------------------
    print("\nLoading LR index from:", args.lr)
    with open(args.lr) as f:
        lr_index = json.load(f)
    print(f"  Pairs in index: {len(lr_index)}")

    # -- Load r_vals ----------------------------------------------------------
    print("\nLoading r_vals from:", args.r_vals)
    r_vals = np.load(args.r_vals)
    if r_vals.ndim != 1 or len(r_vals) < 2:
        raise ValueError(f"r_vals must be 1-D with >=2 entries, got shape {r_vals.shape}")

    # -- Build job list -------------------------------------------------------
    jobs = build_job_list(lr_index, genes_in_data, strips)
    print(f"\nJobs to run: {len(jobs)}")
    print(f"  r_vals: {len(r_vals)} points, {r_vals[0]:.1f} -> {r_vals[-1]:.1f} px")
    print(f"  n_sim:  {args.n_sim}")
    print(f"  window: {args.window}" + (f" (ratio={args.concave_ratio})" if args.window == "concave" else ""))
    print(f"  strips: {strips}")

    if not jobs:
        print("No viable jobs found. Check the LR index and parquet gene names.")
        sys.exit(1)

    # -- Run ------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("Starting batch run...")
    print("=" * 60)

    results = run_batch(
        clean_df      = clean,
        jobs          = jobs,
        r_vals        = r_vals,
        window_type   = args.window,
        concave_ratio = args.concave_ratio,
        n_sim         = args.n_sim,
        seed          = args.seed,
        resolution    = args.resolution,
    )

    # -- Save -----------------------------------------------------------------
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    results_df = pd.DataFrame(results)
    results_df.to_parquet(args.out, index=False)

    # Copy r_vals next to the results so downstream analysis has one bundle
    r_vals_out = args.out + ".r_vals.npy"
    np.save(r_vals_out, r_vals)

    ok_count  = (results_df.get("status", pd.Series()) == "ok").sum()
    err_count = len(results_df) - ok_count

    print("\n" + "=" * 60)
    print(f"Batch complete.")
    print(f"  Successful: {ok_count}")
    print(f"  Errors:     {err_count}")
    print(f"  Saved to:   {args.out}")
    print(f"  r_vals at:  {r_vals_out}  (load with np.load(...))")


if __name__ == "__main__":
    main()
