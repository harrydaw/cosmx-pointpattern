#!/usr/bin/env python
"""
Aggregate v2 per-job parquets into a single results file.

Mirrors scripts/aggregate.py but points at the v2 manifest, v2 shard
directory, and writes to results/lr_panel_results_v2.parquet. The original
results/lr_panel_results.parquet (the CellChatDB 438-job run at n_sim=199)
is never touched.

Usage:
    python scripts/aggregate_v2.py
    python scripts/aggregate_v2.py --keep-shards
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

MANIFEST = "results/hpc_job_manifest_v2.csv"
SHARD_DIR = "results/per_pair_v2"
OUT_PARQUET = "results/lr_panel_results_v2.parquet"
OUT_FAILED = "results/failed_jobs_v2.csv"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--keep-shards", action="store_true",
                        help="Do not print the shard-deletion command at the end.")
    args = parser.parse_args()

    if not os.path.exists(MANIFEST):
        sys.exit(f"Manifest missing: {MANIFEST}")

    manifest = pd.read_csv(MANIFEST)
    print(f"v2 manifest rows: {len(manifest)}")

    rows, failed, r_vals_seen = [], [], None

    for _, m in manifest.iterrows():
        lig, rec, strip = m["ligand"], m["receptor"], m["strip"]
        shard = f"{SHARD_DIR}/{lig}_{rec}_{strip}.parquet"

        if not os.path.exists(shard):
            failed.append({"ligand": lig, "receptor": rec, "strip": strip,
                           "reason": "shard_missing"})
            continue

        try:
            df = pd.read_parquet(shard)
        except Exception as e:
            failed.append({"ligand": lig, "receptor": rec, "strip": strip,
                           "reason": f"read_error: {e}"})
            continue

        if len(df) != 1:
            failed.append({"ligand": lig, "receptor": rec, "strip": strip,
                           "reason": f"expected_1_row_got_{len(df)}"})
            continue

        row = df.iloc[0].to_dict()
        if row.get("status") != "ok":
            failed.append({"ligand": lig, "receptor": rec, "strip": strip,
                           "reason": f"status:{row.get('status')}"})
            continue

        rows.append(row)

        rvf = shard + ".r_vals.npy"
        if os.path.exists(rvf):
            rv = np.load(rvf)
            if r_vals_seen is None:
                r_vals_seen = rv
            elif not np.array_equal(rv, r_vals_seen):
                failed.append({"ligand": lig, "receptor": rec, "strip": strip,
                               "reason": "r_vals_mismatch"})

    print(f"Successful: {len(rows)}")
    print(f"Failed:     {len(failed)}")

    if rows:
        out_df = pd.DataFrame(rows)
        out_df.to_parquet(OUT_PARQUET, index=False)
        print(f"Wrote: {OUT_PARQUET}")

        if r_vals_seen is not None:
            np.save(OUT_PARQUET + ".r_vals.npy", r_vals_seen)
            print(f"Wrote: {OUT_PARQUET}.r_vals.npy")

    pd.DataFrame(failed, columns=["ligand", "receptor", "strip", "reason"]).to_csv(
        OUT_FAILED, index=False
    )
    print(f"Wrote: {OUT_FAILED}")

    if failed:
        print("\nRe-run failing tasks by computing manifest row numbers:")
        print("  python -c \"import pandas as pd; "
              "m = pd.read_csv('results/hpc_job_manifest_v2.csv'); "
              "f = pd.read_csv('results/failed_jobs_v2.csv'); "
              "ids = [m.index[(m.ligand==r.ligand)&(m.receptor==r.receptor)&(m.strip==r.strip)][0]+1 "
              "for _, r in f.iterrows()]; "
              "print(','.join(map(str, ids)))\"")
        print("Then: sbatch --array=<comma-separated-ids> scripts/run_array_v2.sh")
    elif not args.keep_shards:
        print("\nAll jobs ok. Once verified, prune shards to free inodes:")
        print(f"  rm -r {SHARD_DIR}")


if __name__ == "__main__":
    main()
