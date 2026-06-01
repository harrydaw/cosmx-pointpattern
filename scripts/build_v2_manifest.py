#!/usr/bin/env python
"""
Build the v2 (extended panel) job manifest for the broadened HPC re-run.

Hard-codes the 81 pre-registered pairs from
``notes/extended_panel_rationale.md`` (kept here for reproducibility — the
markdown is the rationale, this file is the executable source of truth).

Pipeline:
  1. Load ``data/processed/s1_all_strips_cleaned.parquet`` and apply the
     same boolean filters as ``batch_k_analysis.py`` (drop noise / small
     cluster / manually excluded).
  2. For each (ligand, receptor) pair, drop the pair if either gene is
     absent from the panel.
  3. For each (pair, strip), drop the (pair, strip) row if either gene
     has fewer than ``MIN_TRANSCRIPTS_PER_GENE_PER_STRIP`` (default 50)
     transcripts in that strip.
  4. Write surviving rows to ``results/hpc_job_manifest_v2.csv`` with
     columns ``ligand,receptor,strip`` matching the v1 manifest schema.
  5. Print a verbose drop log + the exact sbatch command to use.

Usage:
    python scripts/build_v2_manifest.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

MIN_TRANSCRIPTS_PER_GENE_PER_STRIP = 50
STRIPS = ["strip_1", "strip_2", "strip_3"]
DATA = Path("data/processed/s1_all_strips_cleaned.parquet")
OUT = Path("results/hpc_job_manifest_v2.csv")

# ---------------------------------------------------------------------------
# Pre-registered pair list (81 pairs)
# Source: notes/extended_panel_rationale.md
# ---------------------------------------------------------------------------

EXTENDED_PAIRS: list[tuple[str, str, str]] = [
    # (ligand, receptor, group)
    # Group 1 — Cytokeratin co-expression anchors (n=10)
    ("KRT8",  "KRT18",  "G1_cytokeratin"),
    ("KRT5",  "KRT14",  "G1_cytokeratin"),
    ("KRT5",  "KRT15",  "G1_cytokeratin"),
    ("KRT8",  "KRT19",  "G1_cytokeratin"),
    ("KRT18", "KRT19",  "G1_cytokeratin"),
    ("KRT14", "KRT15",  "G1_cytokeratin"),
    ("KRT5",  "KRT17",  "G1_cytokeratin"),
    ("KRT14", "KRT17",  "G1_cytokeratin"),
    ("KRT19", "EPCAM",  "G1_cytokeratin"),
    ("KRT8",  "EPCAM",  "G1_cytokeratin"),
    # Group 2 — Epithelial cell-type co-markers (n=10)
    ("SCGB1A1", "SCGB3A1", "G2_epi_celltype"),
    ("SCGB1A1", "SCGB3A2", "G2_epi_celltype"),
    ("MUC5AC",  "MUC5B",   "G2_epi_celltype"),
    ("MUC5AC",  "TFF3",    "G2_epi_celltype"),
    ("FOXJ1",   "TUBB4B",  "G2_epi_celltype"),
    ("FOXJ1",   "DNAH5",   "G2_epi_celltype"),  # verify DNAH5 in panel
    ("SFTPC",   "SFTPB",   "G2_epi_celltype"),
    ("SFTPC",   "SFTPA1",  "G2_epi_celltype"),
    ("AGER",    "HOPX",    "G2_epi_celltype"),
    ("KRT5",    "TP63",    "G2_epi_celltype"),
    # Group 3 — Epithelial niche overlap (n=8)
    ("KRT8",  "SCGB1A1", "G3_epi_niche"),
    ("KRT18", "SCGB1A1", "G3_epi_niche"),
    ("KRT8",  "MUC5B",   "G3_epi_niche"),
    ("KRT18", "MUC5AC",  "G3_epi_niche"),
    ("KRT19", "SCGB3A1", "G3_epi_niche"),
    ("EPCAM", "MUC5B",   "G3_epi_niche"),
    ("EPCAM", "SCGB3A1", "G3_epi_niche"),
    ("EPCAM", "FOXJ1",   "G3_epi_niche"),
    # Group 4 — Immune cell-type co-markers (n=10)
    ("CD3D",  "CD3E",    "G4_immune_celltype"),
    ("CD3D",  "CD4",     "G4_immune_celltype"),
    ("CD3D",  "CD8A",    "G4_immune_celltype"),
    ("CD8A",  "CD8B",    "G4_immune_celltype"),
    ("CD4",   "FOXP3",   "G4_immune_celltype"),  # verify FOXP3
    ("CD68",  "CD163",   "G4_immune_celltype"),
    ("CD68",  "CD14",    "G4_immune_celltype"),
    ("CD19",  "MS4A1",   "G4_immune_celltype"),
    ("ITGAX", "HLA-DRA", "G4_immune_celltype"),
    ("NCAM1", "NKG7",    "G4_immune_celltype"),  # verify NCAM1, NKG7
    # Group 5 — Stromal co-expression (n=8)
    ("COL1A1", "COL3A1", "G5_stromal"),
    ("COL1A1", "COL1A2", "G5_stromal"),
    ("ACTA2",  "TAGLN",  "G5_stromal"),
    ("ACTA2",  "MYH11",  "G5_stromal"),
    ("PECAM1", "VWF",    "G5_stromal"),
    ("PECAM1", "CDH5",   "G5_stromal"),
    ("VWF",    "CDH5",   "G5_stromal"),
    ("DCN",    "LUM",    "G5_stromal"),
    # Group 6 — Negative anchors (n=10)
    ("MALAT1", "KRT18",  "G6_negative_anchor"),
    ("MALAT1", "CD3D",   "G6_negative_anchor"),
    ("MALAT1", "COL1A1", "G6_negative_anchor"),
    ("KRT8",   "CD3D",   "G6_negative_anchor"),
    ("KRT8",   "COL1A1", "G6_negative_anchor"),
    ("KRT8",   "PECAM1", "G6_negative_anchor"),
    ("CD3D",   "COL1A1", "G6_negative_anchor"),
    ("MUC5AC", "CD68",   "G6_negative_anchor"),
    ("SFTPC",  "ACTA2",  "G6_negative_anchor"),
    ("SFTPC",  "COL1A1", "G6_negative_anchor"),
    # Group 7 — Infection-response candidates (n=15)
    ("ISG15",   "IFIT1",  "G7_infection"),
    ("IFIT1",   "IFIT3",  "G7_infection"),
    ("ISG15",   "MX1",    "G7_infection"),
    ("OAS1",    "ISG15",  "G7_infection"),
    ("ISG15",   "KRT8",   "G7_infection"),
    ("IFIT1",   "KRT18",  "G7_infection"),
    ("MX1",     "KRT5",   "G7_infection"),
    ("CXCL10",  "CXCL11", "G7_infection"),
    ("CXCL8",   "IL6",    "G7_infection"),
    ("IL6",     "IL1B",   "G7_infection"),
    ("TNF",     "IL1B",   "G7_infection"),
    ("ISG15",   "CD68",   "G7_infection"),
    ("CXCL10",  "CD3D",   "G7_infection"),
    ("HLA-DRA", "CD3D",   "G7_infection"),
    ("HLA-DRA", "CD68",   "G7_infection"),
    # Group 8 — Exploratory / paracrine signalling (n=10)
    ("TGFB1",  "TGFBR2",   "G8_paracrine"),
    ("CXCL12", "CXCR4",    "G8_paracrine"),
    ("CCL2",   "CCR2",     "G8_paracrine"),
    ("CCL5",   "CCR5",     "G8_paracrine"),
    ("IL1B",   "IL1R1",    "G8_paracrine"),
    ("TNF",    "TNFRSF1A", "G8_paracrine"),
    ("WNT5A",  "FZD2",     "G8_paracrine"),
    ("FGF7",   "FGFR2",    "G8_paracrine"),
    ("HGF",    "MET",      "G8_paracrine"),
    ("EGF",    "EGFR",     "G8_paracrine"),
]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    if not DATA.exists():
        sys.exit(f"ERROR: data file not found: {DATA}")

    print(f"Loading: {DATA}")
    df = pd.read_parquet(DATA)
    clean = df[
        ~df["is_noise"]
        & ~df["is_small_cluster"]
        & ~df["manually_excluded"]
    ].copy()
    print(f"  clean transcripts: {len(clean):,}")

    available_genes = set(clean["target"].unique())
    print(f"  unique genes:      {len(available_genes):,}")

    # Per-gene per-strip transcript counts
    counts = (
        clean.groupby(["strip", "target"])
        .size()
        .unstack(fill_value=0)
    )

    print(f"\nPre-registered pairs: {len(EXTENDED_PAIRS)}")
    print(f"Strips: {STRIPS}")
    print(f"Filter: each gene must have >= {MIN_TRANSCRIPTS_PER_GENE_PER_STRIP} transcripts in the given strip.\n")

    surviving: list[tuple[str, str, str]] = []
    dropped_pairs: list[tuple[str, str, str]] = []  # (lig, rec, reason)
    dropped_rows: list[tuple[str, str, str, str]] = []  # (lig, rec, strip, reason)

    for lig, rec, group in EXTENDED_PAIRS:
        # Gene-availability check (drops the whole pair across all strips)
        missing = [g for g in (lig, rec) if g not in available_genes]
        if missing:
            dropped_pairs.append((lig, rec, f"gene(s) absent from panel: {missing}"))
            continue

        # Per-strip abundance check
        for strip in STRIPS:
            n_lig = int(counts.loc[strip].get(lig, 0))
            n_rec = int(counts.loc[strip].get(rec, 0))
            if n_lig < MIN_TRANSCRIPTS_PER_GENE_PER_STRIP:
                dropped_rows.append((lig, rec, strip,
                                     f"{lig} n={n_lig} < {MIN_TRANSCRIPTS_PER_GENE_PER_STRIP}"))
                continue
            if n_rec < MIN_TRANSCRIPTS_PER_GENE_PER_STRIP:
                dropped_rows.append((lig, rec, strip,
                                     f"{rec} n={n_rec} < {MIN_TRANSCRIPTS_PER_GENE_PER_STRIP}"))
                continue
            surviving.append((lig, rec, strip))

    # Write manifest
    OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(surviving, columns=["ligand", "receptor", "strip"]).to_csv(
        OUT, index=False
    )

    print(f"--- Dropped pairs (whole pair, gene-availability) ---")
    if dropped_pairs:
        for d in dropped_pairs:
            print(f"  {d[0]:>8} x {d[1]:<8}  {d[2]}")
    else:
        print("  (none)")

    print(f"\n--- Dropped (pair, strip) rows (abundance < {MIN_TRANSCRIPTS_PER_GENE_PER_STRIP}) ---")
    if dropped_rows:
        for d in dropped_rows:
            print(f"  {d[0]:>8} x {d[1]:<8}  {d[2]}  {d[3]}")
    else:
        print("  (none)")

    n_jobs = len(surviving)
    print(f"\n--- Summary ---")
    print(f"  Pre-registered pairs:     {len(EXTENDED_PAIRS)}")
    print(f"  Pairs dropped (genes):    {len(dropped_pairs)}")
    print(f"  Rows dropped (abundance): {len(dropped_rows)}")
    print(f"  Surviving manifest rows:  {n_jobs}")
    print(f"  Wrote:                    {OUT}")

    if n_jobs == 0:
        sys.exit("ERROR: no surviving jobs.")

    print(f"\nOn HPC, submit with:")
    print(f"  sbatch --array=1-{n_jobs}%2 scripts/run_array_v2.sh")


if __name__ == "__main__":
    main()
