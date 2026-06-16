# NoSegs — Segmentation-Free Spatial Co-Expression Analysis

**MSc Bioinformatics Research Project** | King's College London | 2026

A research codebase for **segmentation-free ligand–receptor co-localisation analysis** in spatial transcriptomics. NoSegs treats each gene as a 2D point pattern and applies bivariate Ripley's K-function with permutation envelopes directly to raw transcript coordinates — bypassing the cell segmentation step entirely.

> This is a dissertation-track research codebase, not a maintained Python library. The notebooks reproduce the analysis end-to-end.

---

## Why this exists

Cell segmentation is a fragile prerequisite for almost every downstream analysis in imaging-based spatial transcriptomics. On the CosMx sample analysed here, the vendor segmentation drops ~42% of transcripts and misses most subcellular compartment labels, rendering segmentation-dependent tools (CellChat, LIANA, Bento, Squidpy) unusable for ligand–receptor inference.

NoSegs sidesteps the problem: each gene becomes a planar point pattern, and bivariate Ripley's K with label-permutation envelopes asks whether two genes' transcripts spatially co-localise more than expected under a tissue-structure-preserving null. The pipeline scales to panel-level screening (146 ligand–receptor pairs × 3 tissue compartments via a 438-job HPC array).

## Where this sits in the existing toolset

| Tool | Input | Output | Segmentation needed? |
|---|---|---|---|
| Bento (Mah et al. 2024) | transcripts + cell + nuclear masks | subcellular localisation patterns | Yes (mandatory) |
| Squidpy `gr.ripley` | cell coordinates + cluster labels | cluster-level co-localisation | Yes (cells exist) |
| CellChat / LIANA / commot | cell-type expression matrix | inferred L–R signalling | Yes (cells + clusters) |
| Baysor | transcripts | soft cell assignments | Soft segmentation |
| FICTURE (Si et al. 2023) | transcripts | pseudo-tissue factors | No, but doesn't do L–R |
| spatstat (R) | generic 2D point pattern | K, L, etc. | N/A — generic stats |
| **NoSegs** | **raw transcripts only** | **L–R co-localisation + communities** | **No.** |

## Approach

1. **Load and QC** raw CosMx transcript tables (millions of reads with global x/y, gene labels, FOV assignments).
2. **Assign tissue strips** within each FOV via Gaussian Mixture Models on the PCA-rotated x-coordinate (control / infected / control layout).
3. **Clean rogue transcripts** outside tissue boundaries via density-based filtering (DBSCAN) plus a manual cluster-exclusion list, so observation windows reflect actual tissue geometry.
4. **Fit observation windows** to cleaned strips — rectangular, concave hull, or custom polygon — for edge correction in the K-function.
5. **Compute bivariate Ripley's K and L functions** per gene pair, with polygon-aware Shapely-based edge correction.
6. **Build permutation envelopes** to distinguish gene-specific co-localisation from tissue-structure confounding.
7. **Screen ligand–receptor pairs** (CellChatDB intersection, 146 viable pairs) across 3 strips, via SLURM array on KCL CREATE.
8. **Network and community detection** on the screened pairs (modularity-based) to identify co-localised signalling programmes.

## The Data

| Property | Value |
|---|---|
| Platform | NanoString CosMx Spatial Molecular Imager (SMI) |
| Tissue | Severe asthmatic, RSV-infected human lung |
| Layout | 4 stitched slides (S1–S4); each slide has 3 tissue strips: control / infected / control |
| Total transcripts | ~2.66 million across slides; 703,000 in S1 (primary analysis slide) after QC |
| Gene panel | 1,000 targets + 78 SystemControl probes |
| FOVs analysed | 6 FOVs (4, 5, 8, 9, 10, 11) on S1; 3 strips per FOV |
| Coordinates | Global pixel coordinates (`x_global_px`, `y_global_px`) |

## Repository layout

```
cosmx-pointpattern/
├── data/
│   ├── raw/                                       # CosMx zarr (never modified)
│   └── processed/                                 # Cleaned parquets, panel index
├── notebooks/
│   ├── 00_functions.ipynb                         # Shared function library
│   ├── 01_load_data.ipynb                         # Raw zarr ingest
│   ├── 02_fov3_exploration.ipynb                  # GMM strip assignment (FOV 3)
│   ├── 02c_fov_review_and_expansion.ipynb         # 6-FOV expansion
│   ├── 03_K_function.ipynb                        # Bivariate K + L + envelope
│   ├── 04_real_analysis.ipynb                     # KRT8×KRT18 positive control
│   ├── 05_negative_control_validation.ipynb       # MALAT1×KRT18, KRT8×SCGB3A1
│   ├── 06_LR_checks.ipynb                         # 146 viable CellChatDB pairs
│   ├── 07_expanded_controls.ipynb                 # Controls on full 6-FOV dataset
│   ├── 08_improved_QC.ipynb                       # DBSCAN noise-removal method development
│   ├── 08b_manual_qc_cleanup.ipynb                # QC dev (superseded by nb10 for the final cleaning)
│   ├── 09_improved_windows_and_edge_correction.ipynb   # Concave hull + Shapely
│   ├── 09a_unified_window_api.ipynb               # get_window dispatcher
│   ├── 09c_window_comparison.ipynb                # rect vs hull vs custom
│   ├── 10_poc_end_to_end.ipynb                    # POC + CANONICAL cleaner (writes s1_all_strips_cleaned.parquet; 18-cluster manual exclusion)
│   ├── 11_results_overview.ipynb                  # Aggregated panel results
│   ├── 12_network_analysis.ipynb                  # Network + community detection
│   ├── 13_diagnostics.ipynb                       # Sparse-signal diagnostics (WIP)
│   └── archive/                                   # Superseded notebooks (e.g. 09b)
├── scripts/
│   ├── batch_k_analysis.py                        # HPC single-job runner
│   ├── run_array.sh                               # 438-job SLURM array
│   ├── dry_run.sh                                 # Single-pair smoke test
│   ├── aggregate.py                               # per_pair shards → panel parquet
│   └── fig_vendor_dropout.py                      # 42% segmentation-failure figure
├── results/
│   ├── lr_panel_results.parquet                   # Aggregated panel screen
│   ├── hpc_job_manifest.csv                       # 438 (pair, strip) job rows
│   └── figures/                                   # Output figures (numbered by notebook)
├── notes/
│   ├── lab_notebook.md                            # Chronological working log
│   ├── notebook_audit.md                          # Per-notebook status (current)
│   ├── network_revival_plan.md                    # Diagnostic roadmap for sparse signal
│   ├── extended_panel_rationale.md                # Pre-registered anchor pairs
│   ├── intro_methods_reading.md                   # Dissertation reading list
│   ├── future_work.md                             # Post-deadline backlog
│   └── STATUS.md                                  # Current status board
├── planning_docs/                                 # Living planning + decisions log
├── docs/
│   ├── NoSegs_Poster.pptx                        # NeuroMonster A0 poster (active)
│   └── archive/                                   # Old talks
├── README.md
├── requirements.txt
└── .gitignore
```

## Control framework

| Control | Gene Pair | Expectation | Confirmed |
|---|---|---|---|
| Positive | KRT8 × KRT18 | Strongest co-localisation (co-expressed epithelial keratins) | ✅ All 3 strips |
| Negative 1 | MALAT1 × KRT18 | No gene-specific co-localisation (ubiquitous × epithelial) | ✅ All 3 strips |
| Negative 2 | KRT8 × SCGB3A1 | No co-localisation; distinct compartments (epithelial vs secretory airway) | ✅ All 3 strips (concave hull, n_sim=199; re-run 2026-06-16) |

The permutation null shuffles gene labels while preserving spatial locations — it tests whether co-localisation is gene-specific or merely a consequence of shared tissue structure.

> **Window/SES note.** The early control validation (rectangular window, `n_sim=99`) showed KRT8×KRT18 above and the negatives within the envelope. Under the **production concave-hull window** the binary envelope test is under-powered for *all* pairs (the tight window concentrates the label-swap null — see `notes/decision_log.md`, Step A.6), so controls are discriminated by **standardised effect size (SES)**: KRT8×KRT18 sits +2.67 envelope half-widths above MALAT1×KRT18, and KRT8×SCGB3A1 sits furthest below (signed-peak SES −6.92 / −4.81 / −3.02 across strips 1/2/3), i.e. spatial anti-association. Reproduce with `scripts/run_controls.sh`. Full record: `Final_Writeup/METHODS_LOG.md`.

## Key functions (in `notebooks/00_functions.ipynb`)

| Function | Description |
|---|---|
| `get_coords(df, gene)` | Extract (N, 2) coordinate array for a given gene |
| `get_window(df, window_type='hull', custom_path=None)` | Return observation window: tuple (rect) or Shapely Polygon (hull/custom) |
| `bivariate_k(coords_a, coords_b, r_vals, window, resolution=64)` | Cross-type Ripley's K; dispatches internally on window type |
| `k_to_l(k_vals, r_vals)` | Variance-stabilising transform: L(r) = √(K(r)/π) − r |
| `compute_envelope(coords_a, coords_b, r_vals, window, n_sim=99)` | Pointwise permutation envelope |
| `run_pair_analysis(strip_df, gene_a, gene_b, r_vals, window_type='hull')` | Full pipeline wrapper; single entry point |
| `fraction_inside_hull(point, r, hull, resolution=64)` | Edge correction weight for a disc of radius r at a query point |

**Edge correction note.** `fraction_inside_hull` introduces a small systematic negative bias in absolute K(r). The bias is consistent across observed and permuted realisations and cancels in the permutation test. Absolute L(r) values are not interpreted directly; only envelope-relative significance is used.

## Status

The pipeline is technically complete end-to-end and validated against positive and negative controls. The 438-job HPC array (146 pairs × 3 strips) at `n_sim=199` completed with no failed jobs.

**Open scientific question.** Only 1 of 438 pairs (SPP1 × CD44, strip 2) passes the `n_consec=3` significance filter at production envelope width. Diagnostic work in `nb13` (Step A.1–A.7 of [`notes/network_revival_plan.md`](notes/network_revival_plan.md)) is in flight to distinguish two hypotheses: (a) the production `n_sim=199` envelope is over-strict compared to the `n_sim=99` validation regime; (b) the CellChatDB panel excludes the strongest co-localising pairs by construction. See [`notes/STATUS.md`](notes/STATUS.md) for the live status.

## Reproducing the analysis

1. Create the environment: `python -m venv .venv && pip install -r requirements.txt`.
2. Run notebooks `01 → 13` in order. Notebooks `00` (functions) is loaded by other notebooks via `nbformat`.
3. For the HPC L-R panel screen: edit `scripts/run_array.sh` to point at your venv and submit via `sbatch`. Aggregate with `python scripts/aggregate.py`.

## Acknowledgements

Dataset provided by King's College London. Ligand-receptor annotations from CellChatDB (Jin et al., *Nature Communications*, 2021). HPC compute on KCL CREATE (`msc_appbio` partition).

## Citation

This work is a 2026 MSc Bioinformatics dissertation. Citation block will be added on submission.
