# NoSeggs — Segmentation-Free Spatial Co-Expression Analysis

**MSc Bioinformatics Research Project** | King's College London | 2026

A Python package for segmentation-free spatial co-expression analysis of spatial transcriptomic data. NoSeggs applies bivariate point pattern analysis (Ripley's K/L functions) directly to raw transcript coordinates, bypassing the need for cell segmentation — a major bottleneck in current spatial omics workflows.

---

## Motivation

This project originates from a failed cell segmentation attempt on NanoString CosMx SMI data. Poor DAPI staining quality in the dataset rendered both the platform's native segmentation and third-party resegmentation tools (including [MOSAIK](https://github.com/anthbapt/MOSAIK)) insufficient. Approximately 42% of transcripts in the dataset are unassigned to any cell, and the majority of subcellular compartment labels are missing.

Rather than treating this as a data quality dead-end, NoSeggs treats it as a methodological opportunity: if segmentation is unreliable, don't segment. Instead, analyse the spatial relationships between transcript species directly via point pattern statistics, and use these co-expression metrics to build gene interaction networks.

## Approach

The pipeline works as follows:

1. **Load and QC** raw CosMx transcript tables (~millions of reads with global x/y coordinates, gene labels, FOV assignments)
2. **Assign tissue strips** within each FOV via Gaussian Mixture Models on x-coordinate distributions (control / infected / control layout)
3. **Clean rogue transcripts** that fall outside tissue boundaries using density-based filtering (DBSCAN), so that downstream observation windows reflect actual tissue geometry
4. **Fit observation windows** to cleaned tissue strips — supporting rectangular bounding box, convex hull, and custom GeoJSON polygon geometries — to correct the intensity estimate (λ) and edge correction in the K-function
5. **Compute bivariate Ripley's K and L functions** for gene pairs within each strip, with polygon-aware Shapely-based edge correction
6. **Build permutation envelopes** to distinguish gene-specific co-localisation from tissue-structure confounding
7. **Screen ligand-receptor pairs** (from CellChatDB) across 146 viable pairs × 3 strips
8. *(Planned)* Feed screening results into modularity maximisation to identify co-localisation communities and infection-specific signalling programmes

## The Data

| Property | Value |
|---|---|
| **Platform** | NanoString CosMx Spatial Molecular Imager (SMI) |
| **Tissue** | Severe asthmatic, RSV-infected human lung |
| **Layout** | 4 stitched slides (S1–S4); each slide has 3 tissue strips: control / infected / control |
| **Total transcripts** | ~2.66 million across all slides; 710,751 in S1 (primary analysis slide) |
| **Gene panel** | 1,000 targets + 78 SystemControl probes (~1,078 unique per FOV) |
| **FOVs** | 63 total (4 missing); up to 12 usable in S1 (expanded from original 6 in 02c) |
| **Coordinates** | Global pixel coordinates used exclusively (`x_global_px`, `y_global_px`) |
| **Slide quality** | S1 > S3 > S2/S4 |

## Key Functions

### Core statistical functions (notebooks 03, 09a)

| Function | Description |
|---|---|
| `get_coords(df, gene)` | Extract (N, 2) coordinate array for a given gene |
| `get_window(df, window_type='hull', custom_path=None)` | Return observation window; `'rect'` → tuple, `'hull'` → Shapely Polygon, `'custom'` → loaded polygon |
| `load_custom_window(path)` | Load a user-drawn polygon from JSON (saved by `draw_custom_window` in 09b) |
| `bivariate_k(coords_a, coords_b, r_vals, window, resolution=64)` | Cross-type Ripley's K; dispatches internally on window type (tuple → arc correction, Polygon → Shapely intersection correction) |
| `k_to_l(k_vals, r_vals)` | Variance-stabilising transform: L(r) = √(K(r)/π) − r |
| `compute_envelope(coords_a, coords_b, r_vals, window, n_sim=99)` | Pointwise permutation envelope; window type passed through to `bivariate_k` |
| `run_pair_analysis(strip_df, gene_a, gene_b, r_vals, window_type='hull', custom_path=None)` | Full pipeline wrapper; single entry point regardless of window type |
| `plot_diagnostics(coords_a, coords_b, window, gene_a, gene_b)` | Scatter plot with window overlay |

### Low-level building blocks (notebook 09, re-used by 09a)

| Function | Description |
|---|---|
| `get_convex_hull(df)` | Shapely Polygon convex hull of all transcript coordinates in a strip |
| `fraction_inside_hull(point, r, hull, resolution=64)` | Edge correction weight: fraction of disc of radius r inside hull polygon |
| `draw_custom_window(strip_df, save_path, ...)` | Interactive polygon drawing tool; saves result as JSON (notebook 09b) |

**Window dispatch:** `bivariate_k`, `compute_envelope`, and `run_pair_analysis` all accept a `window` object (tuple or Shapely Polygon). `get_window` is the single point of control for window type selection. The `_hull`-suffixed variants from nb09 are superseded by the unified functions in nb09a.

**Edge correction note:** `fraction_inside_hull` introduces a small systematic negative bias in absolute K(r). This bias is consistent across observed and permuted realisations and cancels in the permutation test. Absolute L(r) values are not interpreted directly; only envelope-relative significance is used.

## Repository Structure

```
cosmx_pointpattern/
│
├── data/
│   ├── raw/                        # Original CosMx files (.zarr) — never modified
│   └── processed/
│       ├── fov3_strips.parquet                    # FOV 3 with GMM strip labels (early work)
│       ├── s1_all_strips.parquet                        # 6 FOVs × 3 strips, ~400k transcripts (02b)
│       ├── s1_all_strips_noise_flagged.parquet          # Cleaned data with is_noise flag (nb08)
│       ├── s1_expanded_strips.parquet                   # All usable FOVs, strip-assigned (02c)
│       ├── s1_expanded_strips_noise_flagged.parquet     # Expanded + DBSCAN QC (02c) — primary input from 09a onwards
│       ├── custom_window_strip_[1-3].json               # User-drawn polygons per strip (nb09b)
│       ├── viable_lr_pairs_all_strips.parquet           # 146 viable CellChatDB L-R pairs with strip counts
│       └── s1_lr_screening_results.parquet              # (planned nb11 output) per-pair co-localisation scores
│
├── notebooks/
│   ├── 01_load_data.ipynb
│   ├── 02_fov3_exploration.ipynb
│   ├── 02c_fov_review_and_expansion.ipynb
│   ├── 03_K_function.ipynb
│   ├── 04_real_analysis.ipynb
│   ├── 05_negative_controls.ipynb
│   ├── 06_LR_checks.ipynb
│   ├── 07_expanded_controls.ipynb
│   ├── 08_improved_QC.ipynb
│   ├── 09_improved_windows_and_edge_correction.ipynb
│   ├── 09a_unified_window_api.ipynb               # unified get_window / bivariate_k / run_pair_analysis
│   ├── 09b_custom_window_drawing.ipynb            # interactive polygon drawing tool
│   ├── 09c_window_comparison.ipynb                # rect vs hull vs custom, visual + L(r) comparison
│   ├── 10_lr_panel_and_network.ipynb              # [planned] L-R screening + network analysis
│   └── 11_network_modularity.ipynb                # [planned] community detection, infection signalling
│
├── src/
│   └── spatialco/                  # Python package (extraction planned after nb12)
│       └── __init__.py
│
├── tests/                          # Unit tests (planned)
├── results/
│   └── figures/                    # Output figures (numbered by notebook)
├── notes/
│   ├── lab_notebook.md             # Chronological working notes and decision log
│   └── intro_methods_reading.md    # Dissertation intro/methods overview and reading list
│
├── .gitignore
├── README.md
└── requirements.txt
```

## Notebook Guide

### Phase 1 — Data Loading & Strip Assignment (01–02b)

**`01_load_data`** loads the full CosMx dataset via SpatialData, characterises the 4-slide structure, identifies S1 as the primary analysis target, and documents the cell assignment failure that motivates the segmentation-free approach.

**`02_fov3_exploration`** narrows to FOV 3 on Slide 1. Visualises x-coordinate distributions to identify the three tissue strips. Initial GMM-based strip separation.

**`02c_fov_review_and_expansion`** revisits the 7 excluded FOVs from 02b. Visual review (scatter + x-histogram) of each excluded FOV with a user-editable `fov_config` dictionary. Re-runs GMM strip assignment for newly included FOVs, merges with existing data and FOV3 (from its separate parquet), then re-runs DBSCAN QC on the full expanded dataset. Outputs: `s1_expanded_strips.parquet` and `s1_expanded_strips_noise_flagged.parquet`.

### Phase 2 — K-Function Development & Validation (03–07)

**`03_K_function`** defines all core statistical functions. Includes identification and fix of a critical edge correction bug (additive arc fractions → correct Ripley's isotropic correction). Validated on synthetic CSR: max |L(r)| < 13 post-fix vs. 195–445 pre-fix.

**`04_real_analysis`** — first positive control: KRT8 × KRT18 on FOV 3.

**`05_negative_controls`** — MALAT1 × KRT18 and KRT8 × SCGB3A1. Establishes that raw L(r) magnitude is dominated by strip geometry; permutation envelope is the only discriminator.

**`06_LR_checks`** cross-references the CosMx panel against CellChatDB (1,912 pairs). 146 pairs viable across all three strips. Output: `viable_lr_pairs_all_strips.parquet`.

**`07_expanded_controls`** re-runs control framework on the full 6-FOV dataset. Identifies interpretability problem with rectangular windows, motivating Phase 3.

### Phase 3 — Pipeline Improvement (08–09)

**`08_improved_QC`** ✅ Removes rogue transcripts via adaptive DBSCAN (1-NN p97, clipped to [20, 30] px). 27,313 flagged as noise (6.8%), 372,925 retained. Output: `s1_all_strips_noise_flagged.parquet`.

**`09_improved_windows_and_edge_correction`** ✅ Replaces rectangular bounding box with convex hull windows (43% of bounding box area for control strips, 27% for infected strip). Implements Shapely-based polygon edge correction. Re-runs three-part control framework on cleaned data — controls behave as expected under permutation envelope. Includes documented limitation on absolute K(r) calibration (bias cancels in permutation test). Note: this notebook defines `_hull`-suffixed variants that are superseded by the unified API in 09a.

**`09a_unified_window_api`** ✅ Consolidates the `_hull`-suffixed functions from nb09 into a unified API. `get_window(df, window_type, custom_path)` returns a tuple (rect) or Shapely Polygon (hull/custom). `bivariate_k`, `compute_envelope`, and `run_pair_analysis` dispatch internally on window type — no `_hull` variants. Uses `s1_expanded_strips_noise_flagged.parquet`. Re-runs controls to confirm parity with nb09.

**`09b_custom_window_drawing`** ✅ Interactive polygon drawing tool (`draw_custom_window`) using `matplotlib.widgets.PolygonSelector`. User draws a polygon on top of the transcript scatter for each strip; result auto-saved as JSON. Per-strip cells for strip_1, strip_2, strip_3. Requires `ipympl` (`pip install ipympl`).

**`09c_window_comparison`** ✅ Three-way visual and L(r) comparison: rectangular bounding box vs convex hull vs user-drawn custom polygon. Overlay figures, area table, and L(r) comparison plots for the positive and negative controls.

### Phase 4 — Screening & Network (10–11) *[Planned]*

**`10_lr_panel_and_network`** Loops over 146 viable L-R pairs using the unified `run_pair_analysis` with `window_type='hull'`, saves per-pair co-localisation scores and envelope significance. Builds co-localisation graph and runs community detection.

**`11_network_modularity`** Detailed network analysis: Louvain/Leiden community structure, comparison of infected vs control strip communities, identification of infection-specific signalling programmes.

## Control Framework

| Control | Gene Pair | Expected Result | Confirmed |
|---|---|---|---|
| **Positive** | KRT8 × KRT18 | L(r) above envelope — co-expressed epithelial keratins | ✅ All 3 strips (rect + hull windows) |
| **Negative 1** | MALAT1 × KRT18 | L(r) inside envelope — ubiquitous housekeeping × epithelial | ✅ All 3 strips (rect + hull windows) |
| **Negative 2** | KRT8 × SCGB3A1 | L(r) inside envelope — different epithelial lineages | ✅ All 3 strips (rect + hull windows) |

The permutation null shuffles gene labels while preserving spatial locations. It tests whether co-localisation is gene-specific or merely a consequence of shared tissue structure.

## Milestones

| Date | Milestone | Status |
|---|---|---|
| Mar 31 | End-to-end pipeline on S1 | ✅ Done |
| Apr 13 | DBSCAN QC + convex hull windows + controls re-run | ✅ Done |
| Apr 14 | FOV expansion (02c) + unified window API (09a) + custom drawing tool (09b) + window comparison (09c) | ✅ Done |
| Apr 30 | L-R panel screening + network analysis (nb10–11) | In progress |
| May 15 | FEATURE FREEZE — background + methods written | |
| May 31 | Packaging, documentation, results finalised | |
| Jun 3 | CODE FROZEN | |
| Jul 16 | SUBMISSION | |

## Requirements

See `requirements.txt`. Core dependencies: `numpy`, `scipy`, `pandas`, `matplotlib`, `scikit-learn`, `shapely`, `spatialdata`, `pyarrow`, `networkx`, `liana`.

## Acknowledgements

Dataset provided by King's College London. Ligand-receptor annotations from [CellChatDB](https://github.com/sqjin/CellChat) (Jin et al., *Nature Communications*, 2021), accessed via [LIANA](https://github.com/saezlab/liana-py) (Dimitrov et al., *Nature Communications*, 2022).
