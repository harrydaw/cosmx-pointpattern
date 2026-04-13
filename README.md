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
| **FOVs** | 63 total (4 missing); 6 retained in S1 after QC (~400,000 transcripts) |
| **Coordinates** | Global pixel coordinates used exclusively (`x_global_px`, `y_global_px`) |
| **Slide quality** | S1 > S3 > S2/S4 |

## Key Functions

### Rectangular window (notebook 03)

| Function | Description |
|---|---|
| `get_coords(df, gene)` | Extract (N, 2) coordinate array for a given gene |
| `get_window(df, window_type='rect')` | Compute observation window; dispatches to rect, hull, or custom |
| `bivariate_k(coords_a, coords_b, r_vals, window)` | Cross-type Ripley's K with rectangular isotropic edge correction |
| `k_to_l(k_vals, r_vals)` | Variance-stabilising transform: L(r) = √(K(r)/π) − r |
| `compute_envelope(coords_a, coords_b, r_vals, window, n_sim=99)` | Pointwise permutation envelope |
| `run_pair_analysis(strip_df, gene_a, gene_b, r_vals)` | Full pipeline wrapper for one gene pair on one strip |
| `plot_diagnostics(coords_a, coords_b, window, gene_a, gene_b)` | Scatter plot with window overlay |

### Convex hull window (notebook 09)

| Function | Description |
|---|---|
| `get_convex_hull(df)` | Shapely Polygon convex hull of all transcript coordinates in a strip |
| `fraction_inside_hull(point, r, hull, resolution=64)` | Edge correction weight: fraction of disc at radius r inside hull polygon |
| `bivariate_k_hull(coords_a, coords_b, r_vals, hull)` | K-function with Shapely-based polygon edge correction; normalises by hull.area |
| `compute_envelope_hull(coords_a, coords_b, r_vals, hull, n_sim=99)` | Permutation envelope using hull-aware K estimator |
| `run_pair_analysis_hull(strip_df, gene_a, gene_b, r_vals)` | Full pipeline wrapper using convex hull window |

**Edge correction note:** `fraction_inside_hull` introduces a small systematic negative bias in absolute K(r). This bias is consistent across observed and permuted realisations and cancels in the permutation test. Absolute L(r) values are not interpreted directly; only envelope-relative significance is used.

## Repository Structure

```
cosmx_pointpattern/
│
├── data/
│   ├── raw/                        # Original CosMx files (.zarr) — never modified
│   └── processed/
│       ├── fov3_strips.parquet                    # FOV 3 with GMM strip labels (early work)
│       ├── s1_all_strips.parquet                  # 6 FOVs × 3 strips, ~400k transcripts
│       ├── s1_all_strips_noise_flagged.parquet    # Cleaned data with is_noise flag (nb08 output)
│       ├── viable_lr_pairs_all_strips.parquet     # 146 viable CellChatDB L-R pairs with strip counts
│       └── s1_lr_screening_results.parquet        # (planned nb11 output) per-pair co-localisation scores
│
├── notebooks/
│   ├── 01_load_data.ipynb
│   ├── 02_fov3_exploration.ipynb
│   ├── 02b_strip_assignment_all_fovs.ipynb
│   ├── 03_K_function.ipynb
│   ├── 04_real_analysis.ipynb
│   ├── 05_negative_controls.ipynb
│   ├── 06_LR_checks.ipynb
│   ├── 07_expanded_controls.ipynb
│   ├── 08_improved_QC.ipynb
│   ├── 09_improved_windows_and_edge_correction.ipynb
│   ├── 10_window_comparison_all_fovs.ipynb        # [planned] rect vs hull vs custom, all FOVs
│   ├── 11_lr_screening.ipynb                      # [planned] loop over 146 L-R pairs
│   └── 12_network_modularity.ipynb                # [planned] community detection, infection signalling
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

**`02b_strip_assignment_all_fovs`** extends strip assignment to all 13 S1 FOVs via a two-pass QC pipeline. 6 FOVs retained (~400k transcripts), 7 excluded. Output: `s1_all_strips.parquet`.

### Phase 2 — K-Function Development & Validation (03–07)

**`03_K_function`** defines all core statistical functions. Includes identification and fix of a critical edge correction bug (additive arc fractions → correct Ripley's isotropic correction). Validated on synthetic CSR: max |L(r)| < 13 post-fix vs. 195–445 pre-fix.

**`04_real_analysis`** — first positive control: KRT8 × KRT18 on FOV 3.

**`05_negative_controls`** — MALAT1 × KRT18 and KRT8 × SCGB3A1. Establishes that raw L(r) magnitude is dominated by strip geometry; permutation envelope is the only discriminator.

**`06_LR_checks`** cross-references the CosMx panel against CellChatDB (1,912 pairs). 146 pairs viable across all three strips. Output: `viable_lr_pairs_all_strips.parquet`.

**`07_expanded_controls`** re-runs control framework on the full 6-FOV dataset. Identifies interpretability problem with rectangular windows, motivating Phase 3.

### Phase 3 — Pipeline Improvement (08–09)

**`08_improved_QC`** ✅ Removes rogue transcripts via adaptive DBSCAN (1-NN p97, clipped to [20, 30] px). 27,313 flagged as noise (6.8%), 372,925 retained. Output: `s1_all_strips_noise_flagged.parquet`.

**`09_improved_windows_and_edge_correction`** ✅ Replaces rectangular bounding box with convex hull windows (43% of bounding box area for control strips, 27% for infected strip). Implements Shapely-based polygon edge correction. Re-runs three-part control framework on cleaned data — controls behave as expected under permutation envelope. Includes documented limitation on absolute K(r) calibration (bias cancels in permutation test).

### Phase 4 — Screening & Network (10–12) *[Planned]*

**`10_window_comparison_all_fovs`** Three-way window comparison (rect, hull, custom GeoJSON) across all 6 FOVs × 3 strips. Adds `get_window(window_type=...)` unified entry point and `load_custom_window()` GeoJSON importer.

**`11_lr_screening`** Loops over 146 viable L-R pairs, runs `run_pair_analysis_hull` per strip, saves co-localisation scores and envelope significance. Output: `s1_lr_screening_results.parquet`.

**`12_network_modularity`** Builds co-localisation graph (nodes = genes, edges = significant pairs), runs Louvain/Leiden community detection, compares infected vs control strip community structure to identify infection-specific signalling programmes.

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
| Apr 30 | Window finalisation + L-R screening + network development | In progress |
| May 15 | FEATURE FREEZE — background + methods written | |
| May 31 | Packaging, documentation, results finalised | |
| Jun 3 | CODE FROZEN | |
| Jul 16 | SUBMISSION | |

## Requirements

See `requirements.txt`. Core dependencies: `numpy`, `scipy`, `pandas`, `matplotlib`, `scikit-learn`, `shapely`, `spatialdata`, `pyarrow`, `networkx`, `liana`.

## Acknowledgements

Dataset provided by King's College London. Ligand-receptor annotations from [CellChatDB](https://github.com/sqjin/CellChat) (Jin et al., *Nature Communications*, 2021), accessed via [LIANA](https://github.com/saezlab/liana-py) (Dimitrov et al., *Nature Communications*, 2022).
