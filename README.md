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
4. **Fit observation windows** to cleaned tissue strips — supporting rectangular, convex hull, and alpha shape (concave hull) geometries — to correct the intensity estimate (λ) and edge correction in the K-function
5. **Compute bivariate Ripley's K and L functions** for gene pairs within each strip, with Ripley's isotropic edge correction
6. **Build permutation envelopes** to distinguish gene-specific co-localisation from tissue-structure confounding
7. **Screen gene panels** (including ligand-receptor pairs from CellChatDB) across multiple FOVs
8. *(Planned)* Feed the resulting co-expression matrix into modularity maximisation to identify gene modules

## The Data

| Property | Value |
|---|---|
| **Platform** | NanoString CosMx Spatial Molecular Imager (SMI) |
| **Tissue** | Severe asthmatic, RSV-infected human lung |
| **Layout** | 4 stitched slides (S1–S4); each slide has 3 tissue strips: control / infected / control |
| **Total transcripts** | ~2.66 million across all slides; 710,751 in S1 (primary analysis slide) |
| **Gene panel** | 1,000 targets + 78 SystemControl probes (~1,078 unique per FOV) |
| **FOVs** | 63 total (4 missing); 6 retained in S1 after QC (~400,000 transcripts) |
| **Coordinates** | Global pixel coordinates used exclusively (local FOV coordinates are meaningless across FOVs) |
| **Slide quality** | S1 > S3 > S2/S4 |

## Key Functions

| Function | Description |
|---|---|
| `get_coords(df, gene)` | Extract (N, 2) coordinate array for a given gene from a transcript DataFrame |
| `get_window(df)` | Compute the observation window (rectangular bounding box) for a transcript DataFrame |
| `bivariate_k(coords_a, coords_b, r_vals, window)` | Cross-type Ripley's K function with isotropic edge correction via `_fraction_inside_rect()` |
| `k_to_l(k_vals, r_vals)` | Variance-stabilising transform: L(r) = √(K(r)/π) − r |
| `compute_envelope(coords_a, coords_b, r_vals, window, n_sim=99)` | Pointwise permutation envelope (label shuffling, not Poisson simulation) |
| `run_pair_analysis(strip_df, gene_a, gene_b, r_vals, n_sim=99)` | Full pipeline wrapper for one gene pair on one strip |
| `plot_diagnostics(coords_a, coords_b, window, gene_a, gene_b)` | Scatter plot with bounding-box overlay |

## Repository Structure

```
cosmx_pointpattern/
│
├── data/
│   ├── raw/                        # Original CosMx files (.zarr) — never modified
│   └── processed/                  # Pipeline outputs
│       ├── fov3_strips.parquet         # FOV 3 with GMM strip labels (early work)
│       ├── s1_all_strips.parquet       # 6 FOVs × 3 strips, ~400k transcripts
│       └── s1_all_strips_cleaned.parquet  # (08_ output) cleaned, with noise flag column
│
├── notebooks/
│   ├── 01_load_data.ipynb              # Initial data exploration and discovery
│   ├── 02_fov3_exploration.ipynb       # Focused EDA on S1, FOV 3 (most abundant)
│   ├── 02b_strip_assignment_all_fovs.ipynb  # GMM strip assignment across all S1 FOVs
│   ├── 03_K_function.ipynb             # Core K/L function definitions and edge correction bug fix
│   ├── 04_real_analysis.ipynb          # First positive control (KRT8 × KRT18)
│   ├── 05_negative_controls.ipynb      # Negative controls (MALAT1 × KRT18, KRT8 × SCGB3A1)
│   ├── 06_LR_checks.ipynb             # CellChatDB L-R panel cross-reference, multi-FOV expansion
│   ├── 07_expanded_controls.ipynb      # Extended controls + synthetic CSR validation
│   ├── 08_improved_QC.ipynb            # [WIP] Density-based rogue transcript removal (DBSCAN)
│   └── 09_window_functions.ipynb       # [WIP] Polygonal/alpha-shape windows + updated edge correction
│
├── src/
│   └── spatialco/                  # Future Python package
│       └── __init__.py
│
├── tests/                          # Unit tests (planned)
├── results/                        # Output figures and summary CSVs
├── notes/
│   └── lab_notebook.md             # Chronological working notes and decision log
│
├── .gitignore
├── README.md
└── requirements.txt
```

## Notebook Guide

### Phase 1 — Data Loading & Strip Assignment (01–02b)

**`01_load_data`** loads the full CosMx dataset via SpatialData, characterises the 4-slide structure, identifies S1 as the primary analysis target, and documents the cell assignment failure that motivates the segmentation-free approach.

**`02_fov3_exploration`** narrows to FOV 3 on Slide 1 (the most transcript-dense FOV). Visualises x-coordinate distributions to identify the three tissue strips. Initial GMM-based strip separation.

**`02b_strip_assignment_all_fovs`** extends strip assignment to all 13 S1 FOVs via a two-pass QC pipeline: visual inspection of each FOV, then GMM assignment with per-FOV decisions recorded in a `fov_config` dictionary. 6 FOVs retained (~400k transcripts), 7 excluded. Output: `s1_all_strips.parquet`.

### Phase 2 — K-Function Development & Validation (03–05)

**`03_K_function`** defines all core statistical functions. Includes the identification and fix of a critical edge correction bug (additive arc fractions → correct Ripley's isotropic correction). Validated on synthetic CSR data: max |L(r)| < 13 post-fix vs. 195–445 pre-fix.

**`04_real_analysis`** runs the first positive control: KRT8 × KRT18 (co-expressed epithelial keratins). L(r) exceeds permutation envelope in strip 2, confirming gene-specific co-localisation.

**`05_negative_controls`** runs two negative controls. MALAT1 × KRT18: high raw L(r) but inside envelope (ubiquitous housekeeping gene, tissue geometry only). KRT8 × SCGB3A1: same result (different epithelial lineages, but spatially interleaved). Key insight: in narrow strips, raw L(r) magnitude is dominated by geometry; only the permutation envelope discriminates gene-specific signal.

### Phase 3 — Panel Expansion & Feasibility (06–07)

**`06_LR_checks`** cross-references the CosMx gene panel against CellChatDB (1,912 curated L-R pairs). Single-FOV analysis is infeasible at the n≥50 threshold. Multi-FOV expansion to 6 FOVs yields 199/203 L-R pairs viable in at least one strip, 146 across all three strips.

**`07_expanded_controls`** re-runs the control framework on the expanded 6-FOV panel and includes synthetic data validation. Identifies the core interpretability problem: bounding-box windows violate stationarity/isotropy assumptions, making results difficult to explain. Motivates the rollback to improve QC and window functions before proceeding to panel-wide screening.

### *IN PROGRESS*
### Phase 4 — Pipeline Improvement (08–09)

**`08_improved_QC`** addresses rogue transcripts (non-tissue-associated points that inflate the bounding box). Uses DBSCAN density-based filtering per strip to classify tissue-core vs. noise transcripts. Includes before/after visualisation, parameter sensitivity analysis, and gene-level bias checks. Outputs `s1_all_strips_cleaned.parquet` with a noise flag column (no data deletion). [See scaffold below.](#08-improved-qc-scaffold)

**`09_window_functions`** replaces the rectangular bounding box with geometry-aware observation windows (convex hull, then alpha shape). Requires updating the edge correction in `bivariate_k()` for polygonal windows (Monte Carlo arc-in-polygon or analytical intersection). Includes regression tests against the rectangular case and re-running the three-part control framework with the new windows. [See scaffold below.](#09-window-functions-scaffold)

## 08 Improved QC — Scaffold

1. **Visualise the problem.** Load `s1_all_strips.parquet`, plot raw scatter per FOV×strip, characterise rogue transcript patterns (isolated points, inter-strip scatter, FOV-edge artefacts).
2. **Density-based filtering.** DBSCAN per strip on 2D coordinates. `eps` derived from nearest-neighbour distance distribution within tissue (e.g., a quantile of the kNN distances). `min_samples` tuned per-strip. Alternative: simple kNN distance threshold.
3. **Before/after comparison.** Scatter overlays, transcript counts, bounding-box change. Check for gene-level bias in removed transcripts (critical: low-abundance signalling genes must not be disproportionately filtered).
4. **Sensitivity analysis.** Vary `eps` ± 20%, report transcript loss. Stable results strengthen the case.
5. **Save output.** `s1_all_strips_cleaned.parquet` with `is_noise` flag column. Document all parameters.

## 09 Window Functions — Scaffold

1. **Implement convex hull window** via `scipy.spatial.ConvexHull`. Refactor `get_window()` to return a window object with type and geometry metadata.
2. **Implement alpha shape window** via `alphashape` or Delaunay triangulation. Tune alpha parameter per strip.
3. **Polygonal edge correction.** Monte Carlo approach first (sample points on circle, count fraction inside polygon via `shapely` point-in-polygon). Validate against `_fraction_inside_rect()` on rectangular test case.
4. **Update `bivariate_k()`** to dispatch on window type. Polygon area via `shapely` for intensity λ.
5. **Regression test.** Rectangular polygon matching old bounding box → L(r) curves must match.
6. **Synthetic CSR validation.** Uniform points inside alpha-shape polygon → L(r) ≈ 0.
7. **Re-run controls.** KRT8×KRT18, MALAT1×KRT18, KRT8×SCGB3A1 with cleaned data + polygon windows. Side-by-side comparison with old rectangular results.

## Control Framework

| Control | Gene Pair | Expected Result | Confirmed |
|---|---|---|---|
| **Positive** | KRT8 × KRT18 | L(r) above envelope (co-expressed epithelial keratins) | Yes (strip 2) |
| **Negative 1** | MALAT1 × KRT18 | L(r) inside envelope (ubiquitous housekeeping × epithelial) | Yes |
| **Negative 2** | KRT8 × SCGB3A1 | L(r) inside envelope (different epithelial lineages) | Yes |

The permutation null shuffles gene labels while preserving spatial locations. It tests whether co-localisation is gene-specific or merely a consequence of shared tissue structure. This is distinct from a Poisson/CSR null, which would flag all gene pairs as significant (transcripts cluster because they come from cells, not because of gene identity).

## Milestones

| Date | Milestone | Status |
|---|---|---|
| Mar 31 | End-to-end pipeline on S1 | Done - crude but functional |
| Apr 30 | Looping across datasets + network development | |
| May 15 | FEATURE FREEZE - background + methods written | |
| May 31 | Packaging, documentation, results finalised | |
| Jun 3 | CODE FROZEN | |
| Jul 16 | SUBMISSION | |

## Requirements

See `requirements.txt`. Core dependencies include: `numpy`, `scipy`, `pandas`, `matplotlib`, `scikit-learn`, `spatialdata`, `pyarrow`, `shapely`, `alphashape`, `liana`.

## Acknowledgements

Dataset provided by the King's College London. Ligand-receptor annotations from [CellChatDB](https://github.com/sqjin/CellChat) (Jin et al., *Nature Communications*, 2021), accessed via [LIANA](https://github.com/saezlab/liana-py) (Dimitrov et al., *Nature Communications*, 2022).
