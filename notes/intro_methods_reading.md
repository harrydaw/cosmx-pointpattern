# Introduction & Methods Overview — Reading List

*Working document for dissertation write-up. Last updated 13/04/26.*

---

## Introduction: What to Cover

### 1. Spatial transcriptomics — the field
Single-cell transcriptomics (scRNA-seq) revolutionised our understanding of cell types but discards spatial context by dissociating tissue. Spatial transcriptomics retains the physical location of each transcript within the tissue section, enabling spatial gene expression analysis at or near single-cell resolution.

Key platforms to describe: **NanoString CosMx SMI** (used here), 10x Genomics Xenium, MERFISH, Slide-seq. Compare: targeted panel vs. transcriptome-wide; single-molecule resolution vs. pixel-resolution.

### 2. The segmentation bottleneck
Most spatial transcriptomics analysis pipelines rely on cell segmentation to assign transcripts to cells before any downstream analysis. Segmentation depends on nuclear staining (DAPI) or membrane markers. In low-quality or clinically-derived tissue, staining can be patchy, leading to high transcript dropout rates (~42% unassigned in this dataset) and unreliable compartment labels.

This motivates segmentation-free approaches: analyse spatial relationships between transcript species directly, without the intermediate cell-assignment step.

### 3. Point pattern analysis
A spatial point pattern is a set of events (here: transcript detections) with 2D coordinates observed within a study region. The core question: are the events of two types (gene A, gene B) more spatially associated than expected by chance?

Introduce: intensity (λ), complete spatial randomness (CSR / homogeneous Poisson process), first-order vs. second-order properties, the K-function as a second-order summary statistic.

### 4. Ripley's K-function
K(r) is the expected number of extra events within distance r of a typical event, normalised by intensity. Under CSR, K(r) = πr². Co-localisation → K(r) > πr²; inhibition → K(r) < πr².

The bivariate (cross-type) K-function K_AB(r) measures the spatial relationship between two different event types (gene A and gene B). L(r) = √(K(r)/π) − r is the variance-stabilised transform, centred at zero under CSR.

Edge correction: transcript detections near the tissue boundary have fewer observable neighbours on one side, biasing K downward. Must correct. Describe Ripley's isotropic correction for rectangular windows, and the Shapely polygon intersection extension used here.

Permutation envelope: the appropriate null for tissue data is label permutation (not Poisson simulation), which preserves spatial intensity structure while destroying gene-specific associations.

### 5. Cell-cell communication
Cells communicate via ligand-receptor (L-R) interactions. If a ligand (secreted by cell type A) and its receptor (expressed on cell type B) are spatially co-localised, this is evidence of active signalling between those cells. Spatial co-localisation of transcript species is a proxy for this.

Reference: CellChatDB (curated human L-R pairs), LIANA (Python interface).

### 6. Network analysis and modularity
A co-localisation matrix (one entry per L-R pair, per strip) can be encoded as a weighted graph. Community detection (modularity maximisation) finds groups of genes that are mutually co-localised — putative co-regulated signalling programmes.

Modularity Q measures the density of edges within communities relative to a random graph null. Louvain and Leiden algorithms optimise Q greedily.

### 7. Lung biology context (brief)
RSV (respiratory syncytial virus) infection in the context of severe asthma. Tissue contains a mix of epithelial, immune, and stromal cells. Three-strip layout: control tissue / RSV-infected tissue / control tissue within each slide section. Goal: identify spatial signalling programmes that are specific to the infected strip.

---

## Methods: What Has Been Implemented

### Data loading and preprocessing
CosMx data loaded via **SpatialData** (Marconato et al., 2024) from `.zarr` archives. Global pixel coordinates (`x_global_px`, `y_global_px`) used throughout — local FOV coordinates discarded as meaningless across stitched FOVs. 710,751 total transcripts in S1 across 13 FOVs; 6 FOVs retained after visual QC (~400k transcripts).

### Strip assignment
Tissue strips (control / infected / control) identified by fitting a Gaussian Mixture Model (GMM, 3 components) to the x-coordinate distribution within each FOV. Implemented in `02b_strip_assignment_all_fovs.ipynb`. Per-FOV decisions recorded in `fov_config` dictionary for reproducibility. Output: `s1_all_strips.parquet`.

### Rogue transcript removal (DBSCAN)
Transcripts far outside the tissue boundary inflate the observation window and bias the K-function. Removed using DBSCAN (density-based spatial clustering) per strip. `eps` derived adaptively from the 97th percentile of 1-nearest-neighbour distances, clipped to [20, 30] px. 27,313 transcripts flagged as noise (6.8%), retained in dataset with `is_noise` flag. Implemented in `08_improved_QC.ipynb`.

### Bivariate Ripley's K-function
Core estimator implemented from scratch in `03_K_function.ipynb`. Per-pair enumeration via `scipy.spatial.cKDTree`. Edge correction: Ripley's isotropic arc-fraction method for rectangular windows (corrected from an initial additive-fraction bug; validated: max |L(r)| < 13 on synthetic CSR). L-transform and permutation envelope (99 label-shuffling permutations).

### Convex hull observation windows
Rectangular bounding boxes include large amounts of empty space (hull = 27–43% of bounding box area in these strips). Replaced with convex hull polygons (Shapely `MultiPoint.convex_hull`) per strip in `09_improved_windows_and_edge_correction.ipynb`. Edge correction updated to use Shapely polygon intersection (`fraction_inside_hull`). Intensity λ_B now estimated over hull area rather than bounding box area.

**Known limitation:** The Shapely intersection edge correction introduces a small systematic negative bias in absolute K(r). This bias cancels in the permutation test (hull is fixed; both observed and permuted use the same biased estimator). Absolute L(r) values are not interpreted directly.

### L-R pair feasibility
146 L-R pairs from CellChatDB viable across all three strips (both ligand and receptor ≥ 50 transcripts per strip). Identified in `06_LR_checks.ipynb`.

---

## Reading List

### Spatial statistics — foundations

| Reference | Why you need it |
|-----------|----------------|
| Ripley, B.D. (1976). "The second-order analysis of stationary point processes." *Journal of Applied Probability*, 13(2), 255–266. | Original K-function paper. Cite when introducing K(r). |
| Diggle, P.J. (2003). *Statistical Analysis of Spatial Point Patterns* (2nd ed.). Arnold. | The definitive textbook. Read Ch. 3 (summary statistics), Ch. 7 (bivariate). |
| Baddeley, A., Rubak, E., & Turner, R. (2015). *Spatial Point Patterns: Methodology and Applications with R*. CRC Press. | Modern reference with spatstat. Good for methods justification. |
| Ripley, B.D. (1988). *Statistical Inference for Spatial Processes*. Cambridge University Press. | Edge correction details; permutation tests for point patterns. |

### Spatial transcriptomics — platforms and methods

| Reference | Why you need it |
|-----------|----------------|
| He, S. et al. (2022). "High-plex imaging of RNA, proteins, and metabolites at subcellular resolution in fixed tissue by spatial molecular imaging." *Nature Biotechnology*, 40, 1794–1806. | CosMx SMI platform paper. Cite when describing the data. |
| Moses, L. & Pachter, L. (2022). "Museum of spatial transcriptomics." *Nature Methods*, 19, 534–546. | Comprehensive review of all ST platforms. Good for intro context. |
| Marconato, L. et al. (2024). "SpatialData: an open and universal data framework for spatial omics." *Nature Methods*, 21, 1106–1116. | SpatialData framework used for data loading. Cite in methods. |
| Eng, C.L. et al. (2019). "Transcriptome-scale super-resolved imaging in tissues by RNA seqFISH+." *Nature*, 568, 235–239. | seqFISH+: useful contrast platform for introduction. |

### Cell-cell communication

| Reference | Why you need it |
|-----------|----------------|
| Jin, S. et al. (2021). "Inference and analysis of cell-cell communication using CellChat." *Nature Communications*, 12, 1088. | CellChatDB source. Cite when describing L-R pair database. |
| Dimitrov, D. et al. (2022). "Comparison of methods and resources for cell-cell communication inference from single-cell RNA-seq data." *Nature Communications*, 13, 3224. | LIANA paper and method comparison. Cite when describing L-R database access. |
| Armingol, E. et al. (2021). "Deciphering cell–cell interactions and communication from gene expression." *Nature Reviews Genetics*, 22, 71–88. | Review of cell-cell communication inference methods. Good for intro. |

### Network analysis and modularity

| Reference | Why you need it |
|-----------|----------------|
| Blondel, V.D. et al. (2008). "Fast unfolding of communities in large networks." *Journal of Statistical Mechanics*, 2008(10), P10008. | Louvain algorithm. Cite if using Louvain community detection. |
| Traag, V.A., Waltman, L. & van Eck, N.J. (2019). "From Louvain to Leiden: guaranteeing well-connected communities." *Scientific Reports*, 9, 5233. | Leiden algorithm (preferred over Louvain). Cite if using Leiden. |
| Newman, M.E.J. & Girvan, M. (2004). "Finding and evaluating community structure in networks." *Physical Review E*, 69(2), 026113. | Original modularity Q paper. Cite when introducing modularity. |
| Fortunato, S. (2010). "Community detection in graphs." *Physics Reports*, 486(3–5), 75–174. | Comprehensive review of community detection. Good background. |

### Python ecosystem

| Reference | Why you need it |
|-----------|----------------|
| Hagberg, A.A., Schult, D.A. & Swart, P.J. (2008). "Exploring network structure, dynamics, and function using NetworkX." *Proceedings of SciPy 2008*. | NetworkX citation. Required when using networkx. |
| Pedregosa, F. et al. (2011). "Scikit-learn: Machine learning in Python." *JMLR*, 12, 2825–2830. | scikit-learn citation. Required for DBSCAN. |
| Ester, M. et al. (1996). "A density-based algorithm for discovering clusters in large spatial databases with noise." *KDD-96 Proceedings*, 226–231. | Original DBSCAN paper. Cite in QC methods section. |
| Virtanen, P. et al. (2020). "SciPy 1.0: Fundamental algorithms for scientific computing in Python." *Nature Methods*, 17, 261–272. | SciPy citation. Required for cKDTree. |

---

## Key Terms to Define in Methods

- **Spatial point pattern** — realisation of a point process; set of event locations within a study region
- **Intensity (λ)** — expected number of events per unit area; estimated as n / |W| where |W| is the window area
- **Complete Spatial Randomness (CSR)** — homogeneous Poisson process; events are independent and uniformly distributed
- **Observation window (W)** — the region within which events are observable; correct specification is critical for unbiased K estimation
- **Edge correction** — adjustment for the fact that events near the window boundary have truncated neighbourhoods
- **Permutation envelope** — pointwise min/max of L(r) across label-shuffled realisations; used as the null reference in lieu of CSR
- **Convex hull** — smallest convex polygon containing all points; used here as the observation window
