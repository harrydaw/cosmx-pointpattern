### 02/03/26 - EDA
**Data exploration and narrowing**
- Today I set up my github structure and venv
- I then looked into the Varsha_1234.zarr and analysed the structure and counts
- Rough segmentation images show S1 to be the most dense and useable
- Focused on the updated_stitched_S1.zarr and then visualised by the 13 FOVs
- 2 extra columns are added here for "realtive" x and y given the new stitched S1 image
- Visualised X coordinate distribution to find peaks and tissue boundaries
- Separated FOV3 by X bounds to give the 3 tissues
#### Tomorrow:
    - Split the 3 tissues into separate data
    - Analyse gene count differences to double check infected vs. control
    - First point pattern analysis

### 05/03/26 Narrowed data + Functions  
**Narrowed data exploration and validation (FOV3, Slide 1)**
- Loaded the fov3_strip.parquet and confirmed structure
- Set up gene x strip matrix
- Determined minimum transcript count of 50 per gene
    - "sparse patterns produce high variance K estimates"
- Potential data limitation wiith the gmm strip assignment 
    - strip 1 is really two tissue sections that were captured under a single GMM component
    - Shouldn't matter because the infected (strip 2 ) vs controls (1 and 3) doesn't invalidate the K-function
    - Need to make sure everything works across other strips and data though for the future package
- 16 genes usable across strip_2 and strip_3
- KRT8 x KRT18 selected as positive control pair
- Scatter confirms fine-scale spatial interleaving in strip_2
- 
- Defined all of my core functions:
    get_coords()
    get_window()
    bivariate_k()
    k_to_l()
    compute_envelope()

- Some issues with massive L values (should be in the 10's and I was getting 100's or 1000's)
- Potentially a global_px vs global_px_transformed issue?
- Need to look back at the gmm histogram and ensure I used the correct inputs all the way through.

#### Bug identified: edge correction in bivariate_k()
- Symptom: L(r) values in the hundreds/thousands for both real data (KRT8×KRT18) and synthetic CSR data where L(r) should ≈ 0.
- Root cause: the arc_in() helper function was computing the fraction of circumference inside the window relative to each edge independently, then summing four per-edge fractions. For an interior point far from all edges, each edge contributes ≈0.5 (half the circle is "on that side"), so the sum ≈ 2.0, clipped to 1.0. This meant interior points got correct weights (1.0), but edge-adjacent points got systematically wrong fractions — the additive decomposition doesn't account for the geometry correctly.
- The correct approach: compute the total arc angle outside the window by accumulating 2·arccos(d_edge / d_ij) for each edge where d_edge < d_ij, then compute frac_inside = (2π - angle_outside) / 2π. This starts from the assumption that the full circle is inside (fraction = 1) and subtracts the protruding portions, rather than trying to sum "inside" contributions from four independent edges.
- This is Ripley's isotropic edge correction for rectangular windows — the same default used by spatstat::Kcross() in R.
- Minor approximation: the angle subtraction double-counts corner regions where two edges both clip the circle. For points very near corners this slightly overcorrects (inflates the weight). Negligible for practical purposes and consistent with standard implementations.
- A safety floor of max(frac_inside, 0.01) prevents division-by-zero for points exactly on the boundary.
- Validation: after the fix, synthetic CSR test (two independent uniform patterns, n=100 each) gives max |L(r)| < 13 for both a 1000×1000 square window and a 1200×4200 elongated window matching strip geometry. Before the fix, the same test gave L(r) ≈ 195–445.
- The _fraction_inside_rect helper uses a leading underscore to mark it as a private/internal function — not part of the public API when the code is packaged.

- Positive control (KRT8xKRT18) successful, co-localisation seen in both strips, significant in strip_2

**Important distinction is permutation vs. simulation. Permutation (for envelope development) retains the spatial locations of the trancripts but destroys the gene-specific labels, and then sees if the co-localisation is linked specifically to the gene label, or just generally is a result of the spatial distribution of the points. Simulation would generate new points in random space, which would mean our result would always come up as more significant, because the transcripts will co-localise based on the fact that they are from the same cells.**

#### NEXT:
- Negative control, run on two unrelated genes that have no reason to co-occur and confirm L(r) values of ~0.

### 09/03/26 Validation via negative control
#### Content check:
**The permutation null preserves the spatial intensity landscape of all transcripts and asks whether the observed co-localisation is specific to the gene labels, or merely a consequence of shared tissue structure. A Poisson simulation null tests departure from complete spatial randomness, which is inappropriate here because transcript clustering is a biological given — virtually all pairs would exceed it, making it uninformative for discriminating between gene pairs.**

### 10/03/26 More control experimentation 
- MALAT1 × KRT18 negative control initially appeared problematic - raw L(r) values (~120–205) were similar in magnitude to the KRT8 × KRT18 positive control (~131–225). - Realised this is expected: MALAT1 is ubiquitous (present in all cells), so it's spatially close to everything. High raw L(r) reflects shared tissue structure, not gene-specific interaction.
- Key conceptual clarification: raw L(r) magnitude vs envelope significance. Under the permutation null, the relevant test is whether observed L(r) exceeds the permutation envelope, not whether raw L(r) is near zero. L(r) ~0 is the expectation under a CSR/Poisson null, which is a different question. The permutation null asks: "given the overall spatial distribution of transcripts, is there additional gene-specific co-localisation?" This distinction is critical for the viva.
- Added plot_diagnostics() function to notebook 03 — scatter plot of both transcript patterns with bounding-box overlay. Accepts an ax parameter for composable multi-panel layouts. Controlled by a DIAGNOSTICS = True/False flag in the analysis loop. Will live in spatialco/plotting.py at packaging stage.
Selected KRT8 × SCGB3A1 as a second negative control pair. Biological rationale: KRT8 marks simple epithelial cells (alveolar), SCGB3A1 marks secretory airway cells (club/goblet). Different epithelial lineages. However, diagnostic scatter plots revealed both transcript types are spatially interleaved in these strips — the expected spatial separation doesn't hold in this FOV.
- Raw L(r) is high for all gene pairs tested due to narrow strip geometry forcing spatial proximity. The permutation envelope is therefore not just a significance test but the only way to distinguish gene-specific co-localisation from tissue-structure confounding.

**Control framework now three-part:**
    Positive: KRT8 × KRT18 — high L(r), above envelope (gene-specific)
    Negative 1: MALAT1 × KRT18 — high L(r), inside envelope (tissue structure only)
    Negative 2: KRT8 × SCGB3A1 — high L(r), inside envelope (tissue structure only)

- Potentially limited atm with a narrow FOV, single strip at a time, and low transcript counts
- Will reevaluate results when I loop this across more data

- Added run_pair_analysis() function 
    - Just runs all of my custom functions on specific inputs for 1 gene pair on a strip dataframe
    - Kept the strip by strip loop outside of this to allow me to call the function for per FOV/strip/plate... analysis

**L-R panel check**
- Loaded CellChatDB (Jin et al., Nat Comms 2021) via the liana Python package (Dimitrov et al., Nat Comms 2022). CellChatDB contains 1,912 curated human ligand-receptor pairs.
- Cross-referenced with the 1,123 unique gene targets in the CosMx panel. Many L-R genes are present in the panel but none reach the 50-transcript threshold in any single strip of FOV3. Highest counts: VEGFA (41 in strip_2, 32 in strip_3), CCL21 (48 in strip_3).
- Signalling molecules are expressed at much lower per-FOV levels than structural/housekeeping genes. This is a genuine data characteristic, not a pipeline issue.
- Key decision: Single-FOV L-R analysis is not feasible at the 50-transcript threshold. Multi-FOV pooling (aggregating results across FOVs, not concatenating coordinates) is required to either boost per-gene counts or to lower the threshold with replication across FOVs as statistical support.

**Multi-FOV strip assignment**
- Loaded full S1 zarr (updated_stitched_S1.zarr) via spatialdata. 13 FOVs, 710,751 total transcripts.
- Created notebook 02b_strip_assignment_all_fovs.ipynb. Two-pass approach:

    Pass 1: Plotted raw spatial scatter + x-coordinate histogram side-by-side for all 13 FOVs.
    Pass 2: Recorded decisions in a fov_config dictionary (FOV ID → n_strips, which strips to keep). Ran GMM strip assignment for all usable FOVs in a single processing loop. Visualised each FOV's GMM result (histogram + coloured scatter) for confirmation and sanity check

    FOVs excluded: 1, 2, 3, 6, 7, 12, 13 [7] — generally due to missing strips, additional structures or bad strip mapping
    FOVs retained: 4, 5, 8, 9, 10, 11 [6] — all with 3-strip structure, all strips kept.

- FOV 3 excluded from this batch as it already has a separate fov3_strips.parquet from earlier work and has a fourth tissue section that is undetermined.

- Final overview plot confirmed consistent strip colouring across all retained FOVs.
- Saved s1_all_strips.parquet - ~400,000 transcripts across 6 FOVs. Had to clear attrs metadata (SpatialData Affine transform objects) before parquet serialisation.
- fov_config dictionary serves as a reproducible record of all FOV selection decisions for the Methods section.

**Design decisions**
- Approach A (per-FOV analysis, not coordinate pooling): Run run_pair_analysis() separately per FOV per strip, then aggregate L(r) curves. This is methodologically cleaner than pooling coordinates across physically non-contiguous FOVs, which would violate the stationarity assumption of the K-function. Also provides replication (one L(r) curve per FOV) rather than a single pooled estimate.

- run_pair_analysis() function added to notebook 03: Wraps the entire analysis loop (coord extraction → K/L computation → envelope) into a single function call. Returns a results dict. This was motivated by having run the same pattern three times (KRT8×KRT18, MALAT1×KRT18, KRT8×SCGB3A1) and prepares for the multi-FOV × multi-pair screening.

- Strip labels (strip_1 = leftmost by x-coordinate) are assigned consistently via GMM ordering

#### Up next:
- Run LR checks across all 6 successful FOVs and rerun previous checks


### 11/03/26 - Expanded L-R validation and controls
- Extended L-R feasibility check to all_labelled (400,238 transcripts, 6 FOVs). Threshold remains n≥50 per gene per strip. Results: 199/203 L-R pairs viable in ≥1 strip, 146 viable across all three strips. Single-FOV sparsity was the limiting factor. Pair selection for notebook 08 to be driven by biological hypothesis, not data availability

### 16/03/26
- Following on from last week, I have clearly been able to incorporate accurate statistical anaylsis and generate a crude end-to-end pipeline but there are some clear violations of the basic assumptions of the model
- Before I expand to wrap the gene x gene checks to search over a whole panel of interactions, I need to go back and tidy up the QC and data handling to be able to get stronger results
- It's very hard to explain the figures as it stands, and hard to  make any conclusions, so I am reluctant to pass them onto my netowork analysis/modularity-maximisation steps before I sort out step 1
- Current plan is to incorporate more options for window generations like spatstat
- Current I just assign columns based on X coordinate distributions, but I want to expand to have polygonal and binary mask options
- Better cleaning of the sparse regions + better tissue mapping would give the results a lot more weight, and mean that (at the very least) the homogeneity assumption wouldn't be so clearly violated
- ^^ This is the goal for the week
    - Get to the point where my figures and K and L-function are more explainable and the controls show what they are supposed to in a more obvious way


### 30/03/26 - Documentation reset and planning

**Project review and scaffolding**
- Reviewed full project state: weeks 1–5 summary, all notebook purposes, function inventory, outstanding assumptions violations
- Core issue remains: bounding-box window inflates |W|, distorts λ, and makes L(r) results hard to interpret. Rogue transcripts outside tissue compound this by pulling the box further from true tissue geometry. Must fix QC (08) before fixing windows (09) - cleaning is a prerequisite for window fitting.

**Scaffolded notebooks 08 and 09:**
- 08_improved_QC: density-based rogue transcript removal via DBSCAN per strip. Key steps: visualise the noise, filter, before/after comparison, sensitivity check, save with `is_noise` flag column. Must check for gene-level bias in removed transcripts (if low-abundance signalling genes are disproportionately removed, we've introduced systematic bias into the L-R screening).
- 09_window_functions: replace rectangular bounding box with convex hull (first) then alpha shape, then push for a binary mask via kernel density. Requires new edge correction for polygonal windows - Monte Carlo arc-in-polygon approach first, validate against rectangular ground truth, then synthetic CSR on polygon. Re-run three-part control framework to compare old vs new windows. Ideally upgrade to arc-intersection analytics for optimum accuracy. 

**README overhaul**
- Rewrote README.md with full project context, notebook guide, function reference, updated file tree, control framework table, and 08/09 scaffolds
- Old README was a placeholder from week 1 — new version reflects actual project state and is suitable for supervisor/examiner review

#### Up next (tomorrow):
- Begin 08_improved_QC implementation
- Step 1: load s1_all_strips.parquet, raw scatter plots per FOV×strip to characterise the noise
- Step 2: DBSCAN implementation - need to decide on eps derivation strategy (kNN distance quantile vs. manual)

### 31/03/26 - Notebook 08: Rogue transcript QC via DBSCAN
**The core problem**
Rogue transcripts far from the tissue inflate this box, which violates the homogeneity assumption of K/L and makes my edge correction unreliable - the window includes empty space that isn't really observation area, deflating my intensity estimate λ and biasing K(r) upward. Ensuring that the window is tight to the tissue border, and includes only informative gene transcripts will better produce interpretable results.

- Loaded s1_all_strips.parquet (400,238 transcripts, 6 FOVs × 3 strips)
- Visualised raw scatter — noise is diffuse/uniform across all FOVs, not concentrated
- Computed 1-NN distance distributions per strip. Medians 1.4–2.0 px (tissue bulk), long tail to ~150 px (noise)
- Tested fixed eps (15, 30, 45) — eps=30 performed best across strip densities
- Tested adaptive eps via p95, p97, p99 of 1-NN. p95 too aggressive (shredded tissue edges on dense strips), p99 too lenient (missed noise on sparse strips), p97 best compromise
- Tested k=5 NN for eps derivation (matching min_samples) — eps inflated to 70+, <1% removal. Rejected.
- Final method: adaptive 1-NN p97, clipped to [20, 30]. Floor prevents over-trimming dense strips, ceiling prevents leniency on sparse strips
- Built `dbscan_noise_filter()` — generalisable function with tunable percentile, floor, ceiling, group_cols, diagnostic toggle. Flags points (`is_noise`) without removing. Records `eps_raw` and `eps_used` per group for audit trail.
- Result: 27,313 flagged as noise (6.8%), 372,925 retained
- Saved `s1_all_strips_noise_flagged.parquet`

#### Tomorrow:
- Begin notebook 09: improved window bounding (polygonal/convex hull) on cleaned data
- Re-run controls on filtered data to check impact on L(r) interpretability

### 13/04/26 — Notebook 09: Convex hull windows and Shapely edge correction

**Implemented three core functions:**
- `get_convex_hull(df)` — returns Shapely Polygon from all transcript coords in a strip. Raises on <3 points or collinear degenerate case.
- `fraction_inside_hull(point, r, hull, resolution=64)` — approximates fraction of a disc of radius r centred at point that lies inside the hull polygon via Shapely intersection. Fast-path for interior points. Denominator is disc.area (polygon approximation) for internal consistency. Resolution parameterised.
- `bivariate_k_hull(coords_a, coords_b, r_vals, hull, resolution=64)` — replaces rectangular arc correction with `fraction_inside_hull` per pair; normalises by hull.area rather than bounding box area.
- Supporting: `compute_envelope_hull`, `run_pair_analysis_hull` — structural mirrors of nb03 equivalents, hull held fixed across permutations.

**Window comparison (09_window_comparison.png):**
- Hull = 43% of bounding box for strips 1 and 3, 27% for strip 2 (infected). 
- Large improvement: the old rectangular window was including 57–73% empty space in the intensity estimate and edge correction. The convex hull tightly follows the actual transcript distribution.

**Validation 1 — Rectangular regression (09_validation_rect_regression.png):**
- Hull and arc methods broadly track each other on K(r) but diverge on L(r): arc method sits near zero (expected for CSR), hull method gives systematically negative L (≈ −4 to −5). Both methods applied to identical synthetic CSR data in a rectangle — they should agree. Indicates the hull correction under-corrects for edge effects, under-weighting boundary pairs and pulling K down.

**Validation 2 — CSR on hexagon (09_validation_csr_polygon.png):**
- Large systematic negative bias: L(r) drops to ≈ −20 at r=300. This is the priority problem. Pattern is consistent with `fraction_inside_hull` over-estimating the fraction inside for boundary points (weight too small → K too small). Suspected cause: the Shapely `contains()` fast-path may be firing when it shouldn't, or the polygon intersection is not accurate enough at the disc/hull boundary. Needs investigation before controls can be trusted.

**Controls (run on noise-cleaned data, convex hull window):**
- KRT8 × KRT18: L(r) far above permutation envelope in all three strips, as expected for a positive control. ✓
- MALAT1 × KRT18: elevated raw L(r) but observed tracks at or within the envelope — consistent with tissue-structure confounding, same pattern as notebook 07.
- KRT8 × SCGB3A1: same as above — elevated L but not clearly above envelope.
- Negative controls are consistent with the permutation null doing its job, but hard to interpret confidently while the CSR validation is failing.

**Decision on CSR calibration bias:**
After review, decided not to fix the absolute calibration and instead document it as a known limitation. Rationale: the permutation test is robust to consistent estimator bias — both observed K and all permuted realisations use the same fixed hull and same edge correction, so the bias cancels in the relative comparison. Controls confirm the test is discriminating correctly. A limitation note has been added to the bottom of nb09. Absolute L(r) values are not to be interpreted; only envelope-relative significance is used.

#### Tomorrow:
- Refactor `get_window()` (in nb03) to accept `window_type` argument: `'rect'` (existing), `'hull'` (wraps get_convex_hull), `'custom'` (imports GeoJSON polygon via new `load_custom_window()` helper)
- Build `10_window_comparison_all_fovs.ipynb`: clean three-way window comparison across all 6 S1 FOVs × 3 strips, side-by-side L(r) curves for rect vs hull on KRT8 x KRT18
- After that: proceed directly to L-R screening (nb11) and network construction (nb12)