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