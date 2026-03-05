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