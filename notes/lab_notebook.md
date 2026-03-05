### 02/03/26
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

### 05/03/26
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