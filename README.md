# STIPPLE

**S**patial, **T**ranscript-level **I**nference of **P**airwise **P**roximity in **L**igand-receptor **E**xpression

A segmentation-free pipeline that screens ligand-receptor pairs for spatial co-localisation **directly on raw transcript coordinates**, using a bivariate, edge-corrected Ripley's K analysis with a label-permutation envelope. No cell segmentation is required at any stage.

*MSc Bioinformatics dissertation, King's College London, 2026.*

> Research codebase, not a maintained library. The `stipple` package holds the reusable
> statistical tools; the notebooks reproduce the dissertation analysis end to end.

---

## Why this exists

Cell segmentation is a fragile prerequisite for almost every downstream analysis in imaging-based spatial transcriptomics. On the CosMx sample analysed here, vendor segmentation left **43.6% of 702,873 transcripts unassigned** to any cell, which makes segmentation-dependent ligand-receptor tools unreliable on exactly the kind of clinical tissue where they are most needed.

STIPPLE sidesteps the problem. Each gene becomes a planar point pattern, and a bivariate Ripley's K with a label-permutation envelope asks whether two genes' transcripts co-localise more than expected under a null that preserves tissue structure and each gene's abundance. Pairs are ranked by a continuous standardised effect size (SES), and the screen scales to panel-level via an HPC job array.

## Install

```bash
git clone https://github.com/harrydaw/cosmx-pointpattern.git
cd cosmx-pointpattern
python -m venv .venv && source .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install -e .                                     # installs the `stipple` package + core deps
```

The core tools depend only on numpy, pandas, scipy, shapely, scikit-learn, matplotlib and networkx. To also load raw CosMx zarr stores and build the CellChatDB panel, install the optional extra: `pip install -e ".[data]"`. To reproduce the exact dissertation environment, use the pinned `requirements.txt` lock instead.

## Quick start

`run_pair_analysis` is the single entry point. Give it a transcript table (a `pandas` DataFrame with a `target` gene column and `x_global_px` / `y_global_px` coordinate columns), two gene names and a radius grid, and it returns the observed L(r), its permutation envelope and everything needed to compute the SES.

```python
import numpy as np, pandas as pd, stipple

# a minimal example: two synthetic, co-localised genes
rng = np.random.default_rng(0)
a = rng.uniform(0, 1000, size=(200, 2))
b = a[:150] + rng.normal(0, 8, size=(150, 2))         # gene B sits near gene A
df = pd.concat([
    pd.DataFrame({"target": "GENEA", "x_global_px": a[:, 0], "y_global_px": a[:, 1]}),
    pd.DataFrame({"target": "GENEB", "x_global_px": b[:, 0], "y_global_px": b[:, 1]}),
], ignore_index=True)

r_vals = np.linspace(5, 100, 20)
res = stipple.run_pair_analysis(df, "GENEA", "GENEB", r_vals, window_type="hull", n_sim=99)

mid  = 0.5 * (res["l_hi"] + res["l_lo"])              # envelope midline
half = 0.5 * (res["l_hi"] - res["l_lo"])              # envelope half-width
ses  = (res["l_obs"] - mid) / half
print("peak SES:", round(np.nanmax(ses), 2))          # strongly positive for a co-localised pair
```

## Running it on your own data

STIPPLE needs only a table of transcript points. To apply it to a new dataset:

1. **Shape your data** as a DataFrame with a gene-label column (`target` by default; override with `x_col` / `y_col` on any function) and two coordinate columns in a consistent unit.
2. **Define an observation window** with `get_window(df, window_type="hull")`. A concave hull is recommended so the window follows the tissue outline rather than its bounding box; `"rect"` and `"custom"` (a saved polygon) are also available.
3. **Screen a pair** with `run_pair_analysis(...)`, or call `bivariate_k` / `compute_envelope` directly for full control.
4. **Rank pairs by peak SES** (the maximum of the SES over radius). Only the envelope-relative position of the curve is interpreted, never the absolute L(r) (see the edge-correction note below).

The preprocessing helpers (`add_rotated_coords`, `reassign_strips_gmm`, `dbscan_noise_filter`, `apply_cleanup`) are specific to the strip-arranged CosMx layout used here and are optional for other datasets.

## Where this sits in the existing toolset

| Tool | Input | Output | Segmentation needed? |
|---|---|---|---|
| Bento (Mah et al. 2024) | transcripts + cell/nuclear masks | subcellular localisation patterns | Yes (mandatory) |
| Squidpy `gr.ripley` | cell coordinates + cluster labels | cluster-level co-localisation | Yes (cells exist) |
| CellChat / LIANA | cell-type expression matrix | inferred L-R signalling | Yes (cells + clusters) |
| Giotto (spatial L-R) | segmented cells + types | spatially informed L-R | Yes |
| FICTURE (Si et al. 2024) | transcripts | pixel-level spatial domains | No, but not pairwise L-R |
| spatstat (R) | generic 2D point pattern | K, L, etc. | N/A (generic stats) |
| **STIPPLE** | **raw transcripts only** | **pairwise L-R co-localisation (SES)** | **No** |

The statistic itself is established prior art; the contribution is the application and the integrated pipeline, running the cross-type Ripley's K on transcripts, with a tissue-fitted window and panel-scale screening.

## Pipeline

1. **Load and QC** raw CosMx transcript tables (global x/y, gene labels, FOV assignments).
2. **Assign tissue strips** per FOV via a Gaussian mixture on the PCA-rotated coordinate (control / infected / control).
3. **Clean off-tissue transcripts** with a density filter (DBSCAN) plus a manual exclusion list, so windows reflect real tissue geometry.
4. **Fit an observation window** to each strip (concave hull by default) for edge correction.
5. **Compute bivariate Ripley's K and the L-transform** per pair, with Shapely polygon-aware edge correction.
6. **Build a label-permutation envelope**, then a continuous SES (deviation from the envelope midline in half-widths).
7. **Screen ligand-receptor pairs** (CellChatDB intersection) across strips via a SLURM array on KCL CREATE.
8. **Assemble per-strip networks** from the surviving pairs.

## The data

| Property | Value |
|---|---|
| Platform | NanoString CosMx Spatial Molecular Imager (SMI) |
| Tissue | RSV-infected, severe-asthmatic human lung |
| Layout | Slide S1, three tissue strips: control / infected / control |
| Transcripts | ~2.66M across four slides; 702,873 on S1 (12 FOVs) after outlier exclusion |
| Gene panel | 1,000 targets + 78 SystemControl probes |
| Coordinates | Global pixel positions (`x_global_px`, `y_global_px`); spatially uncalibrated |

The raw data are patient-derived and are **not** included in this repository (`data/` is gitignored).

## Repository layout

```
cosmx-pointpattern/
├── src/stipple/          # the importable library (pip install -e .)
│   ├── __init__.py       # public API
│   └── core.py           # K/L, envelope, SES, windows, edge correction, preprocessing
├── notebooks/            # 00–17: the analysis, in order (00 is now superseded by src/stipple)
│   └── archive/          # superseded notebooks
├── scripts/              # HPC screen: batch runner, SLURM arrays, aggregation
│   ├── batch_k_analysis.py   run_array.sh   aggregate.py   (+ v2 variants)
│   └── writeup_figs/     # dissertation figure builders
├── results/              # HPC job manifests (screen outputs are gitignored)
├── docs/                 # poster and talk assets
├── pyproject.toml        # packaging + core dependencies
├── requirements.txt      # full pinned lock (exact reproduction environment)
└── README.md
```

## Key functions (`stipple`)

| Function | Description |
|---|---|
| `run_pair_analysis(strip_df, gene_a, gene_b, r_vals, window_type='hull', n_sim=99)` | Full per-pair pipeline; single entry point |
| `bivariate_k(coords_a, coords_b, r_vals, window, resolution=64)` | Cross-type Ripley's K, edge-corrected |
| `k_to_l(k_vals, r_vals)` | Variance-stabilising transform L(r) = sqrt(K/pi) - r |
| `compute_envelope(coords_a, coords_b, r_vals, window, n_sim=99)` | Label-permutation envelope |
| `get_window(df, window_type='hull')` | Observation window (rect tuple, or Shapely polygon) |
| `fraction_inside_hull(point, r, hull, resolution=64)` | Edge-correction weight for a radius-r disc |
| `get_coords(df, gene)` | Extract an (N, 2) coordinate array for one gene |

**Edge-correction note.** `fraction_inside_hull` introduces a small systematic negative bias in absolute K(r). Because that bias is identical across observed and permuted patterns, it cancels in the permutation test, so absolute L(r) is never interpreted, only envelope-relative significance.

## Controls

| Control | Gene pair | Expectation | Peak SES (strips 1/2/3) |
|---|---|---|---|
| Positive | KRT8 × KRT18 | co-expressed keratins, near midline | +0.20 / +0.23 / +0.16 |
| Negative | KRT8 × SCGB3A1 | distinct compartments, anti-associated | −3.63 / −2.41 / −1.59 |
| Negative | MALAT1 × KRT18 | ubiquitous × epithelial, anti-associated | −2.03 / −2.69 / −1.54 |

Under the production concave-hull window the binary envelope test is under-powered for all pairs (the tight window concentrates the label-swap null), so controls are discriminated by the continuous SES rather than by envelope crossings. This window-as-lever finding is the study's main methodological result.

## Results summary

The screen ran 614 (pair, strip) tests (605 unique) with no failed jobs, at `n_sim=199` (v1 panel) and `n_sim=99` (v2), `seed=42` throughout. It calibrates correctly on the controls and orders the pre-registered biological groups as designed. Strict detections are sparse: **1 pair, SPP1 × CD44 in the infected strip, is the single strict envelope hit** and the only surviving network edge. The infected strip is not globally elevated; its co-localisation structure differs from the controls at specific, biologically coherent pairs (immune, co-stimulatory, interferon and growth-factor signalling). The results are a proof of concept on a single slide; the associations are candidate co-localisations, not confirmed communication.

## Acknowledgements

CosMx data provided by the Immunomodulation lab at Guy's and St Thomas' Hospital. Supervision and guidance from Dr Dan Nicolau. Ligand-receptor annotations from CellChatDB (Jin et al., *Nature Communications*, 2021). HPC compute on the KCL CREATE cluster.

## License

Released under the MIT License. See [`LICENSE`](LICENSE).
