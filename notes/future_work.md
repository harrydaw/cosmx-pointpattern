# Future work

Items that fell outside the scope of the 16 July 2026 submission but are worth
recording so the obvious next moves aren't lost. Each entry lists *why* it
matters, a sketched *approach*, rough *effort*, and a pointer to the *current
code* that would change. Items 1 and 2 also have a detailed implementation
appendix at the end of this document.

---

## 1. Vectorise the Shapely hull edge correction in `bivariate_k`

**Why.** The per-pair Shapely call `Point(p).buffer(r).intersection(hull)` is
pure-Python and dominates `bivariate_k` cost for any gene with more than a few
thousand transcripts. The MALAT1 × KRT18 diagnostic in nb13 made this concrete:
a single `bivariate_k` call at native N (16,808) takes ~7 min, and 99
permutations would take ~11.5 h. The tactical fix landed in this session was a
`max_n` subsample cap (defaults `None`, used in nb13 §7 at `max_n=3000`). That
mitigation works but costs a footnote in the methods section — "we cap abundant
genes at N=3000". A vectorised edge correction removes the need for the cap
entirely and lets future selected-pair HPC batches run on full N.

**Approach (sketch).** Rasterise the concave hull once into a high-resolution
binary mask. Precompute each point's distance-to-boundary. The disc-fraction-
inside-window then becomes a function of `(boundary_distance, r)` — implementable
as a 2D lookup table, or as a closed-form approximation for points further than
`r_max` from the boundary (most of them, on tissue interiors). Validate the new
path against the current Shapely implementation on KRT8 × KRT18 strip_2 at
production resolution, and against the rectangular arc-fraction path for
consistency. Hold the API stable: `bivariate_k`'s signature can keep `max_n` as
a vestigial knob (default `None`) for graceful backward compatibility with any
saved diagnostic results.

**Effort.** ~1–2 days including a side-by-side validation notebook.

**Current code.** `notebooks/00_functions.ipynb` → `bivariate_k`, the
`else:` (Shapely Polygon) branch that calls `fraction_inside_hull`.

**Detailed proposal:** see Appendix B.

---

## 2. Leiden over Louvain for community detection

**Why.** Louvain has documented failure modes — badly-connected communities and
occasional disconnected partitions. Leiden (Traag, Waltman & van Eck, 2019)
fixes these with a guaranteed-connected refinement step. If any downstream
biological claim leans on community structure (which it does, via the nb12
network), the modularity-optimiser choice should be the defensible one.

**Approach.** Drop-in replacement using `leidenalg` (with `python-igraph`). Same
resolution-parameter semantics as Louvain. On the existing nb12 graph, compare
partition modularity, number of disconnected components per community, and
whether the dominant communities are stable across optimisers. If the numbers
move materially, the network figure needs regenerating.

**Effort.** ~half day including a comparison cell.

**Current code.** `notebooks/12_network_analysis.ipynb` — find the Louvain call
site.

**Detailed proposal:** see Appendix A.

---

## 3. Pseudo-segmentation via marker-gene clustering

**Why.** CosMx is segmentation-free, which the current pipeline treats as a
virtue (no segmentation bias). The cost of that choice is that L(r) is computed
on all transcripts globally, with no concept of cell-type — so the test answers
"do these two genes co-cluster across the whole tissue?" rather than the much
stronger "do these two genes co-cluster within macrophages?". A
pseudo-segmentation step would unlock per-cell-type co-localisation analysis
without committing to a hard segmentation primitive.

**Approach (sketch).** Spatially cluster transcripts at a cell-scale bandwidth
(DBSCAN at ~5–15 px is a reasonable starting point, reusing the primitive
already used for noise removal in nb08). For each cluster, build a marker-gene
composition vector. Classify the cluster against a known marker panel
(epithelial / immune / fibroblast / endothelial) — supervised if a panel is
available, or unsupervised via Gaussian mixture in marker-composition space
otherwise. Every transcript inherits the label of its containing cluster. Then
`bivariate_k` can be restricted to within-label pairs, giving per-type L(r)
curves.

This is its own research subproject and is unlikely to fit before the
submission deadline. Worth scoping as a follow-up paper rather than a thesis
section.

**Effort.** ~1–2 weeks for a defensible implementation.

**Current code.** DBSCAN is already wired up in `notebooks/08_improved_QC.ipynb`
— that's the entry point for the spatial clustering primitive.

---

## Suggested sequencing (items 1 & 2)

**Do Leiden first** (smaller, isolated, immediately useful for nb12). **Do
raster second** (longer, needs its own validation notebook, but the upside is
the extended-panel HPC batch).

If only one fits before the writeup begins: Leiden. The network finding is the
primary deliverable; a defensible community-detection choice matters more than
the per-call runtime of the K function (which `max_n=3000` already handles for
any practical pair).

---

## Appendix A — Detailed proposal: Leiden as `community_method='leiden'`

### API design

Move the community-detection call into a small helper in
`notebooks/00_functions.ipynb`, then make it accept a `method` kwarg. Proposed
signature:

```python
def detect_communities(
    G: nx.Graph,
    method: str = 'louvain',         # 'louvain' | 'leiden' | 'greedy'
    weight: str = 'weight',
    seed: int = 42,
    resolution: float = 1.0,         # gamma for both Louvain and Leiden
) -> list[set]:
    """
    Returns a list of node sets, sorted largest-first.
    For 'louvain'/'greedy' uses NetworkX directly. For 'leiden' uses
    leidenalg + python-igraph with ModularityVertexPartition (same Q
    objective as Louvain, so resolution semantics line up).
    """
```

Resolution semantics: `leidenalg.find_partition` with
`ModularityVertexPartition` accepts `resolution_parameter=gamma` — same gamma
as NetworkX `louvain_communities` `resolution=`. Document this in the docstring
so future-Harry doesn't need to re-derive it.

### Implementation steps

1. Add `leidenalg` and `python-igraph` to `requirements.txt`. Check the HPC venv
   can install them on `python/3.11.6-gcc-13.2.0`.
2. Write `detect_communities` in `00_functions.ipynb` next to the other network
   helpers (if any exist yet) or in a new "Network" section. Internal dispatch:
   - `'louvain'` → existing `louvain_communities(G, weight=weight, seed=seed, resolution=resolution)`
   - `'leiden'` → convert `G` to an `igraph.Graph` (preserving edge weights),
     call `leidenalg.find_partition(g, leidenalg.ModularityVertexPartition,
     weights='weight', resolution_parameter=resolution, seed=seed)`, then map
     igraph vertex IDs back to networkx node names. Return list of sets sorted
     by length descending.
   - `'greedy'` → existing `greedy_modularity_communities(G, weight=weight)`
     (no seed needed; deterministic).
3. Replace the direct `louvain_communities(...)` calls in nb12 with
   `detect_communities(G, method='louvain', ...)` so the swap is one line.

### Validation

Add a side-by-side cell to nb12 after the community-assignment step:

| Strip | Method | Q | n_communities | n_disconnected_components |
|---|---|---|---|---|
| strip_1 | louvain | … | … | … |
| strip_1 | leiden  | … | … | … |
| … | … | … | … | … |

`n_disconnected_components` is the count, summed over communities, of extra
connected components within each community (0 = all communities are internally
connected). Leiden should give 0 by construction; Louvain will sometimes give
>0. This is the headline argument for Leiden being more defensible.

### Methodological flags

- Both Leiden and Louvain are stochastic — report at least one seed, ideally a
  stability check across 5 seeds (already on the wishlist per
  `notes/notebook_audit.md` "Outstanding Issues" item 3).
- Leiden has `iterations` and `n_iterations` knobs in leidenalg. Default is
  fine, but document the choice.
- If a reviewer asks "why Leiden": it guarantees internally connected
  communities and avoids Louvain's badly-connected-community failure mode.
  Cite Traag et al. 2019.

---

## Appendix B — Detailed proposal: rasterised edge correction as `edge_correction='raster'`

### API design

Add an `edge_correction` parameter to `bivariate_k` (and pass-through from
`compute_envelope` and `run_pair_analysis`). Proposed signature:

```python
def bivariate_k(
    coords_a, coords_b, r_vals, window,
    resolution: int = 64,
    max_n: int | None = None,
    seed: int | None = None,
    edge_correction: str = 'shapely',     # 'shapely' | 'raster'
    raster_resolution: int = 1024,        # pixels along the longer axis
) -> np.ndarray:
```

- `edge_correction='shapely'` is the existing path (default).
- `edge_correction='raster'` triggers the new path. Only meaningful for
  Shapely-Polygon windows; for tuple (rectangle) windows the existing
  arc-fraction analytical correction is already fast — fall through to it
  regardless of the kwarg.
- `raster_resolution` controls the mask resolution. 1024 along the longer axis
  is the proposed default — tunable.

Helper to add to `00_functions.ipynb`:

```python
class HullRasterEdgeCorrection:
    def __init__(self, hull, raster_resolution=1024):
        ...  # build binary mask + boundary-distance map

    def fraction(self, point, r):
        ...  # return disc-fraction-inside-hull for a single (point, r)

    def fraction_batch(self, points, r):
        ...  # vectorised version for many points at one r
```

The class form makes it trivial for `bivariate_k` to build the corrector once
and reuse it across the whole simulation loop (saves ~100× the build cost on a
99-permutation envelope).

### Implementation steps

1. **Mask + distance transform.** Rasterise the hull at `raster_resolution`
   pixels along the longer axis. Use `rasterio.features.rasterize` or
   `shapely`+`numpy` (Shapely 2.0 has a `shapely.contains_xy` vectorised
   predicate that takes an array of points → boolean mask in one call). Then
   compute the distance transform of the *outside* of the mask using
   `scipy.ndimage.distance_transform_edt`. This gives, for every pixel inside
   the hull, its Euclidean distance to the nearest boundary pixel.
2. **Per-point boundary distance.** For each query point, map (x, y) → (row,
   col) → lookup the distance from the transform. Vectorisable.
3. **Disc-fraction lookup.** For each (boundary_distance, r):
   - if boundary_distance ≥ r: weight = 1.0 (disc entirely inside hull —
     trivial early exit, applies to most points on the tissue interior).
   - else: integrate the in-mask area of a disc of radius r centred at the
     point, using the mask itself. For a single point this is one 2D numpy
     slice + sum — fast.
4. **Wire into `bivariate_k`.** In the `else:` (polygon) branch, build the
   corrector once before the per-point loop, then call
   `corrector.fraction((x_i, y_i), d_ij)` in place of
   `fraction_inside_hull(...)`. Leave the `frac = max(frac, 0.01)` guard and
   the rest of the loop unchanged.

### Validation

The raster path is an approximation and must be shown to agree with the Shapely
ground truth before any production use. Create
`notebooks/14_edge_correction_validation.ipynb` (or similar):

1. Load the standard test case: KRT8 × KRT18, strip_2, concave hull
   ratio=0.1, production r_vals.
2. Compute `l_obs` and a small envelope (n_sim=9) twice — once with
   `edge_correction='shapely'`, once with `'raster'` at
   `raster_resolution=1024`.
3. Plot the two `l_obs` curves overlaid. They should be visually
   indistinguishable at production radii.
4. Quantify: max absolute and relative deviation in `l_obs` and in envelope
   width across r_vals. Set a tolerance (e.g. `<1%` relative deviation at all
   r ≤ r_max) and assert it.
5. Repeat at `raster_resolution=512` and `2048` to show the result is stable
   across mask resolutions — this is the methods-section defence.
6. Time both paths on MALAT1 × KRT18 strip_2 (the original pathological case).
   Report the speedup.

Only after this notebook is green should the raster path become the default —
and even then I'd argue for keeping `'shapely'` as the documented exact
reference and `'raster'` as a fast approximation, so any single result can be
re-run for verification.

### Methodological flags

- The disc-fraction-via-mask is an **approximation**, with error bounded by the
  pixel size. At 1024-px raster resolution on a ~5000 px tissue strip, pixel
  size ≈ 5 px. For radii ≥ 25 px this is a sub-% area error; for the smallest
  production radius (5 px) it's larger. The validation notebook needs to
  quantify the small-r regime explicitly.
- For points within `r` of the boundary but not in any pathological pocket, the
  mask integration is faithful. For points in concave pockets of the hull, the
  mask captures the geometry better than the current Shapely call would
  (Shapely uses an exact intersection but with a 64-sided polygonal disc
  approximation — also an approximation, just a different one).
- `max_n` stays in the function signature as a documented fallback. Phrasing
  for the methods section: "We default to a rasterised edge correction
  (validated against the exact Shapely intersection in nb14); for production
  results requiring the exact correction at large N, the `max_n` parameter
  caps per-call cost at the price of an explicit subsample."

### Open questions to resolve before implementation

1. **Does the HPC venv allow `leidenalg`?** Need to install it on the KCL
   CREATE venv (`/scratch/users/k25115761/NoSegs/.../.venv`) and confirm
   imports work in a SLURM job. `python-igraph` ships as a wheel for Linux
   x86_64 so should be fine.
2. **What's the right default `raster_resolution`?** Validation notebook will
   answer this — start at 1024 and shift up or down based on the small-r
   error.
3. **Should `run_pair_analysis` expose `edge_correction` and
   `raster_resolution`?** Yes — same pattern as `max_n` would have been
   threaded through. Plan to add them at the same time as the raster path is
   wired into `bivariate_k`.
4. **Reproducibility note.** The raster path is deterministic given the mask
   resolution and the hull. No new seed needed. Document this so the docstring
   doesn't accidentally suggest otherwise.
