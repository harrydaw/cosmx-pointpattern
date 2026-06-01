# Future implementation plan — Leiden + rasterised edge correction

**Status:** Drafted before nb12 (network analysis) is finalised. Revisit
after the network notebook lands, since both items touch code that nb12
will exercise heavily and we want the API to fit what nb12 actually needs.

**Scope.** Two short-term improvements, each gated on a function-level
parameter so the existing behaviour stays the default:

1. **Leiden** as an alternative to Louvain for modularity-based community detection.
2. **Rasterised edge correction** in `bivariate_k` as an alternative to the per-pair Shapely intersection.

In both cases the existing behaviour is the default, the new behaviour is
opt-in via a kwarg, and `max_n` (the subsample fallback added in the
MALAT1 session) is retained — it remains the methodologically explicit
safety net for cases where the rasterised approximation is judged
unacceptable, or where a reviewer asks for a no-approximation result.

---

## Item 1 — Leiden as a `community_method='leiden'` option

### Why
Louvain has documented failure modes — badly-connected communities, and
occasional disconnected partitions (a single community spanning two
graph components linked only via the meta-graph). Leiden (Traag, Waltman
& van Eck, 2019) fixes both with a guaranteed-connected refinement step.
For the network notebook, having Leiden available as a one-flag swap
makes the community detection result more defensible without committing
to abandoning Louvain — both can be reported side-by-side.

### API design
Move the community-detection call into a small helper in
`notebooks/00_functions.ipynb`, then make it accept a `method` kwarg.
Proposed signature:

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
`ModularityVertexPartition` accepts `resolution_parameter=gamma` — same
gamma as NetworkX `louvain_communities` `resolution=`. Document this in
the docstring so future-Harry doesn't need to re-derive it.

### Implementation steps
1. Add `leidenalg` and `python-igraph` to `requirements.txt` /
   `pyproject.toml` (whichever the repo uses). Check the HPC venv can
   install them on `python/3.11.6-gcc-13.2.0`.
2. Write `detect_communities` in `00_functions.ipynb` next to the other
   network helpers (if any exist yet) or in a new "Network" section.
   Internal dispatch:
   - `'louvain'` → existing `louvain_communities(G, weight=weight, seed=seed, resolution=resolution)`
   - `'leiden'` → convert `G` to an `igraph.Graph` (preserving edge weights), call `leidenalg.find_partition(g, leidenalg.ModularityVertexPartition, weights='weight', resolution_parameter=resolution, seed=seed)`, then map igraph vertex IDs back to networkx node names. Return list of sets sorted by length descending.
   - `'greedy'` → existing `greedy_modularity_communities(G, weight=weight)` (no seed needed; deterministic).
3. Replace the direct `louvain_communities(...)` calls in nb12 with
   `detect_communities(G, method='louvain', ...)` so the swap is one line.

### Validation
Add a side-by-side cell to nb12 after the community-assignment step:

| Strip | Method | Q | n_communities | n_disconnected_components |
|---|---|---|---|---|
| strip_1 | louvain | … | … | … |
| strip_1 | leiden  | … | … | … |
| … | … | … | … | … |

`n_disconnected_components` is the count, summed over communities, of
extra connected components within each community (0 = all communities
are internally connected). Leiden should give 0 by construction;
Louvain will sometimes give >0. This is the headline argument for
Leiden being more defensible.

### Effort
Half a day, including the comparison cell. Most of the time is on the
igraph ↔ networkx node-name round-trip and on a sanity-check that
`resolution_parameter` semantics actually match (verify by checking
that with `resolution=1.0`, two algorithms with the same partition
yield the same Q).

### Methodological flags
- Both Leiden and Louvain are stochastic — report at least one seed,
  ideally a stability check across 5 seeds (already on the wishlist
  per `notes/notebook_audit.md` "Outstanding Issues" item 3).
- Leiden has `iterations` and `n_iterations` knobs in leidenalg. Default
  is fine, but document the choice.
- If a reviewer asks "why Leiden": it guarantees internally connected
  communities and avoids Louvain's badly-connected-community failure
  mode. Cite Traag et al. 2019.

---

## Item 2 — Rasterised edge correction as `edge_correction='raster'`

### Why
The per-pair Shapely call `Point(p).buffer(d_ij).intersection(hull)` in
`fraction_inside_hull` is the cost driver for `bivariate_k` on any
gene with >~3k transcripts. Rasterising the hull once and doing
per-point boundary-distance lookups eliminates the Shapely call from
the inner loop. This is what enables the extended-panel HPC run to use
full N — no more `max_n=3000` footnote in the methods section.

### API design
Add an `edge_correction` parameter to `bivariate_k` (and pass-through
from `compute_envelope` and `run_pair_analysis`). Proposed signature:

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
  arc-fraction analytical correction is already fast — fall through to
  it regardless of the kwarg.
- `raster_resolution` controls the mask resolution. 1024 along the
  longer axis is the proposed default — tunable.

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

The class form makes it trivial for `bivariate_k` to build the
corrector once and reuse it across the whole simulation loop (saves
~100× the build cost on a 99-permutation envelope).

### Implementation steps
1. **Mask + distance transform.** Rasterise the hull at
   `raster_resolution` pixels along the longer axis. Use
   `rasterio.features.rasterize` or `shapely`+`numpy` (Shapely 2.0 has
   a `shapely.contains_xy` vectorised predicate that takes an array of
   points → boolean mask in one call). Then compute the distance
   transform of the *outside* of the mask using
   `scipy.ndimage.distance_transform_edt`. This gives, for every pixel
   inside the hull, its Euclidean distance to the nearest boundary
   pixel.
2. **Per-point boundary distance.** For each query point, map (x, y) →
   (row, col) → lookup the distance from the transform. Vectorisable.
3. **Disc-fraction lookup.** For each (boundary_distance, r):
   - if boundary_distance ≥ r: weight = 1.0 (disc entirely inside hull —
     trivial early exit, applies to most points on the tissue interior).
   - else: integrate the in-mask area of a disc of radius r centred at
     the point, using the mask itself. For a single point this is one
     2D numpy slice + sum — fast.
4. **Wire into `bivariate_k`.** In the `else:` (polygon) branch, build
   the corrector once before the per-point loop, then call
   `corrector.fraction((x_i, y_i), d_ij)` in place of
   `fraction_inside_hull(...)`. Leave the `frac = max(frac, 0.01)`
   guard and the rest of the loop unchanged.

### Validation
This is the non-negotiable part — the raster path is an approximation
and must be shown to agree with the Shapely ground truth before any
production use.

Create `notebooks/14_edge_correction_validation.ipynb` (or similar):
1. Load the standard test case: KRT8 × KRT18, strip_2, concave hull
   ratio=0.1, production r_vals.
2. Compute `l_obs` and a small envelope (n_sim=9) twice — once with
   `edge_correction='shapely'`, once with `'raster'` at
   `raster_resolution=1024`.
3. Plot the two `l_obs` curves overlaid. They should be visually
   indistinguishable at production radii.
4. Quantify: max absolute and relative deviation in `l_obs` and in
   envelope width across r_vals. Set a tolerance (e.g. `<1%` relative
   deviation at all r ≤ r_max) and assert it.
5. Repeat at `raster_resolution=512` and `2048` to show the result is
   stable across mask resolutions — this is the methods-section
   defence.
6. Time both paths on MALAT1 × KRT18 strip_2 (the original
   pathological case). Report the speedup.

Only after this notebook is green should the raster path become the
default — and even then I'd argue for keeping `'shapely'` as the
documented exact reference and `'raster'` as a fast approximation, so
any single result can be re-run for verification.

### Effort
1–2 days, broken down:
- Mask + distance transform helper: half a day.
- Wiring into `bivariate_k` and `compute_envelope`: a few hours.
- Validation notebook: half to one day.

### Methodological flags
- The disc-fraction-via-mask is an **approximation**, with error
  bounded by the pixel size. At 1024-px raster resolution on a ~5000 px
  tissue strip, pixel size ≈ 5 px. For radii ≥ 25 px this is a sub-%
  area error; for the smallest production radius (5 px) it's larger.
  The validation notebook needs to quantify the small-r regime
  explicitly.
- For points within `r` of the boundary but not in any pathological
  pocket, the mask integration is faithful. For points in concave
  pockets of the hull, the mask captures the geometry better than the
  current Shapely call would (Shapely uses an exact intersection but
  with a 64-sided polygonal disc approximation — also an
  approximation, just a different one).
- `max_n` stays in the function signature as a documented fallback.
  Phrasing for the methods section: "We default to a rasterised edge
  correction (validated against the exact Shapely intersection in
  nb14); for production results requiring the exact correction at
  large N, the `max_n` parameter caps per-call cost at the price of an
  explicit subsample."

---

## Suggested sequencing

**Do Leiden first** (smaller, isolated, immediately useful for nb12).
**Do raster second** (longer, needs its own validation notebook, but
the upside is the extended-panel HPC batch).

If only one fits before the writeup begins: Leiden. The network
finding is the primary deliverable; a defensible community-detection
choice matters more than the per-call runtime of the K function (which
`max_n=3000` already handles for any practical pair).

---

## Open questions to resolve before implementation

1. **Does the HPC venv allow `leidenalg`?** Need to install it on the
   KCL CREATE venv (`/scratch/users/k25115761/NoSeggs/.../.venv`) and
   confirm imports work in a SLURM job. `python-igraph` ships as a
   wheel for Linux x86_64 so should be fine.
2. **What's the right default `raster_resolution`?** Validation
   notebook will answer this — start at 1024 and shift up or down
   based on the small-r error.
3. **Should `run_pair_analysis` expose `edge_correction` and
   `raster_resolution`?** Yes — same pattern as `max_n` would have
   been threaded through. Plan to add them at the same time as the
   raster path is wired into `bivariate_k`.
4. **Reproducibility note.** The raster path is deterministic given
   the mask resolution and the hull. No new seed needed. Document this
   so the docstring doesn't accidentally suggest otherwise.
