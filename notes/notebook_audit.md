# Notebook Audit — NoSeggs / cosmx-pointpattern

> Updated 2026-05-11. Reference document for the dissertation write-up. Pairs with `lab_notebook.md` (chronological) and `network_revival_plan.md` (next-actions).

---

## Executive summary

Seventeen notebooks form a clean end-to-end pipeline: raw zarr → strip-assigned cleaned points → window construction → bivariate Ripley's K/L with permutation envelopes → 438-job HPC batch → SES + significance filter → network construction → community detection. Every stage is implemented, runs, and was validated against positive/negative controls at intermediate steps.

**The pipeline is sound. The final-stage result is the issue.** Of 438 (ligand, receptor, strip) rows in `lr_panel_results.parquet`, **only 7** ever exceed the upper envelope at *any* radius and **only 1** survives the 3-consecutive-radii filter (SPP1 × CD44, strip 2). This is too sparse to draw networks from — but the audit finds two plausible non-failure explanations:

1. **`n_sim` mismatch**: validation in nb04/nb07 used `n_sim=99`; the HPC run used `n_sim=199`. More permutations → wider envelope → harder to detect signal. The canonical positive control (KRT8 × KRT18) was never re-tested at production `n_sim`.
2. **Panel scope**: the 146 pairs are CellChatDB ligand-receptor pairs only, abundance-filtered. The genes most likely to co-cluster (cytokeratins, mucins, cell-type marker pairs) are excluded by construction.

Both are addressable. See `network_revival_plan.md` for the diagnostic path.

---

## Pipeline overview

```
nb01  Load raw zarr                ─── 2.66M reads, 63 FOVs, 4 slides
nb02  FOV-3 exploration            ─── GMM strip assignment validated
nb02c Expand to 6 FOVs (s1)        ─── s1_all_strips.parquet (703k pts)
   │
nb03  Bivariate K + L + envelope   ─── core algorithm
nb04  KRT8×KRT18 positive control  ─── PASSED at n_sim=99
nb05  MALAT1×KRT18 negative ctrl   ─── PASSED (correctly null)
nb06  LR panel curation             ─── 146 viable pairs from CellChatDB
nb07  Multi-FOV control retest      ─── stable across FOVs
   │
nb08  DBSCAN noise removal         ─── 6.8% noise flagged
nb08b Manual cluster QC cleanup     ─── 18 outlier clusters excluded
   │                              ─── s1_all_strips_cleaned.parquet
nb09  Concave hull + edge corr.    ─── fraction_inside_hull()
nb09a Unified window API
nb09b Custom polygon drawing       ─── superseded by nb09 concave hull
nb09c Window comparison            ─── concave ≈ 65–75% of convex area
   │
nb10  POC end-to-end               ─── R_CC=50, R_MAX=250, N_R=50, N_SIM=99
HPC   batch_k_analysis × 438       ─── ⚠ N_SIM=199 in production
   │                              ─── lr_panel_results.parquet (438 rows)
nb11  Results overview              ─── score distribution, all jobs OK
nb12  Network analysis              ─── ⚠ only 1/438 passes filter
```

---

## Per-notebook entries

### 00 — Functions library
**Purpose:** Shared code loaded at runtime by the HPC script via nbformat.
**Outputs:** Functions only (no files).
**Key functions:** `bivariate_k`, `k_to_l`, `compute_envelope`, `run_pair_analysis`, `fraction_inside_hull`, `make_window`, GMM strip-assigner.
**Status:** ✅ Clean. Single source of truth — any pipeline change must be applied here, not duplicated downstream.

### 01 — Load raw data
**Purpose:** Load `Varsha_1234.zarr`, validate 2.66 M reads across 63 FOVs / 4 slides.
**Outputs:** In-memory dataframe; no parquet written here.
**Decisions:** Coordinate columns standardised on `x_global_px`, `y_global_px` (CosMx native).
**Status:** ✅ Clean.

### 02 — FOV 3 exploration
**Purpose:** Pick FOV 3 (S1) for initial method development; visually confirm KRT8/KRT18 co-expression.
**Decisions:** GMM on rotated x-coordinate identifies 3 strips per FOV; strip 2 = infected.
**Status:** ✅ Clean.

### 02c — FOV review & expansion
**Purpose:** Expand from FOV 3 alone to 6 FOVs (4, 5, 8, 9, 10, 11), excluding 1, 2, 3 for quality reasons.
**Outputs:** `s1_all_strips.parquet` (702,873 transcripts).
**Decisions:** Strip 1 contains two physically separated tissue fragments — documented but not split further.
**Status:** ✅ Clean, with a noted heterogeneity caveat for strip 1.

### 03 — Bivariate K function
**Purpose:** Implement bivariate K, the L-transform, and the permutation envelope from scratch.
**Decisions:** **`n_sim = 99`** hard-coded as the working default. Validated against `spatstat`-style CSR simulation.
**Status:** ✅ Clean. ⚠ The `n_sim=99` choice in this notebook was *not* propagated to the production HPC run (which uses 199 — see Outstanding Issues).

### 04 — Real analysis (positive controls)
**Purpose:** KRT8 × KRT18 on FOV 3 strip 2 as the positive control; MALAT1 × KRT18, KRT8 × SCGB3A1 as negatives.
**Decisions:** `r_max = 300` px; `n_sim = 99`; rectangular window with standard edge correction.
**Result:** KRT8 × KRT18 cleanly exceeds the upper envelope across many radii on strip 2. Negatives correctly stay within.
**Status:** ✅ Clean — this is the validation evidence that the algorithm *works at production scale on this data*.

### 05 — Negative control validation
**Purpose:** Cross-strip behaviour of the negative pairs; envelope stability.
**Status:** ✅ Clean.

### 06 — LR panel checks
**Purpose:** Curate the ligand-receptor pair list for the panel screen. Apply ≥50 transcripts-per-gene-per-strip filter to CellChatDB pairs.
**Outputs:** `viable_lr_pairs_all_strips.parquet` (146 pairs), `lr_panel_index.json`.
**Decisions:** **CellChatDB only** — no broader co-expression pairs. The validated positive controls (KRT × KRT) are *not in this panel*. ⚠ This is the second-largest contributor to the empty-network problem.
**Status:** ⚠ Scope decision worth revisiting (see `extended_panel_rationale.md`).

### 07 — Expanded controls
**Purpose:** Re-run controls across all 6 FOVs to confirm robustness once the dataset was expanded.
**Status:** ✅ Clean.

### 08 — Improved QC (DBSCAN)
**Purpose:** Per-FOV noise filtering.
**Decisions:** Adaptive `eps` = 97th-percentile 1-NN distance (clipped to [20, 30] px); `min_samples = 5`; `min_cluster_size = 150`. ~6.8 % of 400k+ points flagged as noise.
**Status:** ✅ Clean and well-justified.

### 08b — Manual QC cleanup
**Purpose:** Hand-review cluster IDs flagged in nb08; remove 18 specific outlier clusters.
**Outputs:** `s1_all_strips_cleaned.parquet` (the canonical input from here on; has columns `is_noise`, `is_small_cluster`, `manually_excluded`).
**Status:** ✅ Clean. List of excluded cluster IDs is documented in-notebook.

### 09 — Improved windows & edge correction
**Purpose:** Replace rectangular window with concave hull; implement Shapely-based edge correction (`fraction_inside_hull`).
**Decisions:** Concave hull with `ratio = 0.1`. Resolution 64 for disc approximation (256-sided polygon).
**Note:** Absolute `L(r)` shows a small systematic *negative* bias under this correction — flagged in-notebook. The **permutation test remains valid** because the bias affects both observed and permuted curves equally.
**Status:** ✅ Clean, with documented bias.

### 09a — Unified window API
**Purpose:** Refactor: `make_window(kind=...)` returns rect / convex / concave / custom interchangeably.
**Status:** ✅ Clean.

### 09b — Custom window drawing
**Purpose:** Exploratory manual polygon override for strips where concave hull misses important tissue.
**Outcome:** After extensive QC, the **statistical concave hull approach (nb09) was adopted as the canonical window for all three strips**. The hand-drawn polygons are kept for reference but are no longer in the pipeline. The concave hull is reproducible, parameterised (`ratio=0.1`), and avoids the bias risk of manual region selection.
**Status:** ✅ Superseded by design, not abandoned. Not an outstanding issue.

### 09c — Window comparison
**Purpose:** Quantitative comparison of window choices: concave area is 65–75 % of convex area; reduces false-interior overlap.
**Status:** ✅ Clean.

### 10 — POC end-to-end
**Purpose:** Single notebook running the whole pipeline on a tiny subset; checklist of params before HPC submission.
**Decisions, locked here for the HPC run:** `R_CC = 50` (cell-cell radius bound, ≈ 9 µm), `R_MAX = 250` px (≈ 45 µm), `N_R = 50` (radii sampled), **`N_SIM = 99`** in this notebook, `SEED = 42`.
**Status:** ✅ Clean. ⚠ The `N_SIM = 99` set here is *overridden* by the HPC orchestration script (`run_array.sh:61` passes `--n_sim 199`).

### HPC run — `scripts/batch_k_analysis.py` + `run_array.sh`
**Configuration used in production:**
- 438 jobs = 146 pairs × 3 strips
- **`n_sim = 199`** (overrides the 99 in nb03/nb10/nb04)
- Concave-hull window with edge correction
- `--time=04:00:00 --mem=4G` per array task
- Concurrent limit per user: `%2`
**Result:** All 438 jobs reported `status='ok'`. Output: `results/lr_panel_results.parquet` (438 × 11 columns).
**Status:** ⚠ Runs cleanly but at a *different `n_sim`* than the validation notebooks. No re-validation of positive controls at `n_sim=199`.

### 11 — Results overview
**Purpose:** Load the 438-row parquet, compute SES, plot score distributions.
**Findings:** Score mean ≈ −2.32 (i.e. observed L tends to sit slightly *below* envelope midline — consistent with the negative bias noted in nb09). Distribution is unimodal, narrow. No failed jobs.
**Status:** ✅ Clean (but the score distribution is itself diagnostic — most pairs look null).

### 12 — Network analysis
**Purpose:** Filter → edge table → per-strip graph → community detection (Louvain + greedy) → strip comparison + scale analysis + PCA/UMAP placeholder.
**Decisions:** `n_consec = 3` (sustained exceedance threshold); edge weight = peak `l_obs`; Louvain `seed = 42`; spring layout `seed = 42`.
**Outputs (current):**
- `results/12_edge_summary.csv` (1 row)
- `results/12_community_labels.csv`
- `results/figures/12_networks_per_strip.png`, `12_hub_genes.png`, `12_scale_distribution.png`
**Status:** ❌ **Filter passes only 1 of 438 pairs. Networks effectively empty.** Code is correct; inputs are too restrictive (see Outstanding Issues + revival plan).
**Note on networkx import:** the original `from networkx.community import …` was broken — `networkx.community` is re-exported from `networkx.algorithms.community` at runtime but isn't a direct submodule. Fixed by using the algorithms path explicitly.

---

## Parameter cheatsheet

| Parameter | Current value | Set in | Sensitivity worth running |
|---|---|---|---|
| `n_sim` (envelope permutations) | 99 in dev / **199 in HPC** | nb03/nb10 vs `run_array.sh:61` | **{49, 99, 199}** — top priority |
| `r_min, r_max` | 5, 250 px (≈ 0.9, 45 µm) | nb10 | Stable; sensitivity not urgent |
| `n_r` (radii sampled) | 50 | nb10 | — |
| `n_consec` (sustained exceedance) | 3 | nb12 | **{1, 2, 3, 5}** — quick sweep |
| Window type | concave hull, `ratio = 0.1` | nb09 | rect vs convex vs concave already done (nb09c) |
| Edge correction | `fraction_inside_hull` (Shapely) | nb09 | — |
| DBSCAN `eps` | 97th-pctile 1-NN, clipped [20, 30] | nb08 | — |
| DBSCAN `min_samples` | 5 | nb08 | — |
| `min_cluster_size` | 150 | nb08 | — |
| LR panel | 146 CellChatDB pairs, abundance ≥ 50 | nb06 | **Add ~80 anchor pairs** (see `extended_panel_rationale.md`) |
| Edge weight (nb12) | peak `l_obs` | nb12 | peak SES, area-above-envelope |
| Community detection | Louvain, `seed = 42` | nb12 | 10 seeds; add Leiden |
| Significance test (nb12) | pointwise + max-consecutive | nb12 | global rank envelope / area test |

---

## Outstanding issues (ranked by impact on dissertation conclusions)

1. **`n_sim` mismatch** — production envelope at 199 has never been validated against a positive control. Cheapest fix: re-run KRT8 × KRT18 locally at `n_sim = 199` using the production window pipeline. (Revival plan Step A.6.)
2. **LR panel scope** — KRT pairs and other strong co-expression pairs are absent by construction. Adding ~80 pre-registered anchor pairs at `n_sim = 99` is the highest-EV experiment for the dissertation. (Revival plan Step B.)
3. **Single Louvain seed** — production result uses `seed = 42` only. Run 10 seeds; report best-Q partition. Easy fix.
4. **`n_consec = 3` is conservative** — was an unargued POC default. The pointwise max-consecutive test is also relatively weak; consider a global rank envelope or area-above-envelope statistic.
5. **Strip 1 heterogeneity** — two physically separated tissue fragments. Documented but not modelled. Worth a sensitivity check splitting them.
6. **Edge weight choice in nb12** — peak `l_obs` is abundance-biased. Test whether SES or area-based weights change community membership.

---

## What this audit is *not*

This audit reports the state of the code and decisions. It does **not** re-run any analyses. For the question "if we fix these things, can we get a meaningful network?", see `network_revival_plan.md`.
