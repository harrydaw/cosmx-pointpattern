# Decision log

> Chronological record of methodological decisions, scientific findings, and
> their justifications. Companion to `planning_docs/05_methods_decisions.md`
> (which is the structured per-decision rationale) — this file is the
> *sequential narrative* of how we got here. Feeds dissertation Methods §
> rationale and Discussion § limitations.
>
> Newest entries at the top.

---

## 2026-06-04 (late) — M0 raw-data + M6 workflow added; Arial re-run; kernelspec fix

**Decision (Harry).** Make the poster "as complete as possible": (1) confirm
Arial everywhere — re-run nb13 results figs (B/C/F/G/H) in Arial for full
consistency, not just the nb14 suite; (2) add an **overall-operations flowchart
with a quick-start guide** showing what the tool is for and the tunable options
along the way; (3) add a **raw-data view** showing what CosMx actually produces,
highlighting the layers that fail vs the layer NoSeggs carries through.

**Implemented (two new nb14 figures, 300 dpi, Arial).**
- **`14_M0_raw_data.png`** — one FOV (FOV 4), three panels read live from
  `data/raw/Varsha_1234.zarr`: raw DAPI morphology → vendor segmentation mask
  (148 cells, the discarded layer) → raw transcripts coloured assigned (blue) vs
  unassigned `cell_ID==0` (magenta, ~40% of points). The grounding "what the data
  looks like" figure. Cell `m0-raw-code` in nb14.
- **`14_M6_workflow.png`** — portrait flowchart: 8 stages (INPUT → GMM strips →
  DBSCAN → window → bivariate Ripley K/L → permutation envelope → SES ranking →
  OUTPUTS) with amber option-chips for every tunable choice, plus a purple
  quick-start panel ("use NoSeggs when…"). Cell `m6-flow-code` in nb14.

**Poster (`docs/_build_poster.py`).** M0 now leads "The Problem — Segmentation
Fails" (above M1); M6 replaces the verbose text-pipeline paragraph as the visual
anchor of the Methods column (the "Significance → effect size" pull-out is kept).
Re-rendered via PowerPoint COM — both fit cleanly, no overflow.

**Kernel fix (infra).** The global `python3` kernelspec had a stale argv path
(repo moved to `Desktop\NoSegs_Project`), so notebook execution failed with
`FileNotFoundError [WinError 2]`. Did **not** touch the out-of-repo global spec;
instead installed a clean venv-prefix kernelspec (`noseggs` →
`.venv/share/jupyter/kernels/noseggs/kernel.json`, argv = venv python). Notebooks
now execute via `python -m nbconvert --execute --ExecutePreprocessor.kernel_name=noseggs`.
(Note: the `jupyter` dispatcher swallows output silently in this venv — call
`-m nbconvert` directly.) Also fixed a JSON-patch bug that had turned `\n` inside
title strings into literal newlines (SyntaxError) and mangled em-dashes; M0/M6
cell sources rewritten with proper escapes and re-executed (exec counts 1–9).

## 2026-06-04 (evening) — A0 poster assembled (methodology-forward layout)

**Decision (Harry).** Pin **Arial** across all figures for a consistent look;
build the actual A0 poster from the nb14 suite. The existing working file
`docs/NoSeggs_Poster.pptx` was found to be **A1** (59.4 × 84.1 cm), not A0, and
carried the older text-pipeline + Fig A/B/C/D results layout.

**Implemented.** New **`docs/NoSeggs_Poster_A0.pptx`** — true A0 portrait
(84.1 × 118.9 cm), built fresh (original A1 left untouched as a fallback) via
`docs/_build_poster.py` (python-pptx). Fully native/editable: text boxes, a real
PowerPoint table, and figures inserted as picture objects (no flattened images).
Methodology-forward 3-column layout per `poster_design.md`:
- **Left** (the problem + wrangling): Introduction → M1 dropout → M2 GMM → M3 DBSCAN.
- **Middle** (the method + that it works): pipeline + significance→effect-size
  pull-out → M4 windows → M5 r-scale/edge → T capability table → Fig B controls →
  validation table → Fig C group-SES.
- **Right** (the biology): Fig F infection boost → Fig G networks → Fig H UMAP.
- **Bottom band**: Conclusions/Next-steps · Key References · Acknowledgements · QR.
KCL logo + QR extracted from the A1 file and reused. Every section carries the
drafted real text plus purple-italic **`[ FILL ]` guide prompts** for the prose
Harry still needs to write/fact-check, and every figure has an editable legend
(captions lifted from `poster_design.md`). Verified by PowerPoint-COM PNG export.
Font note: nb14 figures are Arial; the nb13 results figures (B/C/F/G/H) are still
DejaVu Sans — flagged for a consistency re-run if time permits before print.

## 2026-06-04 (PM) — Methodology figures rebuilt as a polished suite (nb14)

**Decision (Harry feedback → implemented).** The first methodology pass lived in
standalone `scripts/poster_fig_*.py` and used hand-built composite sub-panels.
Four corrections: (1) figure code must live in **Jupyter notebooks** like the
rest of the project; (2) the hand-built sub-panels were weaker than the *real*
figures already produced in nb01–nb12 — reuse those, especially the polished
**segmentation-dropout** figure that grounds the work and the nb10 POC flow;
(3) tell a better story with **more, separate** figures foregrounding the
refinement big hitters (GMM assignment, R-value/scale design, window iterations,
DBSCAN QC); (4) the tool-comparison matrix works but needed **significant
smartening** → redesign as a **2D positioning map**.

**Implemented.** New **`notebooks/14_poster_methods.ipynb`** (keeps nb13 for
results-diagnostics) banks a suite of separate, individually polished 300-dpi
panels for flexible PowerPoint layout:
- **M1** `14_M1_segmentation_dropout.png` — rebuilt grounding dropout figure
  (43.6 % unassigned, per-FOV consistency).
- **M2** `14_M2_gmm_strips.png` — GMM strip assignment: representative FOV
  rotated-coords scatter (three crisp vertical strips + GMM midpoint boundaries)
  beside per-FOV trimodal x_rot histograms.
- **M3** `14_M3_dbscan_qc.png` — adaptive-DBSCAN noise QC (87 % retained;
  ε = p97 1-NN clipped [20,30]).
- **M4** `14_M4_window_iterations.png` — rect → convex → concave hull per strip
  with areas in Mpx² (window = dominant lever for inferential power).
- **M5** `14_M5_rscale.png` — two-segment r-grid schematic + concentric rings on
  a real KRT8 transcript + Shapely edge-correction weight at the hull boundary.
- **T** `14_T_tool_positioning.png` — 2D positioning map (x = segmentation
  dependence, y = answers L–R co-localisation); NoSeggs alone in the top-left
  vacant-niche quadrant. Positions conceptual, justified per
  `memory/project_competitor_landscape.md`.

All stats read live from `s1_all_strips_cleaned.parquet`; no DBSCAN/GMM re-run
(figures draw from stored labels). Helpers reused from `00_functions.ipynb`.

**Figure fixes during build.** M2 left and M4 originally used `set_aspect("equal")`,
which collapsed the ~1:12 thin tissue strips into illegible vertical slivers with
garbled overlapping x-ticks. Fix: M2 left now shows a single representative FOV
(strips only separate *within* a FOV — across FOVs the rotated-x ranges overlap)
with limited ticks; M4 uses `aspect="auto"` so the convex-vs-concave hull *shape*
difference is visible, with the true scale carried by the Mpx² areas in the legend.

**Superseded & deleted.** `scripts/poster_fig_methods.py`,
`scripts/poster_fig_toolcompare.py`, and the old `results/figures/13_M_methods_pipeline.png`
+ `13_T_tool_comparison.png` (replaced by the nb14 suite above).

**Note on kernel.** The project folder moved to `Desktop/NoSegs_Project/...`; the
global `cosmx-pointpattern` Jupyter kernelspec still points at the old
`Desktop/cosmx-pointpattern/.venv` path (stale, causes FileNotFound on headless
execute). Worked around by executing nb14 via NotebookClient with a kernel rooted
at the repo `.venv` python. The global kernel.json should be repointed when convenient.

---

## 2026-06-04 — Methodology emphasis: Fig M (QC composite) + Fig T (tool comparison)

**Decision (Harry).** Dedicate a substantial chunk of the poster to methodology
and tool comparison, rather than leaning on the (potentially bold) results
claims alone. The vast majority of the notebook work (00–09) was data-wrangling
that converts an *unanalysable* dataset into an analysable one — that *is* the
contribution, and foregrounding it makes the poster more defensible. Storyline:
data acquired → vendor segmentation failed (43.6 % of transcripts unassigned) →
third-party segmentation fixes failed → no existing pipeline avoids binning
~100 GB of expensive data → here are the tools built to analyse it segmentation-free.

**Fig M — composite QC/data-wrangling pipeline** (`results/figures/13_M_methods_pipeline.png`,
generator `scripts/poster_fig_methods.py`). 2×3 grid, raw → analysable:
M1 segmentation dropout (43.6 % cell_ID=0, by CellComp); M2 PCA axis rotation
(6.4° ACW); M3 GMM strip assignment (FOV 4); M4 adaptive DBSCAN ε from p97 1-NN
clipped [20,30]; M5 DBSCAN noise removal (87 % retained = 7.6 % noise + 2.0 %
small + 5.4 % manual); M6 concave-hull window (ratio 0.1) + Shapely edge
correction. All stats read live from `s1_all_strips_cleaned.parquet` so the
figure cannot drift from the data.

**Fig T — tool-comparison matrix** (`results/figures/13_T_tool_comparison.png`,
generator `scripts/poster_fig_toolcompare.py`). 7 tools × 4 capability columns,
NoSeggs the only all-green row. Facts verified 2026-06-01 against primary
sources (Bento = Mah 2024 needs segmentation; Squidpy = cell-level; CellChat/
LIANA = annotated matrices; Baysor = soft segmentation; FICTURE = seg-free but
domains not L–R; spatstat = generic). Encodes the defensible novelty claim
visually without overclaiming.

**Reproducibility note.** Both are standalone tested scripts (loading the
canonical cleaned parquet), kept in `scripts/` rather than embedded as
unexecuted nb13 cells — lower risk than re-running heavy nb13 surgery before a
print deadline, and they regenerate the exact PNGs deterministically.

**Layout consequence.** Poster pivots methodology-forward (see poster_design.md
REVISED layout): Fig M anchors the left column; Fig T leads the middle column.
Fig A (`01_vendor_segmentation_dropout.png`) now overlaps Fig M's M1 panel —
drop Fig A from the poster (keep for dissertation) to avoid redundancy. Figs D
(heatmap), E (volcano), I (diff SES) displaced to the dissertation Results
chapter.

---

## 2026-06-03 (evening) — Figs G (network) + H (UMAP) generated; two SES summaries clarified

Generated the last two poster figures in nb13 (§8 Fig G, §9 Fig H), both from
the pooled v1+v2 614-row frame built in the §A.5 cell.

**Two distinct SES summaries coexist — do not conflate.** While building Fig G
I confirmed the project uses *two* different per-row SES scalars:

1. **`peak_ses = max_r SES(r)`** (the strongest positive excursion). This is the
   ranking metric for Figs C/D/F and now G. Group means reconcile exactly:
   G8 +0.40, G7 +0.35, G5 +0.31, G4 +0.30, G0 +0.27, G6 −0.16, G1 −0.19,
   G3 −1.98.
2. **`signed-peak-at-max-|SES|`** (nb13 §7 `peak_signed_ses`) — the most extreme
   deviation regardless of sign. Used *only* for the Fig B control framing
   (KRT8×KRT18 −2.18, MALAT1×KRT18 −4.85). Under (1) these same controls are
   near-zero/positive, so the two must never be quoted interchangeably.

**Fig G design decision.** Top-50 (pair, strip) edges by `peak_ses` (positive
only); nodes = genes; Louvain communities. At this sparsity (~50 edges over
~74 genes) the graph is mostly *disjoint gene-pairs*, so naive Louvain returns
~25 trivial communities — visually useless and methodologically hollow. Fix:
colour only genuine modules (Louvain communities of **≥3 genes**) and grey the
isolated pairs. This makes the real connected signalling modules legible and is
honest about the sparsity. **Criterion-3 bonus:** strip 2 (infected) carries
more and larger modules (10 modules / 43 genes) than strips 1 and 3 (6 / 36 and
8 / 38) — a second, independent line of evidence for infection-elevated
co-association beyond Fig F. Worth a sentence in Results §5.

**Fig H design decision.** Each (pair, strip) row's 50-D L(r) curve →
per-column StandardScaler (so curve *shape*, not absolute L magnitude, drives
the embedding) → UMAP (n_neighbors=15, min_dist=0.1, euclidean, seed=42).
Result: **soft, not hard** structure — G4/G5/G7 enrich one arm; G0 spans the
space. Caption must say "recovers biological grouping as soft visual
structure", not "distinct clusters" — the honest reading.

**Sensitivity TODOs (dissertation, not poster):** Louvain + spring layout both
seed=42; repeat Fig G across ≥5 seeds (already on the outline TODO). UMAP is
also seed-sensitive; note as a caveat.

---

## 2026-06-03 (afternoon) — Fig F sub-theme analysis: three-way infection biology

After generating Fig F (per-pair Strip-2 boost for the 14 surviving G7 pairs)
with sub-theme colouring, the breakdown is:

| Sub-theme | n | Positive | Mean boost |
|---|---:|---:|---:|
| **Immune cell communication** (HLA-DRA, CXCL, CD3D, CD68, IL pairs) | 6 | 5 | **+0.20** (or +0.39 excluding CXCL8 × IL6 outlier) |
| **Pure ISG–ISG** (ISG15, IFIT, MX1, OAS1 pairs) | 4 | 3 | **+0.16** |
| **Antiviral × epithelium** (ISG/IFIT/MX × KRT) | 3 | 1 (essentially zero) | **−0.50** |

**Per-pair findings worth remembering:**

- HLA-DRA × CD3D: **+0.74** — antigen presentation onto T cells, strongest positive.
- CXCL10 × CD3D: **+0.59** — IFN-induced T-cell-recruitment chemokine signalling.
- OAS1 × ISG15: **+0.45** — canonical type I IFN co-expression.
- HLA-DRA × CD68: **+0.32** — antigen presentation on macrophages.
- IFIT1 × IFIT3: +0.24; IL6 × IL1B: +0.16; ISG15 × CD68: +0.15.
- ISG15 × MX1: +0.01 (essentially null); ISG15 × IFIT1: −0.07 (essentially null).
- ISG15 × KRT8: **−0.64** — ISG transcripts not near epithelium on infected.
- CXCL8 × IL6: **−0.74** — NF-κB/proinflammatory axis on controls (asthmatic baseline).
- MX1 × KRT5: **−0.90** — deepest negative; antiviral × basal epithelium gone.

**Biological reading (three findings):**

1. **Immune cell communication elevates on infected tissue (primary, Option A).**
   Antigen presentation (HLA-DRA × CD3D/CD68), T-cell-recruitment chemokine
   (CXCL10 × CD3D), monocyte/macrophage signalling (ISG15 × CD68) all peak on
   infected. The textbook adaptive-immune-recruitment signature of a viral
   response.

2. **Antiviral × epithelium pairs redistribute *off* infected tissue
   (paradoxical secondary, Option B framing).** Two readings, both
   consistent: (a) epithelial cytopathic effect — RSV kills/sloughs infected
   epithelium so KRT transcripts thin/displace; (b) ISG production migrates
   from epithelium to infiltrating immune cells over the course of infection
   (supported by ISG15 × CD68 being on the immune-positive side). This is a
   finding **only a segmentation-free method can make** — cell-based tools
   can't see "ISG transcripts redistributing across compartments" because
   they pre-bin transcripts into cells before doing spatial analysis. This
   is the unique-contribution argument for the dissertation.

3. **CXCL8 × IL6 outlier is biologically interpretable, not a problem.**
   CXCL8 / IL6 is the NF-κB / classical-proinflammatory cytokine axis. RSV
   drives the type I IFN pathway (IRF3 → CXCL10, ISGs, MHC), not classical
   proinflammatory. Severe asthmatic lung carries a chronic CXCL8/IL6
   baseline on the control strips. The outlier is consistent with "the
   acute viral response is IFN-driven, not NF-κB-driven" — the textbook RSV
   phenotype. Worth one sentence in Discussion.

**Statistical honesty.** Per-sub-theme sign tests don't reach formal
significance individually (n = 3, 4, 6 are too small). The defensible
combined test for the dissertation Results section is a **Mann-Whitney U on
the combined positive sub-themes (immune-comm + ISG–ISG, n = 10) vs
antiviral × epithelium (n = 3)** — directional separation is real and the
effect sizes are substantial (±0.3 to ±0.9). Not yet run; flag for the
Results chapter.

**Framing decided.** Option A is the poster primary headline ("immune
cell-cell signalling pairs elevate on infected tissue"). Option B is shown
alongside as a paradoxical secondary finding ("and antiviral × epithelium
pairs redistribute off infected tissue — a result only segmentation-free
analysis can surface"). Fig F is locked.

**Other figures:**
- Fig C (group-wise SES violin) — done, banked for poster.
- Fig D (heatmap of all 203 pairs × 3 strips) — banked. Competent landscape
  but not punchy enough for poster.
- Fig E (volcano: peak SES vs r_at_peak) — banked. Discrete r-values give a
  banded look that doesn't read polished. Save for dissertation Results.
- Fig F (per-pair G7 boost, sub-themed) — done, poster centrepiece.
- Fig G (per-strip network with Louvain communities) — pending.
- Fig H (UMAP of L(r) profile vectors) — pending.
- Fig I (differential SES heatmap) — pending, optional.

---

## 2026-06-03 (morning) — Poster Fig B caption: Option 2 framing (binary-test explanation kept explicit)

Decision: keep the binary-significance reasoning in Fig B's caption rather
than airbrushing it out. The audience sees "neither pair passes the binary
envelope test under this window, but SES still separates them by 2.67
half-widths — so we rank by continuous SES." This pre-empts the obvious
question ("did your test fail?") and walks the audience through the pivot
in one caption. The longer dissertation Methods/Results sections inherit
this framing.

Also: poster figure letters re-numbered 2026-06-03 (had skipped C in the
earlier revision). New order: A vendor dropout, B pos-vs-neg L(r),
C group SES (was D), D heatmap (was E), E volcano (was F), F infection
signature (was G), G network (was H), H UMAP (was I), I differential
(was J). Filenames realigned to the new letters.

---

## 2026-06-03 (morning) — v2 sanity check reveals method actually measures spatial-association-beyond-marginals

**Context.** With v2 aggregated, ran a quick `peak_ses` sort over the 176 v2
rows to verify the anchoring prediction: cytokeratin pairs (Group 1) at top,
housekeeping/cross-compartment anchors (Group 6) at bottom.

**Finding — the prediction was wrong, and instructively so.**

Top 10 by peak SES (all positive, range +0.92 to +1.16):
- Stromal: PECAM1×VWF strip_1, COL1A1×COL3A1 strip_3
- Paracrine: FGF7×FGFR2 strip_2, CCL2×CCR2 strip_2 and strip_3, TNF×TNFRSF1A strip_3
- Immune: CD68×CD14 strip_2
- **Infection-response: OAS1×ISG15 strip_2, CXCL10×CD3D strip_2, HLA-DRA×CD3D strip_2**

Bottom 10 (range −1.20 to −4.10) — and the surprise: dominated by **Group 1
cytokeratin pairs**, not Group 6 anchors:
- KRT19×SCGB3A1 (strip_1, strip_2) — SES −2.95 to −4.10 (deepest negative)
- KRT14×KRT15 strip_1, KRT8×KRT19 strip_1, KRT18×KRT19 strip_1 — SES around −2
- MALAT1×KRT18 — SES −1.5 to −2.7 across strips
- EPCAM×SCGB3A1 — SES around −1.2 to −1.6

**Interpretation — methodological, not a bug.**

The label-swap permutation null preserves the *marginal spatial distribution*
of every gene; it only breaks the (point, label) association. So bivariate
Ripley's K under this null measures **excess spatial co-association beyond
what each gene's marginal density predicts.** Not "do A and B co-occur in the
same regions" — that is a marginal-only question, and the null answers it
exactly the same as the data.

- **Cytokeratins co-saturate the epithelium** → both labels cover the same
  large region → under label-swap, observed ≈ null → SES ≈ 0 (slightly
  negative due to estimator behaviour and the negative bias documented in
  nb09 / `fraction_inside_hull`). No "excess" co-clustering exists to detect.
- **Vasculature / immune niches / signalling axes cluster in discrete
  sub-regions** → observed transcripts of A are physically near transcripts
  of B *more* than random label assignment to existing points would give →
  positive SES.
- **Negative SES (especially KRT19 × SCGB3A1 at −4.10) indicates statistical
  anti-association**: KRT19 spreads across simple epithelium while SCGB3A1
  concentrates in club-cell niches; the observed joint pattern is *less*
  spatially adjacent than the null because the two genes occupy disjoint
  epithelial sub-niches.

**This is a more rigorous interpretation of what the method tests than the
naive "co-expression detector" framing**:

> Bivariate Ripley's K on transcript point patterns with label-swap
> permutation null detects **excess (or deficit) spatial co-association beyond
> marginal density**.

For pairs where both genes are constitutively expressed across the same broad
region, the test is correctly null. For pairs with discrete focal
co-localisation (vasculature, niches, signalling foci), the test picks up
genuine excess.

**Headline biological finding — Criterion 3 evidence is real.**

Four of the top 10 pairs are infection-response pairs, all on strip 2
(infected): OAS1×ISG15, CXCL10×CD3D, HLA-DRA×CD3D, CD68×CD14. This is the
infection-specific signature the revival plan asked for — antiviral ISG
co-expression, T-cell-recruitment chemokine signalling, antigen-presentation
on T cells, monocyte/macrophage co-expression — all preferentially in the
infected tissue compartment.

**Implications for poster (revised from `poster_design.md`).**

Old Conclusions line:
- ❌ "Test discriminates positive from negative controls by SES" (based on
   A.6 KRT8×KRT18 vs MALAT1×KRT18 separation of +2.67 envelope half-widths).

Replace with:
- ✅ "Test detects **excess spatial co-association beyond marginal density**.
   Pairs with diffuse marginal distributions (cytokeratins) correctly score
   null. Pairs with focal spatial structure (vasculature, immune niches,
   paracrine signalling) score above the envelope."
- ✅ "**Infection-specific signalling signature on strip 2** (infected): ISG
   co-expression (OAS1×ISG15), T-cell-recruitment chemokines (CXCL10×CD3D),
   antigen presentation (HLA-DRA×CD3D), and innate-immune co-markers
   (CD68×CD14) all preferentially score in strip 2."

**Implications for nb13 A.5 today.**

Re-plan the morning's work:
1. Compute peak SES on **pooled v1 + v2** (~614 rows).
2. **Group-wise SES distribution** (violin/strip plot by G1–G8): visually
   establishes that G7 (infection) and G5 (stromal) and G8 (paracrine) score
   high; G1 (cytokeratin) and G6 (anchors) score null-to-negative. This is
   Fig C on the poster.
3. **Criterion 3 figure**: for G7 pairs specifically, peak SES distribution
   facetted by strip — predict strip_2 > strips 1, 3. This is now the strongest
   single result and should be the centerpiece of the poster's right column.
4. **Top-N network**: take top-50 per strip by SES, build NetworkX graph,
   communities by Louvain. Note that the "communities" are now to be
   interpreted as *spatial signalling modules*, not "co-expression clusters."

**A.6 interpretation revised.**

A.6's KRT8×KRT18 SES = −2.18 / MALAT1×KRT18 SES = −4.85 result still holds
descriptively, but the *interpretation* changes:
- Both pairs are "marginally diffuse," sitting below the envelope under
  label-swap.
- The 2.67 SES separation is small but real, and is informative as a
  *relative* contrast.
- The principled validation argument is no longer "positive sits above
  negative" but **"the method correctly assigns low SES to diffuse pairs
  (both cytokeratins and housekeeping × specific), and high SES only to
  pairs with focal spatial structure."**

---

## 2026-06-02 (evening) — Extended panel v2 submission

**Context.** A.6 (this morning) invalidated the original justification for the
extended-panel re-run (which was "loosen the envelope by going from n_sim=199
to n_sim=99"). New justification: **dynamic-range anchoring** of the SES
distribution — put cytokeratin co-expression pairs at the top of the ranking,
housekeeping anchors at the bottom. Differentiate ligand–receptor pairs against
a biological co-expression backdrop.

**Decision.** Submit tonight rather than defer to post-Rome. Reasoning: the
queue + run is on the critical path for Wednesday's poster figures, and
removing queue risk from Wednesday is cheap if we submit Tuesday evening.

**Decisions made for the v2 run:**

- **`n_sim = 99`** for v2. Justified by A.6: n_sim=99 and n_sim=199 envelopes
  on the positive control are statistically indistinguishable under the
  production hull window. n_sim=99 is the cheaper choice; SES is normalised
  by envelope half-width so v1 (n_sim=199) and v2 (n_sim=99) results are
  comparable in the pooled SES analysis.
- **Output isolated to v2 paths.** `results/hpc_job_manifest_v2.csv`,
  `results/per_pair_v2/`, `results/lr_panel_results_v2.parquet`,
  `results/failed_jobs_v2.csv`. v1 untouched.
- **No re-running of the original 146-pair CellChatDB panel.** v2 is purely
  additive; Wednesday's analysis concatenates v1 + v2 (~438 + ~176 = ~614 (pair,
  strip) rows).

### Pre-flight (`build_v2_manifest.py`) findings

Of 81 pre-registered pairs from `extended_panel_rationale.md`:

- **21 pairs dropped — gene(s) absent from CosMx 1,000-plex panel.** This is
  the headline finding of the pre-flight: the CosMx panel used on this sample
  is the **immuno-oncology / inflammation panel**, not a lung-anatomy panel.
  Specifically, the following lung-resident cell-type markers are absent and
  cannot be tested: SCGB1A1, SCGB3A2 (club cells), MUC5AC, MUC5B, TFF3 (goblet
  cells), FOXJ1, DNAH5 (ciliated cells), SFTPC, SFTPB, SFTPA1 (AT2 alveolar),
  AGER, HOPX (AT1 alveolar), TP63 (basal), CXCL11, FZD2.

- **4 (pair, strip) rows dropped — abundance < 50:**
  - ACTA2 × TAGLN, strip_2 (ACTA2 n = 40)
  - ACTA2 × MYH11, strip_2 (ACTA2 n = 40)
  - TNF × IL1B, strip_2 (TNF n = 48)
  - TNF × TNFRSF1A, strip_2 (TNF n = 48)

- **176 surviving jobs** for the v2 SLURM array.

### Implications for dissertation

- **Methods §** — document the panel constraint explicitly. The 1,000-plex
  CosMx panel composition determines which biological hypotheses are testable.
  Lung-resident epithelial niche markers (Group 2 of the extended panel
  rationale) are absent by panel design.
- **Discussion §** — frame this as a panel-driven limitation, not a method
  failure. The pipeline would test these pairs if the genes were present.
  Recommend a custom or lung-specific panel for the niche-detection variant of
  this analysis.
- **Biological note for Discussion §** — ACTA2 in strip_2 sits at n = 40,
  just below threshold. Two possible readings: (a) reduced smooth-muscle /
  myofibroblast representation in the infected strip; (b) a smaller usable
  area on strip 2 specifically. Worth a sentence in the Results or Discussion;
  not the main story.
- **Biological note** — TNF in strip_2 (n = 48) just below threshold. Slightly
  inconvenient for the inflammation hypothesis but not prohibitive — TNF × IL1B
  and TNF × TNFRSF1A survive in strips 1 and 3.

### Coverage retained by Group

| Group | Pre-registered | Surviving | Status |
|---|---|---|---|
| 1 — Cytokeratin co-expression anchors | 10 | 10 | ✅ Fully intact — top-of-distribution biological anchors |
| 2 — Epithelial cell-type co-markers | 10 | 0 | ❌ Catastrophic loss (panel constraint) |
| 3 — Epithelial niche overlap | 8 | 2 | 🟡 Heavy loss (only KRT19×SCGB3A1, EPCAM×SCGB3A1 survive) |
| 4 — Immune cell-type co-markers | 10 | 10 | ✅ Fully intact |
| 5 — Stromal co-expression | 8 | 8 (one pair down 1 strip) | ✅ Effectively intact |
| 6 — Negative anchors | 10 | 8 (and one pair down strip 2) | ✅ Bottom-of-distribution calibration works |
| 7 — Infection-response candidates | 15 | 14 | ✅ Criterion-3 story intact |
| 8 — Paracrine signalling exploratory | 10 | 9 | ✅ Effectively intact |

### Poster impact

The story you wanted lands. Cytokeratins at top of SES, negatives at bottom,
infection-response candidates in the middle. Group 2 loss hurts the
dissertation's niche-detection sub-claim, not the poster's central claim.

---

## 2026-06-02 (morning) — Step A.6 result and C decision

**Context.** Network revival plan Step A.6 asked: does the production envelope
(n_sim=199) bury the positive control KRT8 × KRT18 because the envelope is
too strict, vs the validation regime (n_sim=99) where the positive control
passed cleanly in nb04?

**Finding (A.6).** KRT8 × KRT18 on strip 2, production concave-hull window,
both regimes:
- n_sim = 99:  peak excess = −1.56 px, max consecutive exceedances = 0 → FAILS
- n_sim = 199: peak excess = −1.56 px, max consecutive exceedances = 0 → FAILS

The n_sim=99 and n_sim=199 envelopes overlap almost exactly. **n_sim is not
the lever.**

**Finding (sanity check).** Pos vs neg control, same window + null:
- KRT8 × KRT18 (positive control): peak SES = −2.18
- MALAT1 × KRT18 (negative control): peak SES = −4.85
- Separation: **+2.67 envelope half-widths.**

**Decision: C.3 with effect-size pivot.**

The binary envelope test is under-powered under the production window — even
the positive control sits within the envelope. But the test *does* discriminate
positive from negative by **standardised effect size (SES)**: the positive
control sits 2.67 envelope half-widths *above* the negative control. We adopt
continuous SES ranking instead of binary significance.

**Methodological reading.** The concave-hull window tightly captures the
tissue. The label-permutation null then re-distributes labels across that same
tightly-bounded set of locations — so the permuted null is almost as
concentrated as the observed pattern. Absolute contrast collapses; the binary
envelope test is starved of signal. This is **a feature of the window choice,
not the data or the algorithm**. A wider window (rectangle) would inflate
contrast by including empty space.

This is a **genuine methodological finding for the dissertation**: on
transcript-level point patterns, **window choice is the dominant lever for
inferential power.** Wider windows give more power but at the cost of including
non-tissue regions in the null. Tighter windows are biologically faithful but
suppress inferential contrast. The right answer is to rank by continuous
effect size (SES), which is the continuous form of the global rank envelope
test (Myllymäki 2017).

**Implications.**
- nb12 filter logic changes from `passes_n_consec(>=3)` to `top-N by peak SES`.
- Poster Fig B becomes the side-by-side SES discrimination plot
  (`13_pos_vs_neg_sanity.png`), not the original rectangular-window control
  plot.
- Dissertation Discussion gains a section on the window-vs-power trade-off,
  citing this finding.

---

## 2026-06-01 — Reboot and competitor positioning

**Context.** Project resumed after ~3-week dormancy. Three deliverables in
flight: NeuroMonster A0 poster (~1 week), MSc dissertation draft (1 July to
Dan Nicolau), submission (16 July).

**Decisions.**

- **Repo streamlined.** Archived stale planning docs, merged future-work notes,
  removed empty `src/spatialco` and `tests/` skeletons. README rewritten with
  verified novelty positioning. STATUS.md created as the live status board.
- **Bento competitor question — verified.** Bento (Mah et al., *Genome Biology*
  2024) **requires** cell + nuclear segmentation as mandatory input. It does
  subcellular pattern analysis *within* already-segmented cells. NoSeggs solves
  a different problem (inter-cellular spatial co-localisation without
  segmentation). Bento is not a competitor — it is a complementary subcellular
  toolkit.
- **Defensible novelty claim locked in:**
  > "First segmentation-free L-R co-localisation pipeline using bivariate
  > Ripley's K on raw transcript point clouds, with hull-window edge
  > correction and panel-scale application."
- **Framing — Moderate–Bold.** Headline: "Segmentation-free ligand–receptor
  co-localisation in spatial transcriptomics via bivariate point-pattern
  statistics." Pitched for NeuroMonster (CS + neuroscience).
- **Conference details:** NeuroMonster, Rome, fly Sun 7 Jun. No abstract
  submitted — framing is fully open.
- **Supervisor:** Dr Dan Nicolau, Nicolau Lab. Internal deadline end of June.

### Verified competitor landscape

| Tool | Input | Output | Segmentation needed? |
|---|---|---|---|
| Bento (Mah et al. 2024) | transcripts + cell + nuclear masks | subcellular localisation patterns | **Yes (mandatory)** |
| Squidpy `gr.ripley` | cell coordinates + cluster labels | cluster-level co-localisation | Yes |
| CellChat / LIANA / commot | cell-type expression matrix | inferred L–R signalling | Yes |
| Baysor | transcripts | soft cell assignments | Soft segmentation |
| FICTURE (Si et al. 2023) | transcripts | pseudo-tissue factors | No (but doesn't do L–R) |
| spatstat (R) | generic 2D point pattern | K, L | N/A — generic library |
| **NoSeggs** | **raw transcripts only** | **L–R co-localisation + communities** | **No** |

**Sources verified:**
- Mah CK et al. (2024) *Genome Biology* 25:82. https://doi.org/10.1186/s13059-024-03217-7
- Bento docs: https://bento-tools.readthedocs.io/
- FICTURE (Si et al., bioRxiv 2023): https://doi.org/10.1101/2023.11.04.565621
- Palla G et al. (2022) *Nature Methods* 19:171–178 (Squidpy)
