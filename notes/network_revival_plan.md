# Network revival plan

> Companion to `notebook_audit.md`. Goal: extract at least one defensible result from the network-analysis stage, or fail to do so with confidence.

## What we're trying to prove

Any one of the following counts as a meaningful, dissertation-worthy result:

1. **Criterion 1 — Communities without segmentation.** The algorithm detects gene-transcript communities directly from spatial point patterns, without per-cell aggregation.
2. **Criterion 2 — Where and at what scales.** Communities can be localised in tissue space and assigned an interaction radius (juxtacrine / paracrine / tissue-level).
3. **Criterion 3 — Infected vs control.** A detectable difference exists between strip 2 (infected) and strips 1 and 3 (controls), at edge, community, or scale level.

Each step below is tagged `[1]`, `[2]`, `[3]`, or combinations.

## Where we are now

`results/lr_panel_results.parquet` (438 rows = 146 CellChatDB pairs × 3 strips, `n_sim=199`):

- Only 7 rows have `l_obs > l_hi` at any single radius
- Only 1 row passes the `n_consec=3` filter: **SPP1 × CD44, strip 2**, peak excess ≈ 8.3 px
- Mean SES score ≈ −2.32 → the bulk of pairs sit *below* the envelope midline

This is the empty network. The plan below tries to figure out why.

---

## Step A — Cheap local diagnostics (no HPC, ~hours)

Run these in a new notebook `nb13_diagnostics.ipynb` (or add cells to nb12). All work against the existing parquet — no recompute.

### A.1 [1, 3] Inspect the 7 "any-radius exceedance" pairs
Print the (ligand, receptor, strip) tuples plus `max(l_obs - l_hi)`. If they cluster biologically (signalling family, cell-type co-marker) or by strip, that's already informative.

### A.2 [1] Distance-to-envelope grid
Small-multiples plot: 146 panels (one per pair), each showing all three strips' `l_obs` vs the upper envelope. Visual scan for "just below" pairs — these indicate the envelope is too tight rather than no signal.

### A.3 [1] `n_consec` sensitivity sweep
Re-run nb12's filter for `n_consec ∈ {1, 2, 3, 5}`. For each, count pairs and run Louvain. Tabulate (n_consec, n_pairs, n_communities, Q).

### A.4 [1] More powerful significance statistics
Replace pointwise max-consecutive with one of:
- **Global rank envelope test (GET-style)**: rank observed against simulated curves jointly across radii. More powerful, controls family-wise error.
- **Area-above-envelope**: ∫(l_obs - l_hi)⁺ dr as the test statistic; permutation-rank it against the same integral on simulated curves.
- **Studentised max** of `(l_obs - midline) / halfwidth` across radii.

Implement at least *one*. The GET-style test is the published gold standard for spatial point-pattern envelopes.

### A.5 [1, 2] Effect-size-only network
Skip significance entirely. Take the **top 50–100 pairs by peak SES** (or by area-above-midline). Build the network. Report community structure and per-strip Q. This gives a defensible criterion-1 figure even if formal significance is unattainable.

### A.6 [1] Positive control retest at production `n_sim`
The single most informative cheap experiment. Re-run **KRT8 × KRT18** on strip 2 locally:
- Use the *production* concave-hull window, edge correction, and `r_vals`
- Run at `n_sim = 199` (matches production) and `n_sim = 99` (matches nb04 validation) side by side
- Plot both envelopes against the same observed L(r)

If the n_sim=199 envelope buries KRT8 × KRT18, the production envelope is materially over-strict.

### A.7 [2] r_at_peak distribution by strip, top-N effect size
Histogram of `r_at_peak` for the top-N pairs from A.5, faceted by strip. Even without significance, this directly answers criterion 2.

---

## Step B — Broaden the panel and re-run on HPC (~1 day)

**Decision: Option B-shortlist** (~80 pre-registered anchor pairs).

See `extended_panel_rationale.md` for the full list and justification. High-level composition:

| Group | Approx count | Purpose |
|---|---|---|
| Cytokeratin co-expression anchors (KRT family pairs) | ~10 | Strong positive control for [1] |
| Epithelial cell-type pairs (KRT × SCGB / MUC) | ~15 | Co-expression in specific epithelial niches |
| Immune cell-type co-markers | ~10 | Lymphocyte/myeloid niches |
| Stromal pairs (fibroblast, endothelial) | ~10 | Connective-tissue niches |
| Negative anchors (MALAT1×_, housekeeping pairs) | ~10 | False-positive rate calibration |
| Infection-response candidates (ISGs × epithelial) | ~15 | Criterion 3 hypotheses |
| Buffer slots for biologically motivated additions | ~10 | — |

Production settings for the rerun:
- `n_sim = 99` (matches the validation regime, *not* the prior production setting)
- Same window / edge correction / radii as the original
- Save outputs to `results/lr_panel_results_v2.parquet` (do **not** overwrite v1 — diff matters)

---

## Step C — Decision point

Depending on Step A and Step B outcomes:

**C.1** — *KRT8 × KRT18 survives at `n_sim = 199`* → production envelope is fine; the CellChatDB panel just doesn't co-cluster strongly. Pivot dissertation figures to: top-N effect-size network (A.5), r_at_peak distribution (A.7), positive-anchor validation (B). Criterion 1 is achieved descriptively; 2 partially; 3 likely unprovable from inferential testing alone.

**C.2** — *KRT8 × KRT18 fails at `n_sim = 199` but passes at `n_sim = 99`* → the production envelope was over-strict. **Re-run the full 438-pair batch at `n_sim = 99`** (this is the path most likely to deliver criteria 1, 2 *and* 3). Redo nb11 / nb12.

**C.3** — *Even with broadened panel and looser envelope, no communities emerge* → criterion 1 unprovable on this data with this method. Write up as a methodological contribution: pipeline works on positive controls but the lung tissue / panel combination is not amenable to community detection. Specify what data *would* be needed (denser panel, multi-sample replication).

---

## Step D — Network methods upgrades (regardless of A–C outcome)

These should be done before any final dissertation figure:

### D.1 [1] Louvain stability across 10 seeds
Run Louvain with seeds 0..9 on the chosen graph. Plot modularity Q distribution. Report **best-Q partition** as the main result; mention stability range. Cross-tabulate which gene pairs end up in the same community ≥ 9/10 runs — these are the "stable core" of each community.

### D.2 [1] Leiden as a cross-check
`networkx.algorithms.community.leiden_communities` is available in NetworkX 3.6.1. Run it; report whether community membership agrees with Louvain.

### D.3 [1] Edge-weight sensitivity
Build the graph three times: edge weight = peak `l_obs`, peak SES, area-above-midline. Check whether community membership is stable across choices. Stability across reasonable weight choices is itself evidence that the structure is real.

### D.4 [3] Differential edge analysis
The existing nb12 logic for **constitutive / control-only / infection-specific** edges is correct — it just needs non-empty input. Once Step B/C yields a populated network, this analysis runs unchanged.

---

## Step E — Honest verdict (1-page memo)

Once Steps A–D are complete, write a 1-page memo as the final section of the audit. Required content:

- Positive control recovery rate at production settings (KRT × KRT survives? at which `n_sim`?)
- Panel pass rate (under the chosen significance criterion, on the broadened panel)
- Modularity Q vs degree-preserving null (Louvain best-Q vs Q from 100 randomised graphs with the same degree sequence)
- Stability across {Louvain seeds, edge-weight choices, n_consec values, significance tests}
- For each of the three success criteria: **achieved / partial / not achieved**, with the supporting figure ID

Decision question for the memo: *Should this pipeline be carried forward outside this project, or was it a fun experiment?*

---

## Order of operations (suggested next-session ordering)

1. **A.6** — Positive control retest at `n_sim = 199`. *One notebook, no HPC.* Highest information per hour.
2. **A.3, A.5, A.7** — local diagnostics that produce dissertation-grade figures regardless of outcome.
3. **A.1, A.2** — context-gathering; cheap.
4. **A.4** — GET-style test; medium effort.
5. **B** — broadened panel rerun on HPC. Only after A.6 has clarified `n_sim`.
6. **D.1–D.4** — applied to whichever graph survives.
7. **E** — final memo.

If under time pressure, the minimal path is **A.6 → A.5 → D.1 → E**. The plan is intentionally modular: each step produces a standalone figure or table.
