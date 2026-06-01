# Extended LR panel — rationale

> Companion to `network_revival_plan.md` (Step B). ~80 pre-registered gene pairs to add to the 146-pair CellChatDB panel for the broadened HPC run.

## Why pre-register?

For criterion 1 ("the algorithm detects communities without segmentation"), the strongest evidence is: *we chose these pairs before looking at the data, on biological grounds, and the algorithm recovered the expected co-clustering.* Listing the pairs here, with a per-pair rationale, before running them is what makes that argument credible.

## Rules

- **All pairs must pass the same abundance filter** as the original panel (≥ 50 transcripts per gene per strip, see nb06). Each pair below needs verification against `viable_lr_pairs_all_strips.parquet`'s gene universe before commit.
- **All gene names must exist in the CosMx panel** used for this sample — verify against the raw zarr `target` column. Pairs whose genes are absent should be dropped (not silently failed at runtime).
- **Use `n_sim = 99`** for this run, matching the validation regime in nb04 (`run_array.sh` previously used 199 — change explicitly).
- **Save to `results/lr_panel_results_v2.parquet`** — do not overwrite the original. The diff is informative.

## Panel composition (~81 pairs)

Each pair tagged with:
- **Criterion(a)** it serves: `[1]` communities-without-segmentation, `[2]` scale localisation, `[3]` infected-vs-control
- **Expected scale**: J = juxtacrine (< 9 µm), P = paracrine (9–27 µm), T = tissue-level (> 27 µm), N = none expected
- **Expected outcome**: ⊕ co-clustering predicted, ⊖ no signal / depletion predicted, ? exploratory

### Group 1 — Cytokeratin co-expression anchors (n = 10) [1]
Cytokeratins form obligate heterodimers in epithelial cytoskeletons. They are the cleanest positive controls available in this panel. **All ⊕, scale J.**

| Source | Target | Rationale |
|---|---|---|
| KRT8 | KRT18 | Canonical type-I/II pair (validated in nb04). Must recover. |
| KRT5 | KRT14 | Basal epithelial obligate pair. |
| KRT5 | KRT15 | Basal layer co-marker. |
| KRT8 | KRT19 | Simple epithelium pair. |
| KRT18 | KRT19 | Simple epithelium pair. |
| KRT14 | KRT15 | Basal layer co-marker. |
| KRT5 | KRT17 | Basal / stress-induced co-expression. |
| KRT14 | KRT17 | Basal layer. |
| KRT19 | EPCAM | Pan-epithelial markers. |
| KRT8 | EPCAM | Pan-epithelial markers. |

### Group 2 — Epithelial cell-type co-markers (n = 10) [1, 2]
Pairs defining the same epithelial cell type. Expected co-clustering at the cellular scale (J → P boundary). All ⊕, scale J.

| Source | Target | Cell type |
|---|---|---|
| SCGB1A1 | SCGB3A1 | Club cells |
| SCGB1A1 | SCGB3A2 | Club cells |
| MUC5AC | MUC5B | Goblet cells |
| MUC5AC | TFF3 | Goblet cells |
| FOXJ1 | TUBB4B | Ciliated cells |
| FOXJ1 | DNAH5 | Ciliated cells (verify DNAH5 in panel) |
| SFTPC | SFTPB | Type II alveolar |
| SFTPC | SFTPA1 | Type II alveolar |
| AGER | HOPX | Type I alveolar |
| KRT5 | TP63 | Basal cells |

### Group 3 — Epithelial niche overlap (n = 8) [1, 2]
KRT × cell-type marker pairs — expected to co-cluster within specific epithelial niches. All ⊕, scale J/P.

| Source | Target | Rationale |
|---|---|---|
| KRT8 | SCGB1A1 | Simple epithelium hosts club cells |
| KRT18 | SCGB1A1 | Simple epithelium hosts club cells |
| KRT8 | MUC5B | Simple epithelium hosts goblet cells |
| KRT18 | MUC5AC | Simple epithelium hosts goblet cells |
| KRT19 | SCGB3A1 | Bronchiolar epithelium |
| EPCAM | MUC5B | Pan-epithelial × secretory |
| EPCAM | SCGB3A1 | Pan-epithelial × club |
| EPCAM | FOXJ1 | Pan-epithelial × ciliated |

### Group 4 — Immune cell-type co-markers (n = 10) [1, 3]
Same-cell-type markers; should co-cluster wherever immune infiltrate exists. **Criterion 3 candidate** if infiltrate is concentrated in strip 2. All ⊕, scale J.

| Source | Target | Cell type |
|---|---|---|
| CD3D | CD3E | Pan-T |
| CD3D | CD4 | T helper |
| CD3D | CD8A | Cytotoxic T |
| CD8A | CD8B | Cytotoxic T |
| CD4 | FOXP3 | Regulatory T (verify FOXP3 in panel) |
| CD68 | CD163 | M2 macrophage |
| CD68 | CD14 | Monocyte/macrophage |
| CD19 | MS4A1 | B cells |
| ITGAX | HLA-DRA | Dendritic cells |
| NCAM1 | NKG7 | NK cells (verify both in panel) |

### Group 5 — Stromal co-expression (n = 8) [1]
Fibroblast / endothelial / smooth-muscle marker pairs. All ⊕, scale J.

| Source | Target | Compartment |
|---|---|---|
| COL1A1 | COL3A1 | Fibroblast ECM |
| COL1A1 | COL1A2 | Collagen-I subunits |
| ACTA2 | TAGLN | Myofibroblast / smooth muscle |
| ACTA2 | MYH11 | Smooth muscle |
| PECAM1 | VWF | Endothelial |
| PECAM1 | CDH5 | Endothelial |
| VWF | CDH5 | Endothelial |
| DCN | LUM | Fibroblast proteoglycans |

### Group 6 — Negative anchors (n = 10) [1]
Cross-compartment pairs and housekeeping-anchored pairs. Used to calibrate false-positive rate. All ⊖ predicted (some may show weak ⊕ if compartments interleave), scale N.

| Source | Target | Reason |
|---|---|---|
| MALAT1 | KRT18 | Validated negative (nb05) |
| MALAT1 | CD3D | Housekeeping × T cell |
| MALAT1 | COL1A1 | Housekeeping × fibroblast |
| KRT8 | CD3D | Epithelial × T cell |
| KRT8 | COL1A1 | Epithelial × stromal |
| KRT8 | PECAM1 | Epithelial × endothelial |
| CD3D | COL1A1 | Lymphoid × stromal |
| MUC5AC | CD68 | Goblet × macrophage |
| SFTPC | ACTA2 | AT2 × smooth muscle |
| SFTPC | COL1A1 | AT2 × fibroblast |

### Group 7 — Infection-response candidates (n = 15) [3]
The criterion-3 workhorses. Pairs where co-clustering is expected to differ between strip 2 and the controls.

| Source | Target | Hypothesis | Predicted strip-2 vs control |
|---|---|---|---|
| ISG15 | IFIT1 | ISG co-expression in infected tissue | strip 2 ⊕ stronger |
| IFIT1 | IFIT3 | ISG co-expression | strip 2 ⊕ stronger |
| ISG15 | MX1 | ISG co-expression | strip 2 ⊕ stronger |
| OAS1 | ISG15 | ISG co-expression | strip 2 ⊕ stronger |
| ISG15 | KRT8 | ISG response in epithelium | strip 2 ⊕ |
| IFIT1 | KRT18 | ISG response in epithelium | strip 2 ⊕ |
| MX1 | KRT5 | ISG response in basal epithelium | strip 2 ⊕ |
| CXCL10 | CXCL11 | Chemokine co-expression | strip 2 ⊕ |
| CXCL8 | IL6 | Acute-phase inflammatory pair | strip 2 ⊕ |
| IL6 | IL1B | Pro-inflammatory pair | strip 2 ⊕ |
| TNF | IL1B | Pro-inflammatory pair | strip 2 ⊕ |
| ISG15 | CD68 | ISG response in macrophages | strip 2 ⊕ |
| CXCL10 | CD3D | T-cell-recruitment chemokine × target | strip 2 ⊕ |
| HLA-DRA | CD3D | Antigen-presentation × T cell | strip 2 ⊕ |
| HLA-DRA | CD68 | MHC-II on macrophages | strip 2 ⊕ |

### Group 8 — Exploratory / paracrine signalling (n = 10) [1, 2, 3]
Classic ligand-receptor paracrine pairs not always covered well by CellChatDB filtering, where co-clustering at scale P would be biologically meaningful. All ?, scale P expected.

| Source | Target | Hypothesis |
|---|---|---|
| TGFB1 | TGFBR2 | TGF-β paracrine signal → ECM remodelling |
| CXCL12 | CXCR4 | Stromal–lymphocyte chemotaxis |
| CCL2 | CCR2 | Monocyte recruitment |
| CCL5 | CCR5 | T-cell recruitment |
| IL1B | IL1R1 | Inflammatory paracrine |
| TNF | TNFRSF1A | Inflammatory paracrine |
| WNT5A | FZD2 | Stromal–epithelial signalling |
| FGF7 | FGFR2 | Mesenchymal–epithelial growth factor |
| HGF | MET | Stromal–epithelial regeneration signal |
| EGF | EGFR | Growth factor signalling |

---

## Total count and budget

| Group | Pairs |
|---|---|
| 1 — KRT co-expression | 10 |
| 2 — Cell-type co-markers (epithelial) | 10 |
| 3 — Epithelial niche overlap | 8 |
| 4 — Immune co-markers | 10 |
| 5 — Stromal co-expression | 8 |
| 6 — Negative anchors | 10 |
| 7 — Infection-response candidates | 15 |
| 8 — Paracrine signalling exploratory | 10 |
| **Total** | **81** |

At 3 strips × 81 pairs = **243 HPC jobs**. With per-job runtime 9–44 s observed for the previous run (production `n_sim = 199`), 99 permutations should be roughly halved — figure ~5–25 s/job. Total wallclock under the `%2` concurrency limit: a few hours.

## Pre-flight checks before submitting

1. **Verify gene availability**: load the raw zarr `target` column and `.isin()`-filter the panel above. Drop any unrecognised gene names with a printed warning.
2. **Verify abundance**: for each remaining pair, confirm both genes have ≥ 50 transcripts in each of the 3 strips. Drop pairs that fail.
3. **Build the v2 manifest**: write `results/hpc_job_manifest_v2.csv` with the surviving (pair, strip) rows.
4. **Confirm `--n_sim 99`** in the run command — do not inherit `199` from the previous `run_array.sh`.
5. **Output path**: write per-pair shards to `results/per_pair_v2/`, aggregate to `lr_panel_results_v2.parquet`. Never overwrite v1.

## What "success" looks like for this run

The minimum bar for declaring the rerun a success:

- **Group 1 (KRT co-expression)**: at least 7 of 10 pairs show `l_obs > l_hi` for ≥ 3 consecutive radii in at least one strip. If fewer, the envelope or window is the problem, not the panel.
- **Group 6 (negative anchors)**: at most 1 of 10 pairs shows the same. If more, the test has a false-positive problem.
- **Group 4 + 5 (cell-type co-markers)**: ≥ 50 % pass-rate is the threshold for "yes, the algorithm detects communities without segmentation" (criterion 1).
- **Group 7 (infection response)**: strip 2 pass-rate > strips 1, 3 pass-rate would directly evidence criterion 3.

If Groups 1 and 6 behave as predicted but Groups 4/5/7 don't, that's still a meaningful negative result for the dissertation: *the algorithm is correctly calibrated, the lung tissue / panel just doesn't carry the signal.* That conclusion is far more defensible than the current empty network.
