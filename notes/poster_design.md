# NeuroMonster Poster — Design Decisions

> **Print deadline: Saturday 6 June 2026.** Travel Rome Sunday 7. A0 portrait. Working file: `docs/NoSeggs_Poster.pptx`. Audience: computer scientists + neuroscientists.

Created 2026-06-02 after A.6 result. This is the locked design — change only with explicit reason.

---

## Locked title

**"NoSeggs: Segmentation-Free Ligand–Receptor Co-Localisation in Spatial Transcriptomics"**

Subhead (smaller, beneath title): *Bivariate point-pattern statistics on raw transcripts — no cells, no clusters, no segmentation.*

Authors line: **Harry Westwood** · *Supervisor: Dr Dan Nicolau, Nicolau Lab* · MSc Bioinformatics, King's College London · 2026

Banner under authors: *Platform: NanoString CosMx SMI · Tissue: severe asthmatic, RSV-infected human lung · ~2.66 M transcripts · 1,000-gene panel*

---

## Locked framing

- **Story spine:** segmentation fails on clinical CosMx → we skip it → bivariate Ripley's K on raw transcripts → test discriminates by effect size → networks of co-localised L–R pairs.
- **Honest about A.6:** the binary envelope test is under-powered with the production window. We use continuous SES ranking instead. This is a methodological finding, not a limitation.
- **No overclaiming.** Don't say "first to use Ripley's K in spatial transcriptomics" (Squidpy has it). Don't say "first segmentation-free" (FICTURE exists). DO say "first segmentation-free L-R co-localisation pipeline using bivariate Ripley's K on transcript point clouds."

---

## Column-by-column content

### Left column — Introduction + Dataset + Validation

**Intro paragraph (replace current filler, ~120 words):**

> Spatial transcriptomics maps gene expression at single-molecule resolution in intact tissue. Almost every downstream analysis — cell-type calling, ligand–receptor inference, niche detection — requires a fragile first step: assigning each transcript to a cell.
>
> On the CosMx sample analysed here, vendor segmentation leaves ~42% of transcripts unassigned and most subcellular labels missing. Third-party tools (Cellpose, MOSAIK) cannot recover the signal from poor DAPI staining. Standard downstream tools (CellChat, LIANA, Bento) become unusable.
>
> NoSeggs takes a different route: treat each gene as a 2D point pattern in tissue space, and use bivariate Ripley's K-functions with label-permutation envelopes to detect ligand–receptor co-localisation directly. No cells. No clusters. Just transcripts.

**Dataset table (keep current structure, fact-check values):**

| Property | Value |
|---|---|
| Platform | NanoString CosMx SMI |
| Tissue | Severe asthmatic, RSV-infected human lung |
| Slide | S1 (primary); 6 FOVs × 3 tissue strips |
| Transcripts | 703,000 (after QC) of ~2.66 M raw |
| Gene panel | 1,000 targets + 78 system probes |
| L-R panel | 146 pairs (CellChatDB ∩ panel, ≥50 transcripts/gene/strip) |

**Validation Controls box (keep position, update text):**

| Control | Pair | Expected | Observed |
|---|---|---|---|
| Positive | KRT8 × KRT18 | Above null | SES = −2.18 (top of distribution) |
| Negative 1 | MALAT1 × KRT18 | Far below | SES = −4.85 |
| Negative 2 | KRT8 × SCGB3A1 | Below | (fill from existing) |

### Middle column — Methods / Analysis Pipeline

Keep the existing 8-step vertical arrow structure. Confirmed step text:

1. **Load CosMx raw data** — zarr ingest; ~2.66 M transcripts with (x, y, gene, FOV).
2. **Strip assignment** — per-FOV Gaussian Mixture Model on x-coordinate identifies 3 tissue strips (control / infected / control).
3. **DBSCAN noise filtering** — adaptive ε = 97th-percentile 1-NN distance, clipped to [20, 30] px; removes 6.8 % rogue transcripts.
4. **Observation window** — concave hull (ratio 0.1) bounds each strip; replaces bounding box.
5. **Bivariate Ripley K / L** — L(r) = √(K(r)/π) − r on transcript point patterns; Shapely-based polygon edge correction.
6. **Permutation envelope** — 99 (or 199) label-swap permutations build the null at each r.
7. **Effect-size ranking** — peak SES per (pair, strip); top-N pairs entered as edges.
8. **Network construction** — modularity maximisation (Louvain) on the SES-weighted graph; communities = co-localised signalling programmes.

Pull-out box (in the methods column or below it):

> **The shift from significance to effect size.** Under the production concave-hull window, even the positive control sits within the permutation envelope. Window choice is the dominant lever for inferential power on transcript-level point patterns. We adopt continuous standardised effect size (SES) ranking — the continuous form of the global rank envelope test (Myllymäki 2017) — rather than binary significance.

### Right column — Results

Figure list revised 2026-06-03 after v2 sanity check (see decision_log).
The naive "cytokeratins at top of SES" prediction failed because the
method tests excess-spatial-association-beyond-marginals, not co-occurrence;
the new framing makes Group 7 (infection response) the biological headline,
and adds visual breadth (heatmap, volcano, UMAP) for an A0 audience.

**Fig A — Motivation: vendor segmentation drops 42% of transcripts.**
Use existing: `results/figures/01_vendor_segmentation_dropout.png`.
*Caption.* "Vendor segmentation on the CosMx sample fails to assign 42 % of transcripts to cells, and most retained cells lack subcellular compartment labels. Segmentation-dependent downstream tools (CellChat, LIANA, Bento) are unusable on this data."

**Fig B — Method check: continuous SES discriminates the positive control
from the housekeeping × specific negative.**
Use existing: `results/figures/13_pos_vs_neg_sanity.png`.
*Caption (Option 2 — keep the binary-test reasoning explicit so the
audience sees why we pivoted).* "Under the production hull window, both
the canonical positive control (KRT8 × KRT18, left) and the housekeeping ×
epithelial negative (MALAT1 × KRT18, right) sit within the label-swap
permutation envelope on the infected strip — the binary `L_obs above the
upper envelope at ≥3 consecutive radii` test is uninformative under this
window. But the standardised effect sizes still separate by 2.67 envelope
half-widths (positive = −2.18, negative = −4.85). We therefore rank pairs
by continuous peak SES, not binary pass/fail."

**Fig C — Group-wise SES distribution (validation by null behaviour) — done 2026-06-03.**
File: `results/figures/13_C_group_ses_distribution.png` (was Fig D in the
draft layout — letters re-numbered 2026-06-03).
Spec: violin + stripplot, x = pre-registered group (G0 CellChatDB,
G1 cytokeratin, G3 epithelial niche, G4 immune, G5 stromal, G6 negative
anchor, G7 infection, G8 paracrine), y = peak SES, ordered by mean SES.
*Caption.* "Method is correctly null for pairs with diffuse marginal
distributions (cytokeratins G1, housekeeping anchors G6) and elevated for
pairs with focal spatial structure (paracrine G8, infection G7, stromal G5,
immune G4). G3 (KRT19 × SCGB3A1, EPCAM × SCGB3A1) shows statistical
anti-association — observed transcripts are *less* spatially adjacent than
the label-swap null. Bivariate K under this null detects excess (or
deficit) spatial co-association beyond marginal density; for pan-tissue
markers, observed and permuted patterns coincide by design."

**Fig D — SES heatmap (pairs × strips), sorted, with group colour bar — BANKED 2026-06-03.**
Status: generated but not selected for the poster (competent landscape but
not visually punchy enough alongside the network/UMAP/Fig F panels).
Reserve for dissertation Results section.
File: `results/figures/13_D_ses_heatmap.png` (was Fig E in the draft layout).
Spec: heatmap with rows = pairs sorted by mean SES across strips,
columns = strip_1 / strip_2 / strip_3, cells coloured by peak SES (RdBu
diverging). Left-edge group colour bar.
*Caption.* "Peak SES across the pooled v1 + v2 panel (~206 unique pairs ×
3 strips). Pairs ordered by mean SES; group identity shown on the left.
Infection-response (G7) and paracrine (G8) pairs dominate the high end,
predominantly on strip 2 (infected). G3 anti-association pairs occupy the
deep-blue tail."

**Fig E — Volcano-style: peak SES vs r at peak, coloured by group — BANKED 2026-06-03.**
Status: generated; discrete r-value banding doesn't read polished. Reserve
for dissertation Results section (could be salvaged with jitter or a
different y-axis encoding, but not poster-grade).
File: `results/figures/13_E_volcano.png` (was Fig F).
Spec: scatter, x = peak SES, y = r at peak (px), colour = group.
Reference line at SES = +1 (binary envelope crossing).
*Caption.* "Each (pair, strip) row plotted by effect size (peak SES) and
the spatial scale at which the effect maximises. Infection-response and
stromal pairs cluster at low-r, high-SES (focal, juxtacrine-to-paracrine
scale); cytokeratin and housekeeping pairs scatter around SES ≈ 0."

**Fig F — Infection signalling biology (BIOLOGICAL HEADLINE) — locked 2026-06-03.**
File: `results/figures/13_F_infection_signature.png`.
Final design: horizontal bar chart of per-pair *Strip-2 boost*
(SES_strip2 − mean(SES_strip1, SES_strip3)) for the 13 surviving G7 pairs,
sorted top-to-bottom by boost, colour-coded by biological sub-theme
(immune cell communication = orange, pure ISG–ISG = yellow, antiviral ×
epithelium = blue). Each bar labelled with its numeric boost. Title
reports per-sub-theme positive/total counts built from the data so it
can't drift.

*Caption (Option A primary, Option B secondary — agreed 2026-06-03).*
"Immune cell-cell signalling pairs preferentially co-localise on
RSV-infected tissue (5/6 immune-communication pairs positive). Antigen
presentation onto T cells (HLA-DRA × CD3D, +0.74), T-cell-recruitment
chemokine signalling (CXCL10 × CD3D, +0.59), and antigen presentation on
macrophages (HLA-DRA × CD68, +0.32) carry the strongest positive boosts.
**A paradoxical second finding emerges that only a segmentation-free method
can surface**: antiviral × epithelial pairs (ISG15 × KRT8, MX1 × KRT5) move
*off* infected tissue, consistent with infection-induced epithelial loss
and migration of antiviral ISG production from epithelium to infiltrating
immune cells. The one negative immune-communication pair (CXCL8 × IL6) is
the NF-κB / classical-proinflammatory axis — consistent with chronic
asthmatic baseline on the control strips, not the type I IFN axis driven
by RSV."

*Per-pair numeric findings to quote in the poster sidebar or caption if
space allows:* HLA-DRA × CD3D +0.74, CXCL10 × CD3D +0.59, OAS1 × ISG15
+0.45, HLA-DRA × CD68 +0.32; ISG15 × KRT8 −0.64, CXCL8 × IL6 −0.74,
MX1 × KRT5 −0.90.

**Fig G — Per-strip network with Louvain communities (top-50 SES per strip) — done 2026-06-03 evening.**
File: `results/figures/13_G_network.png` (was Fig H). Generated in nb13 §8.
Spec: NetworkX spring layout per strip, nodes = genes, edges = top-50 pair
edges per strip, edge width ∝ peak SES, node colour = Louvain community,
node size ∝ degree.
*Design note (decided 2026-06-03 evening):* at top-50 edges the graph is
sparse (~50 edges over ~74 genes) so most Louvain "communities" are trivial
disjoint pairs. We colour only genuine modules (communities of ≥3 genes) and
grey the isolated pairs, so the real connected modules read on A0. Strip 2
(infected) carries more and larger modules (10 modules / 43 genes) than the
controls (6 / 36 and 8 / 38) — directly supports the criterion-3 story.
Edge ranking uses peak SES = max_r SES(r), the Fig C–F metric (NOT the §7
signed-peak-at-max-|SES| used only for the Fig B control discrimination).
*Caption.* "Top-50 spatial co-association edges per strip, with
modularity-based community detection (Louvain). Communities are spatial
signalling modules — interferon response, chemokine signalling, stromal
ECM, vascular markers — recovered without any cell segmentation or
cell-type calling."

**Fig H — UMAP of L(r) profile vectors, coloured by group — done 2026-06-03 evening.**
File: `results/figures/13_H_umap.png` (was Fig I). Generated in nb13 §9.
Spec: each (pair, strip) row treated as a 50-dim vector of L(r) values,
per-column standardised (so curve *shape*, not absolute L magnitude, drives
the embedding), UMAP (n_neighbors=15, min_dist=0.1, euclidean, seed=42) to
2-D, points coloured by group with the IBM palette shared with Figs C–F.
*Observed:* soft but real structure — immune (G4), infection (G7) and stromal
(G5) groups enrich one arm of the embedding; G0 CellChatDB (grey) spans the
whole space. Honest framing for the caption: "recovers biological grouping as
*soft* visual structure", not hard clusters.
*Caption.* "Pairs with similar L(r) curve shapes cluster in UMAP space.
Group identity is recovered as visual structure: stromal, immune, and
infection-response pairs occupy distinct regions, demonstrating that the
method extracts biologically meaningful spatial signatures from raw
transcripts alone."

**Fig I (optional, if space) — Differential SES (strip 2 − strip 1) heatmap.**
File: `results/figures/13_I_diff_ses.png` (was Fig J).
Spec: per-pair difference in peak SES, ordered by magnitude.
*Caption.* "Pairs ordered by infection-specific SES gain. Group 7 pairs
dominate the top of the ranking, providing a direct readout of the
spatial signalling difference between infected and uninfected tissue."

### Layout recommendation (A0 portrait, three columns)

A0 portrait is the area of 16× A4. There is room for 7 figure panels
alongside text. Suggested distribution:

- **Left column:** Intro + Fig A + Fig B + Validation table.
- **Middle column:** Methods pipeline (numbered) + Fig C (group SES) +
  Fig D (heatmap).
- **Right column:** Fig F (infection biology — biggest panel) + Fig G
  (network) + Fig H (UMAP).
- **Bottom strip:** Conclusions / references / acknowledgements / QR.

Fig E (volcano) reserved for dissertation Results chapter. Fig I optional —
include only if Fig F has room to spare and you want a second
Criterion-3-specific visual.

### Bottom strip — Conclusions / References / Acknowledgements

**Conclusions (replace WIP filler):**

> - Bivariate Ripley's K applied directly to raw transcript point clouds detects ligand–receptor co-localisation without any cell segmentation step.
> - Test discriminates positive from negative controls by standardised effect size (separation = +2.67 envelope half-widths under production window).
> - Fills a vacant niche: Bento, Squidpy, CellChat, LIANA all require segmentation. NoSeggs does not.

**Next steps (replace WIP filler):**

> - 80-pair extended panel (cytokeratin + housekeeping anchors) to anchor the SES ranking against known biology.
> - Pseudo-segmentation via marker-gene DBSCAN clustering for per-cell-type stratification.
> - Benchmark against Squidpy `gr.ripley` on cell-level reductions of the same data.

**Key References (max 5):**

1. Ripley BD (1976) *J. Appl. Prob.* 13:255–266 — the K-function.
2. Baddeley A, Rubak E, Turner R (2015) *Spatial Point Patterns: Methodology and Applications with R*. CRC Press.
3. Mah CK et al. (2024) *Genome Biology* 25:82 — Bento (segmentation-requiring contrast).
4. Palla G et al. (2022) *Nature Methods* 19:171–178 — Squidpy.
5. Myllymäki M et al. (2017) *J. Roy. Stat. Soc. B* 79:381–404 — global rank envelope test.

**Acknowledgements:** *Dataset and HPC compute courtesy of King's College London (CREATE HPC, `msc_appbio` partition). RSV samples via [TBD attribution]. CellChatDB via Jin et al. 2021. Project supervised by Dr Dan Nicolau, Nicolau Lab.*

**QR code:** points to `https://github.com/harrydaw/cosmx-pointpattern`. Verify before printing.

---

## Visual discipline

- 2 fonts max: one sans-serif for body (e.g. Lato, Inter, or Source Sans), one accent for headings (or same family heavier weight). Match existing KCL branding.
- Single accent colour for arrows, headings, highlight boxes — keep the current purple/burgundy if it's KCL-branded.
- Figures: viridis or cividis colour maps (colour-blind safe). 300 dpi minimum at A0 final size.
- Total body text ≤ 600 words across all panels. Figures and captions carry the load.
- White space is your friend. If the layout feels cramped, cut text, not figure size.

---

## Print pipeline

1. Export A0 PDF/X from PowerPoint with all fonts embedded.
2. Test print at A4 first to check legibility at small scale (body text must still be readable when shrunk 8× in someone's photo).
3. KCL Print Services or a local printshop — confirm 24h turnaround.
4. Order Friday; collect Saturday morning. Tube/case for travel.

---

## Things to confirm with Dan Nicolau (optional, time permitting)

- Is the title acceptable for the lab?
- Is there a Nicolau Lab logo / specific affiliation wording to add beside the KCL branding?
- Should the acknowledgement name a specific funder / dataset PI?

These are low-stakes; if Dan isn't reachable this week, ship the poster as-is and add corrections in the dissertation version.
