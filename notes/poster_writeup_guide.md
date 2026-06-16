# Poster write-up guide — reading list, key points, what to write

> **Purpose.** A single working sheet to take to Rome. For every poster block it
> tells you (a) what `[ FILL ]` prompts remain, (b) what to read, (c) the key
> points to extract, and (d) which dissertation section it feeds. Use it to
> finish the poster prose first, then roll the same material straight into the
> 12,000-word draft.
>
> **Companion docs:**
> - `notes/dissertation_outline.md` — section-by-section bullet outline (the
>   spine of the write-up; this guide cross-refs its section numbers).
> - `notes/intro_methods_reading.md` — full annotated reading list + key terms.
> - `notes/decision_log.md` — the *why* behind every methodological call.
> - `notes/poster_design.md` — locked figure list + captions.
>
> **Reference shorthand** used below (full citations in `intro_methods_reading.md`):
> Ripley 1976/1988 · Diggle 2003 · Baddeley 2015 · Myllymäki 2017 ·
> Jin 2021 (CellChatDB) · Dimitrov 2022 (LIANA) · Armingol 2021 ·
> Palla 2022 (Squidpy) · Mah 2024 (Bento) · Si 2023 (FICTURE) ·
> He 2022 (CosMx) · Moses & Pachter 2022 · Marconato 2024 (SpatialData) ·
> Blondel 2008 (Louvain) · Traag 2019 (Leiden) · Newman & Girvan 2004 ·
> Ester 1996 (DBSCAN).

---

## Order of attack (do this, in this order)

1. **Gather the RSV / lung-infection refs** — the one real gap (see §"Reference
   gap" at the foot). ~5 papers; everything else is already in
   `intro_methods_reading.md`.
2. **Fill the poster `[ FILL ]` prompts** — they are listed per block below.
   Mostly fact-checks + sign-offs, not new writing.
3. **Write the poster prose** block by block using the key points below
   (≤600 words total across the poster — figures carry the load).
4. **Roll into the dissertation** — each block maps to a dissertation section;
   expand the same points to the target word counts in `dissertation_outline.md`.

---

## Block-by-block

### A. Title / subhead / author line
- **`[ FILL ]`:** none (locked). Confirm funder/PI attribution wording (see
  Acknowledgements block).
- **Read:** —
- **Key points:** title is locked — *"NoSegs: Segmentation-Free Ligand–Receptor
  Co-Localisation in Spatial Transcriptomics."* Don't broaden the novelty claim.
- **Dissertation:** title page + running head.

### B. Introduction paragraph
- **`[ FILL ]`:** confirm the **42% dropout** headline number is the one you want
  to lead with (M1 panel reports the precise **43.6%**; the intro prose rounds to
  ~42%). Pick one and be consistent across intro + M1 caption + abstract.
- **Read:** Moses & Pachter 2022 (ST field overview); He 2022 (CosMx platform);
  Marconato 2024 (SpatialData, how the data was loaded); Jin 2021 + Armingol 2021
  (why L–R / cell-cell communication matters).
- **Key points to extract:**
  - One sentence on what spatial transcriptomics is + single-molecule resolution.
  - The segmentation bottleneck: nearly every downstream tool needs a per-cell
    matrix; segmentation depends on DAPI/membrane quality; fails on clinical tissue.
  - The hook: on *this* sample vendor segmentation leaves ~42% of transcripts
    unassigned → CellChat/LIANA/Bento become unusable.
  - NoSegs route: each gene = a 2D point pattern; bivariate Ripley's K + label
    permutation detects L–R co-localisation directly. "No cells. No clusters."
- **Dissertation:** Intro §1 (ST + CosMx, ~600w), §2 (segmentation problem,
  ~700w), §4 (L–R biology, ~500w).

### C. "The Problem — Segmentation Fails" (M1 dropout figure)
- **`[ FILL ]`:** none in prose; just keep the dropout % consistent with the intro
  (see B). Optionally one line: "consistent ~38–55% across all 12 FOVs."
- **Read:** the dropout literature note in `dissertation_outline.md` Intro §2
  (12–55% reported across imaging platforms — web-sourced 2026-06-01; find 1–2
  primary refs to cite for that range if you want it in the dissertation).
- **Key points:** vendor segmentation never assigns 43.6% of transcripts; loss is
  consistent across FOVs (not a one-FOV artefact); every cell-level tool silently
  discards them.
- **Dissertation:** Intro §2 + Results §1 framing; Discussion §1 (why standard
  tools fail here).

### D. Methods — the NoSegs pipeline (native flowchart M6 + M2–M5 panels)
- **`[ FILL ]`:** none (flowchart + captions locked). Editable now — tweak any
  bubble's text/colour directly in PowerPoint.
- **Read:** Ripley 1976 (K-function); Baddeley 2015 + Diggle 2003 (bivariate K,
  L-transform, edge correction); Ripley 1988 (permutation/edge); Myllymäki 2017
  (global rank envelope — the justification for continuous SES); Ester 1996
  (DBSCAN); Marconato 2024 (SpatialData ingest).
- **Key points (one per pipeline stage — these are your Methods subsections):**
  - **M2 strip assignment** — PCA aligns tissue axis; per-FOV 3-component GMM on
    rotated-x splits the 3 strips (control/infected/control); manual overrides
    where the mixture is ambiguous.
  - **M3 DBSCAN QC** — adaptive ε = 97th-pct 1-NN distance clipped [20,30] px;
    ~87% of transcripts survive (7.6% noise + 2.0% small + 5.4% manual).
  - **M4 observation window** — concave hull (ratio 0.1) vs convex vs bbox;
    concave shrinks area ~20× (e.g. strip 2: 313→61→14 Mpx²). **Window is the
    dominant lever for inferential power** — this is a genuine methods finding.
  - **M5 r-scale + edge correction** — two-segment r-grid (dense cell-scale
    [5,50]px ×30, coarse tissue (50,250]px ×20); Shapely disc-clipping edge
    correction; note the small negative bias that *cancels in the permutation*.
  - **Bivariate K → L(r)** [teal core] — K_AB(r), L(r)=√(K/π)−r, =0 under CSR.
  - **Label-swap permutation envelope** [teal core] — shuffle gene labels over
    existing locations; n_sim 99/199. **The central conceptual point:** this null
    preserves each gene's marginal distribution, so the test measures *excess
    co-association beyond marginal density*, NOT "do A and B occupy the same
    regions". Develop this carefully — it explains every later result.
  - **SES ranking** — peak SES = max_r (L_obs−mid)/half-width; continuous form of
    Myllymäki's global rank envelope test.
- **Dissertation:** Methods §1–§7 (the heaviest mapping — ~2,500w total). The
  flowchart stages are literally your Methods subsection headings.

### E. Validation & controls (Fig B + table)
- **`[ FILL ]`:** **RESOLVED — Negative-2 row dropped.** KRT8×SCGB3A1 was never
  run at panel scale (absent from both `lr_panel_results*.parquet`), so the table
  now shows Positive (KRT8×KRT18, SES −2.18) vs Negative (MALAT1×KRT18, SES −4.85).
  No further fill needed. If you *want* a third control for the dissertation,
  recompute KRT8×SCGB3A1 live via nb13 §7 (`peak_signed_ses`, the signed-peak-at-
  max-|SES| metric — NOT the parquet peak_ses) on the cleaned transcript parquet.
- **Read:** Myllymäki 2017 (why continuous SES is defensible); Baddeley 2015
  (envelope tests).
- **Key points:** under the tight production hull, even the positive control sits
  inside the envelope → the binary test is under-powered; but SES separates pos
  from neg by **+2.67 envelope half-widths**. This motivates the significance →
  effect-size pivot. Distinguish the two SES metrics (signed-peak for control
  discrimination vs max_r peak_ses for Figs C–H ranking).
- **Dissertation:** Results §1 (validation) + §2 (window as inferential lever);
  Discussion §4 (limitations: window, edge-correction bias).

### F. Results — infection biology (Fig F headline + Fig G network + Fig H UMAP)
- **`[ FILL ]`:** **`[ FILL: confirm headline wording ]`** under the purple
  take-home banner — sign off or tweak. Current draft: *"Immune cell-communication
  pairs preferentially co-localise on infected tissue — recovered with no
  segmentation step,"* with the **5/6** big-number callout and "+2.67 envelope
  half-widths" sub-line.
- **Read (THIS IS THE GAP — see foot of doc):** RSV / lung-infection biology —
  type I IFN response, ISGs (ISG15, MX1, OAS1), antigen presentation (HLA-DRA),
  T-cell recruitment chemokines (CXCL10), monocyte/macrophage markers (CD68),
  IFN vs NF-κB axis (CXCL8/IL6). Plus Blondel 2008 / Traag 2019 (Louvain/Leiden)
  for the network; Newman & Girvan 2004 (modularity Q).
- **Key points:**
  - **Primary (Fig F):** immune cell-cell signalling elevates on infected tissue —
    5/6 immune-communication pairs positive. Top boosts: HLA-DRA×CD3D +0.74,
    CXCL10×CD3D +0.59, HLA-DRA×CD68 +0.32, ISG15×CD68 +0.15. Textbook adaptive
    immune RSV response.
  - **Paradox (the unique contribution):** antiviral×epithelium pairs (ISG15×KRT8
    −0.64, MX1×KRT5 −0.90) move *off* infected tissue — consistent with epithelial
    cytopathic loss and ISG production migrating from epithelium to immune cells.
    **Only a segmentation-free method can see this** (cell-based tools pre-bin
    transcripts and cannot detect cross-compartment redistribution).
  - **CXCL8×IL6 outlier** (−0.74): NF-κB axis, not IFN — consistent with chronic
    asthmatic baseline on controls; RSV is IFN-driven. A *supporting* observation.
  - **Fig G network:** infected strip carries more/larger Louvain modules
    (10/43 genes vs controls 6/36 and 8/38) — a second, independent line of
    infection-difference evidence beyond per-pair SES.
  - **Fig H UMAP:** L(r) curve-shape vectors recover group structure as *soft*
    visual clustering (immune/infection/stromal arms) — honest framing, not hard
    clusters.
  - **Stat caveat:** per-sub-theme sign tests are under-powered (n=3,4,6); the
    coherent direction + effect sizes are the substance. **TODO before submission:**
    Mann-Whitney U on positive sub-themes (immune+ISG, n=10) vs antiviral×epi (n=3).
- **Dissertation:** Results §4 (infection, ~900w — the biological headline),
  §5 (network), §6 (UMAP); Discussion §3 (biological readings).

### G. Group-wise SES (Fig C — sits in the middle column)
- **`[ FILL ]`:** none.
- **Read:** (methods refs above; no new biology).
- **Key points:** method is correctly null on diffuse marginals (G1 cytokeratin
  −0.19, G6 negative anchor −0.16) and elevated on focal sub-themes (G8 paracrine
  +0.40, G7 infection +0.35, G5 stromal +0.31, G4 immune +0.30, G0 CellChatDB
  +0.27). G3 epi-niche −1.98 = real *anti-association* (KRT19/EPCAM × SCGB3A1
  occupy disjoint sub-niches). Bivariate K detects spatial *avoidance* too.
- **Dissertation:** Results §3 (~600w); Discussion §2 (what the test measures).

### H. Tool positioning (T capability matrix)
- **`[ FILL ]`:** none.
- **Read:** Mah 2024 (Bento — needs segmentation); Palla 2022 (Squidpy gr.ripley
  — cell-level); Si 2023 (FICTURE — segmentation-free but tissue domains, not
  L–R); Dimitrov 2022 (LIANA). Cross-check against `decision_log.md` 2026-06-01
  verified landscape before quoting positions.
- **Key points:** no published tool screens L–R co-localisation on raw transcripts
  without segmentation. NoSegs is the only all-✓ row across {segmentation-free,
  raw transcripts, L–R co-localisation, point-pattern statistics}. **Do NOT**
  claim "first segmentation-free ST" (FICTURE) or "first Ripley's K in ST"
  (Squidpy). DO claim the specific defensible sentence.
- **Dissertation:** Intro §5 (where this sits + aims); Discussion §1 (comparison
  with toolset).

### I. Conclusions / Next steps
- **`[ FILL ]`:** confirm the 3 conclusion bullets + 3 next-step bullets read
  correctly (they're drafted in `poster_design.md` Bottom-strip section).
- **Key points:** (1) bivariate K on raw transcripts detects L–R co-localisation
  with no segmentation; (2) discriminates controls by effect size (+2.67 half-
  widths); (3) fills a vacant niche. Next: extended panel anchoring; pseudo-
  segmentation via marker-gene DBSCAN; Squidpy benchmark.
- **Dissertation:** Discussion §5 (future work) — pull from `notes/future_work.md`.

### J. References / Acknowledgements / QR
- **`[ FILL ]`:**
  - **References** — confirm the 5 key refs (Ripley 1976; Baddeley 2015;
    Mah 2024; Palla 2022; Myllymäki 2017). Max 5 on the poster.
  - **Acknowledgements** — **`[ FILL ]` funder / dataset-PI attribution**:
    "RSV samples via [TBD]". Confirm with Dan Nicolau. KCL CREATE HPC
    (`msc_appbio`) + CellChatDB (Jin 2021) + supervisor already drafted.
  - **QR** — confirm target URL `https://github.com/harrydaw/cosmx-pointpattern`
    resolves and the repo is public before printing.
- **Dissertation:** references section (full, not capped at 5).

---

## Reference gap — gather these before writing Results/Discussion

Everything in §Methods/Intro/Tool-positioning is already covered by
`intro_methods_reading.md`. The **one real gap is RSV / lung infection biology**,
needed for Fig F's interpretation (Results §4 + Discussion §3). Gather ~5:

1. **Type I IFN response to RSV** — ISG induction (ISG15, MX1, OAS1) in airway
   epithelium / immune cells.
2. **Antigen presentation + T-cell recruitment** — HLA class II (HLA-DRA),
   CXCL10/CXCL11 as IFN-induced T-cell chemokines, CD3D⁺ T-cell infiltration.
3. **Monocyte/macrophage antiviral role** — CD68⁺ cells, ISG expression in the
   myeloid compartment.
4. **RSV cytopathic effect on epithelium** — sloughing / epithelial loss
   (supports the "antiviral moves off epithelium" reading).
5. **IFN vs NF-κB axis in RSV + asthma** — why CXCL8/IL6 is the chronic-asthma
   baseline, not the acute viral signal.

A good 2023–2025 review of RSV airway immunopathology likely covers 1–4 in one
citation; the IFN-vs-NF-κB point may need a dedicated ref.

---

## Definition of done (poster print + Rome write-up readiness)

**Poster — before print:**
- [ ] Sign off / tweak the Fig F headline wording (`[ FILL: confirm headline ]`).
- [ ] Pick one dropout figure (42% vs 43.6%) and use it consistently.
- [ ] Fill Acknowledgements funder/PI attribution; confirm with Dan Nicolau.
- [ ] Verify QR URL resolves + repo public.
- [x] Negative-2 control row dropped (was never run; not fabricated).
- [ ] Touch up spacing in PowerPoint (flowchart bubbles now individually editable).
- [ ] A4 test print → check legibility at ~8× shrink → A0 order at KCL Print
      Services (confirm 24h turnaround) → collect Sat 6 June.

**Write-up kit — ready to start in Rome:**
- [x] `dissertation_outline.md` (section spine + word budgets).
- [x] `poster_writeup_guide.md` (this file — reading + key points per section).
- [x] `decision_log.md` (the *why* for Methods/Discussion).
- [x] `intro_methods_reading.md` (annotated reading list + key terms).
- [ ] Gather the 5 RSV/lung-infection refs above.
- [ ] (Dissertation-phase TODOs, from `dissertation_outline.md`): Mann-Whitney U
      on G7 sub-themes; Louvain seed stability; Squidpy benchmark; global rank
      envelope supplementary check.
