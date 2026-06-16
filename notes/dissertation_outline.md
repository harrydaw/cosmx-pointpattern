# Dissertation outline & key points

> Living outline to scaffold the 12,000-word MSc dissertation. Bullet form;
> turn into prose during the writing weeks (W3–W4 of the master plan).
> Each section heading lists the key points + which figures/data support them
> + which docs in the repo carry the source material.
>
> Cross-references:
> - `notes/decision_log.md` — chronological methodological + scientific findings (the *why* for each section)
> - `notes/poster_design.md` — current figure list with locked captions (Methods + Results bones)
> - `notes/poster_writeup_guide.md` — per-poster-block reading list + key points + what to write
> - `notes/notebook_audit.md` — per-notebook state
> - `notes/intro_methods_reading.md` — 120+ ref bibliography scaffold + Intro scaffolding + key terms
> - `notes/future_work.md` — Discussion future-work feed
>
> *(The old `planning_docs/04_background_topics.md` + `05_methods_decisions.md`
> no longer exist — Intro scaffolding folded into `intro_methods_reading.md`,
> per-decision rationale into `decision_log.md`.)*
>
> Target: ~12,000 words. Suggested split: 3,000 Intro / 2,500 Methods /
> 3,500 Results / 2,500 Discussion / 250 Abstract / ~250 acknowledgements +
> references inline.

---

## Introduction (target ~3,000 words)

### 1. Spatial transcriptomics and CosMx (~600 words)
- Single-molecule resolution gene expression in intact tissue.
- Imaging-based platforms (CosMx SMI, Xenium, MERSCOPE) vs sequencing-based
  (Visium, Slide-Seq). Why imaging gives single-transcript coordinates.
- CosMx specifics: 1,000-plex panel, sub-cellular resolution, vendor
  segmentation step.
- Reference data: `planning_docs/04_background_topics.md` covers this.

### 2. The cell segmentation problem (~700 words)
- Almost every downstream analysis needs a per-cell expression matrix:
  cell-type calling, ligand-receptor inference, niche detection.
- Segmentation is the bottleneck: relies on DAPI/membrane image quality,
  fails on overlapping cells, autofluorescent regions, noisy clinical samples.
- Quantitative literature: 12–55% transcript dropout reported across
  imaging-based platforms (from web search 2026-06-01).
- **Our dataset: 42% transcript dropout from vendor CosMx segmentation.**
  Third-party tools (Cellpose, MOSAIK) did not recover. → Fig A.
- Implication: segmentation-dependent downstream tools (CellChat, LIANA,
  Bento, CellPhoneDB) become unusable for L-R analysis on this sample.

### 3. Point-pattern statistics as an alternative (~700 words)
- 50 years of statistical point-pattern analysis (Ripley 1976, 1977;
  Diggle 2003; Baddeley, Rubak, Turner 2015).
- Bivariate Ripley's K-function: K(r) = expected number of B points within
  radius r of a randomly chosen A point, normalised by intensity of B.
- L-transform: L(r) = √(K(r)/π) − r, variance-stabilising; L(r) = 0 under CSR.
- Permutation envelopes (label-swap null) for significance.
- **Key insight to develop:** classical point-pattern is applied to *cells
  with type labels*. The contribution here is applying it directly to
  *transcripts with gene labels*, skipping the cell step entirely. This
  requires careful attention to what the label-swap null actually tests
  (excess co-association beyond marginal density — develop this in Methods).

### 4. Ligand-receptor signalling biology (~500 words)
- L-R interactions as the molecular basis of cell-cell communication.
- Curated databases: CellChatDB (Jin et al. 2021), CellPhoneDB, OmniPath.
- Cell-type-resolved L-R analysis tools all require segmentation +
  clustering as upstream input.
- The biological motivation for spatially resolved L-R: location matters
  (paracrine vs juxtacrine vs tissue-level), and physical proximity is
  evidence of plausible interaction.

### 5. Where this work sits + aims (~500 words)
- Verified competitor landscape (cite + summarise from `notes/decision_log.md`
  2026-06-01 entry):
  - **Bento (Mah et al. 2024 *Genome Biology*) — requires cell + nuclear
    segmentation; does subcellular pattern analysis.** Not a competitor.
  - Squidpy `gr.ripley` — operates on cluster-labelled cell coordinates.
    Cell-level, not transcript-level.
  - CellChat / LIANA / commot — all require cell-type expression matrices.
  - Baysor — produces soft cell assignments, still a segmentation product.
  - FICTURE (Si et al. bioRxiv 2023) — segmentation-free but outputs
    pseudo-tissue domains, not L-R.
  - spatstat (R) — generic library, not L-R-specific.
- **Defensible novelty claim:** "First segmentation-free L-R co-localisation
  pipeline using bivariate Ripley's K on raw transcript point clouds, with
  hull-window edge correction and panel-scale application." DO NOT broaden
  this — see decision_log 2026-06-01 for the verified language.
- Three aims (criteria from `notes/network_revival_plan.md`):
  - **A1:** detect gene-transcript communities from spatial point patterns
    without per-cell aggregation;
  - **A2:** localise them in tissue space and assign a characteristic
    interaction scale;
  - **A3:** detect differences between infected and control tissue
    compartments.

---

## Methods (target ~2,500 words)

### 1. Data and preprocessing (~400 words)
- CosMx SMI, severe asthmatic RSV-infected human lung sample S1.
- 6 FOVs (4, 5, 8, 9, 10, 11; expanded from initial 3), 3 tissue strips
  per FOV (control / infected / control). 702,873 transcripts total →
  611,150 (87%) kept after QC.
- Strip assignment via GMM on the PCA-rotated x-coordinate (`x_rot_px`;
  `reassign_strips_gmm`, per-FOV, k=3; rotation from `add_rotated_coords`).
- DBSCAN noise removal (canonical run in nb10, grouped by FOV on `x_rot_px`):
  adaptive ε = 97th-percentile 1-NN distance (sklearn `NearestNeighbors`)
  clipped to [20, 30] px; min_samples = 5; min_cluster_size = 150. 7.6% noise.
- Manual cleanup (nb10 `apply_cleanup`, min_cluster_size=120): 18 hand-picked
  outlier clusters + small clusters → 2.0% small + 5.4% manual.
- Cleaned output: `data/processed/s1_all_strips_cleaned.parquet`.
- NB: the K-function runs on raw `x_global_px/y_global_px` (rotation-invariant).
- Full reproducibility record: `Final_Writeup/METHODS_LOG.md`.
- Source for the rationale of each parameter: `planning_docs/05_methods_decisions.md`.

### 2. Observation window and edge correction (~400 words)
- Long narrow tissue strips → rectangular window inflates apparent area
  with empty space.
- Adopted concave hull (ratio = 0.1) per strip — 65–75% of convex hull area;
  43% of bounding-box area.
- Shapely-based polygon edge correction (`fraction_inside_hull`): for each
  point at radius r, compute fraction of the disc that lies inside the hull.
- Documented limitation (from nb09): introduces a small systematic *negative*
  bias in absolute K(r). Bias cancels in the permutation test because both
  observed and permuted realisations are biased equally — only envelope-
  relative significance is interpreted, not absolute L(r).

### 3. Bivariate Ripley's K and L (~400 words)
- Formal definitions; cite Ripley (1976), Baddeley et al. (2015).
- Implementation in nb00 / nb03; CSR validation in nb03.
- L-transform.

### 4. Label-permutation envelope (~300 words)
- Null hypothesis: gene labels are interchangeable across the existing
  point locations (= no gene-specific spatial association beyond marginal
  density).
- Procedure: shuffle gene labels uniformly across all points in the strip;
  recompute K(r) for the permuted labelling; repeat n_sim times.
- **Crucial conceptual point for the Discussion:** this null preserves each
  gene's marginal spatial distribution. So K (and L) under this null tests
  for *excess (or deficit) of joint clustering beyond what each gene's
  marginals already explain*. It does NOT test "do A and B occupy the same
  regions" — that question is answered identically under the null and the
  data, by construction. This reframes what positive vs negative SES means
  biologically (see Results §3 + Discussion).
- n_sim choice: 99 in validation (nb04), 199 in production v1 HPC run, 99
  in v2 HPC run. Empirically equivalent at the production window (A.6,
  2026-06-02 morning).

### 5. Significance: continuous SES instead of binary envelope (~400 words)
- Originally adopted a binary criterion: `L_obs > L_hi at ≥3 consecutive
  radii`. Yielded 1 of 438 pairs under v1 production (sparse-network problem).
- **Methodological finding (logged in decision_log 2026-06-02 morning):**
  under tight observation windows the binary envelope test is under-powered
  even on positive controls. Window choice is the dominant lever for
  inferential power. A larger window inflates contrast (includes empty
  space); a tight hull suppresses it (label-swap re-distributes labels over
  the same tightly-bounded set of locations the observed pattern already
  occupies).
- Solution: rank pairs by continuous **standardised effect size**:
  SES(r) = (L_obs(r) − midline(r)) / half_width(r). Peak SES = max over r.
- Defensible as the continuous form of the global rank envelope test
  (Myllymäki et al. 2017). Cite explicitly.

### 6. Panel construction (~300 words)
- v1 panel: 146 viable CellChatDB pairs (gene presence + abundance ≥ 50
  transcripts/gene/strip).
- v2 extended panel: 81 pre-registered pairs from
  `notes/extended_panel_rationale.md`, organised into 8 biological groups
  (G0 = v1 CellChatDB; G1 cytokeratin; G2 epi cell-type — wiped out by
  panel constraint; G3 epi niche; G4 immune; G5 stromal; G6 negative
  anchor; G7 infection-response; G8 paracrine).
- **Panel constraint to flag in Discussion:** the CosMx 1,000-plex is the
  immuno-oncology / inflammation panel, not a lung-anatomy panel. 21 of 81
  v2 pairs dropped because lung-resident epithelial markers (SCGB, MUC,
  SFTP, AGER, HOPX, FOXJ1, TP63) are absent. Method would test them if
  present.
- Pooled v1 + v2 = 614 (pair, strip) rows after gene/abundance filtering.
- HPC: 438 jobs v1 at n_sim = 199; 176 jobs v2 at n_sim = 99; KCL CREATE
  `msc_appbio` partition; SLURM array with concurrency cap 2; scripts in
  `scripts/`.

### 7. Network construction (~300 words)
- Pooled (pair, strip) rows ranked by peak SES; top-50 per strip as edges.
- NetworkX graph; edge weight = peak SES; nodes = genes.
- Louvain modularity-based community detection (single seed = 42; report
  stability across multiple seeds in supplementary).
- Leiden as cross-check (future work, `notes/future_work.md` Appendix A).

---

## Results (target ~3,500 words)

### 1. Pipeline validation on positive and negative controls (~500 words)
- nb04 / nb07 / Fig B.
- Under rectangular window + n_sim = 99 (nb04 baseline), KRT8 × KRT18
  cleanly exceeds envelope on all three strips. MALAT1 × KRT18 and
  KRT8 × SCGB3A1 correctly null. **Algorithm works at production scale.**
- Under production concave hull window + label-swap null, the binary
  envelope test is no longer met by the positive control. **But SES
  discriminates** — KRT8 × KRT18 peak SES = −2.18, MALAT1 × KRT18 peak
  SES = −4.85, separation = +2.67 envelope half-widths (Fig B, the side-
  by-side L(r) plot with the Option-2 caption from `poster_design.md`).
- This reframes the central inferential question (forward into §3).

### 2. Window choice as the principal inferential lever (~400 words)
- Window comparison from nb09c (rect vs hull); A.6 evidence.
- Tight hull preserves tissue geometry but suppresses absolute K(r)
  contrast because the label-swap null re-distributes across the same
  high-density locations.
- Recommend continuous-SES interpretation under tight windows.

### 3. Group-wise SES: method is correctly null on diffuse marginals (~600 words)
- Fig C.
- Pooled panel (v1 + v2, 614 (pair, strip) rows, 8 groups).
- Mean peak SES per group, table form. Quote from `decision_log.md`
  2026-06-03 morning:
  - G8 Paracrine +0.40
  - G7 Infection +0.35
  - G5 Stromal +0.31
  - G4 Immune +0.30
  - G0 CellChatDB +0.27
  - G6 Negative anchor −0.16
  - G1 Cytokeratin −0.19
  - G3 Epi niche −1.98
- Interpretation:
  - **Pan-tissue / diffuse marginal pairs (G1, G6) score null.** This is
    not a failure — it is the correct behaviour of bivariate K under a
    label-swap null (Methods §4). Cytokeratin co-expression is a *marginal*
    property; the bivariate K asks an *interaction* question.
  - **Focal / discrete-sub-region pairs (G4, G5, G7, G8) score positive.**
    These have spatial structure beyond marginals.
  - **G3 anti-association is a real signal.** KRT19 × SCGB3A1, EPCAM ×
    SCGB3A1 occupy disjoint epithelial sub-niches. SES at −4.10 is the
    largest single effect in the panel, just in the negative direction.
    Bivariate K detects *spatial avoidance*, not just association.
  - G0 CellChatDB scoring +0.27 (higher than G1) — incidental finding worth
    noting: the curated L-R database is enriched for biologically focal
    signalling pairs vs constitutive co-expression.

### 4. Infection biology (BIOLOGICAL HEADLINE) (~900 words)
- Fig F. Source for everything below: `decision_log.md` 2026-06-03 afternoon.
- Per-pair Strip-2 boost = SES_strip2 − mean(SES_strip1, SES_strip3) for
  the 13 surviving G7 pairs.
- Sub-theme breakdown:

  | Sub-theme | n | Positive | Mean boost |
  |---|---:|---:|---:|
  | Immune cell communication | 6 | 5 | +0.20 (or +0.39 excl. CXCL8 × IL6) |
  | Pure ISG–ISG | 4 | 3 | +0.16 |
  | Antiviral × epithelium | 3 | 1 (zero) | −0.50 |

- **Primary finding (Option A): immune cell-cell signalling elevates on
  infected tissue.** Top 4 by boost: HLA-DRA × CD3D (+0.74), CXCL10 × CD3D
  (+0.59), HLA-DRA × CD68 (+0.32), ISG15 × CD68 (+0.15). Textbook
  type I IFN-driven adaptive-immune recruitment: antigen presentation onto
  T cells, IFN-induced T-cell-recruitment chemokine signalling, antiviral
  activity in macrophages.
- **Paradoxical secondary finding (Option B): antiviral × epithelium pairs
  move OFF infected tissue.** ISG15 × KRT8 (−0.64), MX1 × KRT5 (−0.90).
  Two consistent biological readings: (a) RSV cytopathic effect — infected
  epithelium dies/sloughs, epithelial transcript density drops on infected
  strip, ISG–KRT joint pattern weakens; (b) ISG production migrates from
  epithelial cells to infiltrating immune cells — supported by ISG15 × CD68
  being on the positive side. **This is a finding only a segmentation-free
  method can make** (cell-based tools pre-bin transcripts before doing
  spatial analysis and cannot see redistribution across compartments).
  This is the unique-contribution argument.
- **CXCL8 × IL6 outlier interpretation (~150 words):** the one immune-
  communication negative pair. CXCL8 / IL6 is the NF-κB / classical-
  proinflammatory cytokine axis. RSV drives the type I IFN pathway
  (IRF3 → CXCL10, ISGs, MHC) more than NF-κB. Severe asthmatic lung
  carries a chronic CXCL8/IL6 baseline on the control strips. The outlier
  is consistent with "the acute viral response is IFN-driven, not
  NF-κB-driven" — textbook RSV phenotype.
- Statistical caveat: per-sub-theme sign tests don't reach formal
  significance individually (n = 3, 4, 6). **Run Mann-Whitney U on
  combined positive sub-themes (immune-comm + ISG–ISG, n = 10) vs
  antiviral × epithelium (n = 3)** — directional separation expected to
  be significant given the effect sizes. (TODO before submission.)

### 5. Network and community structure (~500 words) — Fig G done 2026-06-03
- Per-strip network from top-50 peak-SES edges; Louvain communities, only
  modules of ≥3 genes coloured (the rest are disjoint gene-pairs — a real
  sparsity property at top-50, stated honestly).
- **Result:** strip 2 (infected) carries more and larger connected modules
  (10 modules / 43 genes) than the controls (strip 1: 6 / 36; strip 3:
  8 / 38). A second, independent line of criterion-3 evidence beyond Fig F:
  infection elevates *connected* co-association structure, not just per-pair
  SES.
- Louvain + spring-layout both seed=42; seed stability across ≥5 seeds →
  supplementary (TODO).

### 6. L(r) profile structure recovers biological grouping (~300 words) — depends on Fig H
- UMAP of 50-dim L(r) profile vectors per (pair, strip).
- Expected: groups occupy distinct UMAP regions, demonstrating curve-level
  biological signal beyond peak SES.

### 7. Pipeline benchmark (Squidpy comparison) (~300 words) — currently TBD
- Deferred during poster sprint; ideally completed for dissertation Results.
- Reproduce KRT8 × KRT18 + 4 representative pairs with `sq.gr.ripley`,
  compare peak L(r), envelope width, significance calls.
- Honest reporting if results differ.

---

## Discussion (target ~2,500 words)

### 1. Comparison with the existing toolset (~500 words)
- The verified competitor landscape (Bento, Squidpy, CellChat/LIANA,
  Baysor, FICTURE, spatstat). All require either segmentation or do not
  do L-R; NoSegs fills a vacant slot.
- Quote the table from `decision_log.md` 2026-06-01 — verified positions.
- Explicit credit to Bento for subcellular work and to Squidpy for cell-
  level Ripley; this work complements rather than displaces them.

### 2. What the test actually measures (the central methods point) (~500 words)
- Bivariate K + label-swap null detects *excess spatial co-association
  beyond marginal density*. Not "do A and B co-occur in the same regions".
- This explains the cytokeratin "null" result (Methods §4, Results §3).
- This explains the antiviral × epithelium "negative" result (Results §4):
  observed joint pattern is *less* spatially adjacent than the null because
  the genes occupy disjoint sub-compartments in infected tissue.
- This is the conceptual reframing the dissertation contributes — the
  method does a more sophisticated job than naive "co-occurrence detector".

### 3. Biological readings (~500 words)
- Immune cell-cell signalling on infected tissue: antigen presentation,
  chemokine-mediated T-cell recruitment, monocyte/macrophage antiviral
  activity. The textbook adaptive-immune RSV response.
- Antiviral programme redistribution: ISG transcripts migrate from
  epithelium to immune compartments. Cytopathic effect contributes.
- Type I IFN vs NF-κB axis: RSV drives IFN preferentially. CXCL8 × IL6
  outlier supports this.
- Caveat: single sample, single platform, single donor; needs replication.

### 4. Limitations (~500 words)
- **Window choice as inferential lever** (Methods finding) — wider windows
  inflate contrast, tighter suppress. Recommend continuous SES not binary
  significance under tight windows.
- **Hull edge correction negative bias** on absolute K(r) — cancels in the
  permutation test, so only envelope-relative claims interpreted. Document
  explicitly.
- **Single Louvain seed** (production), single edge-weight choice (peak SES,
  not area-above-midline). Stability checks for supplementary.
- **Pointwise envelope inflation** of Type I error across multiple radii;
  preferable to switch to global rank envelope (Myllymäki 2017) for any
  per-pair significance claim in supplementary.
- **Panel constraint:** CosMx 1,000-plex is immuno-oncology, lung-resident
  epithelial markers absent (SCGB, MUC, SFTP, AGER, HOPX, FOXJ1, TP63).
  Method would test them if present.
- **Single sample (S1).** Replication across S2–S4 and external CosMx
  datasets is future work (`notes/future_work.md`).
- **Statistical power for per-sub-theme tests is low** (n = 3, 4, 6 per
  G7 sub-theme) — directional patterns are coherent and effect sizes are
  substantial, but formal sign tests don't reach α = 0.05 individually.

### 5. Future work (~500 words)
- Pull straight from `notes/future_work.md`:
  - Vectorise Shapely edge correction → rasterised mask (Appendix B).
  - Leiden over Louvain for guaranteed-connected communities (Appendix A).
  - Pseudo-segmentation via marker-gene DBSCAN clustering for per-cell-type
    stratification.
- External dataset validation (public CosMx if available; otherwise Xenium
  with similar tissue).
- Squidpy benchmark on shared pairs (if not in Results §7).
- Multi-sample replication across S2–S4.

---

## Abstract (target ~250 words)

Draft this last. Should cover:
- Problem: segmentation fails on this CosMx sample (42% dropout); standard
  L-R tools become unusable.
- Method: bivariate Ripley K + label-swap null on raw transcripts; SES
  ranking; panel-scale screen (614 pair-strip rows).
- Validation: positive controls discriminate from negatives; method is
  correctly null on diffuse marginals (cytokeratins) and elevated on focal
  sub-themes (immune, stromal, paracrine, infection).
- Biology: immune cell-cell signalling elevates on infected tissue;
  antiviral pairs redistribute off epithelium onto immune cells — a
  finding only segmentation-free analysis can make.
- Contribution: first segmentation-free L-R co-localisation pipeline using
  bivariate Ripley K on raw transcript point clouds, with hull-window edge
  correction and panel-scale application.

---

## Style notes
- Avoid the bold novelty claim "first segmentation-free spatial
  transcriptomics" (FICTURE exists; not true). Use the more specific
  defensible claim quoted in §Introduction.5 above.
- Avoid the binary-significance "passes / fails the envelope test" framing.
  Lead with continuous SES throughout Results and Discussion.
- Each Results subsection should end with one sentence connecting back to
  one of the three Aims (A1 / A2 / A3 from Introduction.5).
- Distinguish *statistical* from *biological* significance throughout.
  Honest about both. The directional coherence of the G7 sub-themes is
  the substance; the per-sub-theme p-values are weak due to small n.

---

## TODO checklist for the dissertation phase (W3–W6 of master plan)

- [ ] Run Mann-Whitney U: positive sub-themes (immune + ISG, n = 10) vs
      antiviral × epithelium (n = 3) on Strip-2 boost values.
- [ ] Run Louvain stability across ≥ 5 seeds; report Q distribution.
- [ ] Edge-weight sensitivity: peak SES vs area-above-midline (3 graphs,
      compare community memberships).
- [ ] Global rank envelope test (Myllymäki) as supplementary significance
      check on a small subset.
- [ ] Squidpy benchmark on 5 representative pairs.
- [ ] Generate UMAP (Fig H) and per-strip networks (Fig G) if not done
      during poster sprint.
- [ ] Format institutional template; word-count check per section.
- [ ] First full draft to Dan Nicolau by end of June 2026.
