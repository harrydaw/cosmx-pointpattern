# Decision log

> Chronological record of methodological decisions, scientific findings, and
> their justifications. Companion to `planning_docs/05_methods_decisions.md`
> (which is the structured per-decision rationale) — this file is the
> *sequential narrative* of how we got here. Feeds dissertation Methods §
> rationale and Discussion § limitations.
>
> Newest entries at the top.

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
