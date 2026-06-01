# NoSeggs — Status

> Single living status board. Update at the end of each work session. Anything older than ~3 days here is stale.

**Last update:** 2026-06-02 (after Step A.6)

**Supervisor:** Dr Dan Nicolau (Nicolau Lab, KCL).

---

## Deliverables in flight

| # | Deliverable | Deadline | State |
|---|---|---|---|
| 1 | NeuroMonster A0 poster — **PRINT** | **Sat 6 June 2026** (fly Rome Sun 7) | Draft layout in place; text is filler; figures placeholder; no abstract submitted, framing open |
| 2 | MSc dissertation draft → Dan Nicolau | End of June 2026 | 0 / 12,000 words |
| 3 | MSc dissertation submission | 16 July 2026 | — |
| 4 | Full-length presentation | Built in Week 5 (post-supervisor feedback) | — |

---

## What's done

- Full pipeline coded across 17 notebooks (00 → 13).
- Positive controls (KRT8 × KRT18) and negative controls validated at n_sim=99.
- 438-job HPC array on KCL CREATE at n_sim=199 — all jobs OK. Aggregated to `results/lr_panel_results.parquet`.
- DBSCAN noise removal + manual cluster QC.
- Concave hull windows + Shapely edge correction.
- Vendor segmentation-failure figure: `results/figures/01_vendor_segmentation_dropout.png`.
- Strong writeup scaffolding: `planning_docs/04_background_topics.md`, `notes/intro_methods_reading.md` (120+ refs), `planning_docs/05_methods_decisions.md` (15 documented decisions).
- Reboot streamline (2026-06-01): repo cleaned, `02_work_to_do.md` archived, `future_work.md` merged, nb09b + old pptx archived, empty `src/` and `tests/` removed, README rewritten with verified novelty positioning.

## What's next (3 bullets)

1. **nb13 Step A.5 (Wed AM):** rank all 438 (pair, strip) rows by peak SES; build effect-size network from top-N per strip; run Louvain → headline poster figure.
2. **nb13 Step A.7 (Wed PM):** r_at_peak histogram by strip for top-N SES pairs.
3. **Poster sprint (Thu–Sat):** see `notes/poster_design.md` for the locked decisions and column-by-column content plan.

## C decision (network revival) — DECIDED 2026-06-02

**Outcome: C.3 with effect-size pivot.**

A.6 result (production hull window, KRT8×KRT18 strip 2):
- n_sim=99: peak excess = −1.56 px, max consecutive exceedances = 0 → FAILS
- n_sim=199: peak excess = −1.56 px, max consecutive exceedances = 0 → FAILS
- The two envelopes are nearly identical. **n_sim is not the lever.**

A.6 sanity check (pos vs neg, same production window):
- KRT8 × KRT18 (positive) peak SES = −2.18
- MALAT1 × KRT18 (negative) peak SES = −4.85
- **Separation = +2.67 envelope half-widths.** The test discriminates by effect size, even though neither pair passes the binary envelope test.

Interpretation: the concave-hull window tightly captures tissue geometry, leaving the permutation null almost as concentrated as the observed pattern. This suppresses absolute contrast and makes the binary envelope test under-powered. The **window choice is the principal lever for inferential power** under this method — a genuine methodological finding for the dissertation.

Pivot: **adopt continuous SES ranking instead of binary significance.** Top-N by peak SES becomes the edge set for the network. Defensible because it is the continuous form of the global rank envelope test (Myllymäki 2017).

## Extended panel decision

**Decision (2026-06-02):** Option B — defer until after Rome.

The extended-panel run remains valuable, but for dynamic-range reasons (anchoring the SES ranking with known cytokeratin positives + housekeeping negatives) rather than the original n_sim=99 reason. Submitting before the poster risks Wednesday HPC queue delays blowing Thursday's layout day. Schedule for after 7 June; result feeds the dissertation Results chapter, not the poster.

## Open decisions / blockers awaiting input

- (Resolved 2026-06-01: poster prints Sat 6 Jun; supervisor Dan Nicolau; no abstract submitted; poster draft reviewed.)
- TBD: whether to defer Squidpy benchmark to post-poster (Week 2 onwards) given compressed schedule.

## Live framing (Moderate–Bold, NeuroMonster-tuned)

**Headline:** Segmentation-free ligand–receptor co-localisation in spatial transcriptomics via bivariate point-pattern statistics.

**30-second pitch (post-A.6):** "Most tools for finding cell communication signals in spatial transcriptomics first segment the data into cells — but segmentation fails on noisy clinical samples, dropping 42% of measurements in ours. I built a pipeline that skips segmentation entirely: each gene becomes a 2D point pattern, and bivariate Ripley's K with permutation tests detects co-localised ligand–receptor pairs directly. The test discriminates positive from negative controls by effect size, and a panel-scale screen of 146 pairs is in progress."

**Do NOT claim** (verified against literature 2026-06-01): "first to use Ripley's K in spatial transcriptomics" (Squidpy has it); "segmentation-free spatial transcriptomics" as a category (FICTURE exists). DO claim: "first segmentation-free L-R co-localisation pipeline using bivariate Ripley's K on transcript point clouds with hull-window edge correction and panel-scale application." See `~/.claude/plans/hello-mr-claude-i-polished-manatee.md` for full landscape.

## Operational plan reference

Current 6-week plan: [`~/.claude/plans/hello-mr-claude-i-polished-manatee.md`](file:///C:/Users/Harry/.claude/plans/hello-mr-claude-i-polished-manatee.md)

Supersedes the now-archived `planning_docs/02_work_to_do.md` (whose 15 May / 31 May milestones slipped).

## Decision history

Live chronological log of methodological decisions and scientific findings:
[`notes/decision_log.md`](decision_log.md). This is what feeds the dissertation
Methods rationale and Discussion limitations sections — keep it updated whenever
we make a methodological call.
