# NoSeggs — Status

> Single living status board. Update at the end of each work session. Anything older than ~3 days here is stale.

**Last update:** 2026-06-01 (reboot day)

---

## Deliverables in flight

| # | Deliverable | Deadline | State |
|---|---|---|---|
| 1 | NeuroMonster A0 poster | ~early June (~1 week) | Draft `docs/NoSeggs_Poster.pptx` from 8 May; needs final figures + new framing |
| 2 | MSc dissertation draft → supervisor | 1 July 2026 | 0 / 12,000 words |
| 3 | MSc dissertation submission | 16 July 2026 | — |
| 4 | Full-length presentation | TBD (built off dissertation in Week 5) | — |

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

1. **nb13 Step A.6** — KRT8 × KRT18 strip 2 at n_sim=199 vs n_sim=99 side-by-side. Saved figure to `results/figures/13_A6_KRT_nsim_comparison.png`.
2. **nb13 Steps A.3, A.5, A.7** — n_consec sweep, effect-size network, r_at_peak by strip.
3. **Decision gate Thursday 4 June** — C.1 / C.2 / C.3 (see `network_revival_plan.md` Step C). Log decision below.

## C decision (network revival)

*To be filled in on or before Thursday 4 June 2026.*

- [ ] C.1 — KRT survives at n_sim=199 → effect-size pivot, no HPC re-run
- [ ] C.2 — KRT fails at 199, passes at 99 → schedule HPC re-run at n_sim=99 with extended panel
- [ ] C.3 — no signal under either → methodological-contribution framing

## Open decisions / blockers awaiting input

- NeuroMonster exact poster deadline + abstract (if already submitted).
- Supervisor name + feedback cadence (single round 1 July vs rolling).
- Current contents of `docs/NoSeggs_Poster.pptx` — describe or screenshot for targeted critique.

## Live framing (Moderate–Bold, NeuroMonster-tuned)

**Headline:** Segmentation-free ligand–receptor co-localisation in spatial transcriptomics via bivariate point-pattern statistics.

**30-second pitch:** "Most tools for finding cell communication signals in spatial transcriptomics first segment the data into cells — but segmentation fails on noisy clinical samples, dropping up to 42% of measurements in ours. I built a pipeline that skips segmentation entirely: it treats each gene as a 2D point pattern and uses Ripley's K-function with permutation tests to detect co-localised ligand–receptor pairs directly. Validated on positive controls and applied to RSV-infected human lung."

**Do NOT claim** (verified against literature 2026-06-01): "first to use Ripley's K in spatial transcriptomics" (Squidpy has it); "segmentation-free spatial transcriptomics" as a category (FICTURE exists). DO claim: "first segmentation-free L-R co-localisation pipeline using bivariate Ripley's K on transcript point clouds with hull-window edge correction and panel-scale application." See `~/.claude/plans/hello-mr-claude-i-polished-manatee.md` for full landscape.

## Operational plan reference

Current 6-week plan: [`~/.claude/plans/hello-mr-claude-i-polished-manatee.md`](file:///C:/Users/Harry/.claude/plans/hello-mr-claude-i-polished-manatee.md)

Supersedes the now-archived `planning_docs/02_work_to_do.md` (whose 15 May / 31 May milestones slipped).
