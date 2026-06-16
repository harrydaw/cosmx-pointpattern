# NoSegs — Status

> Single living status board. Update at the end of each work session. Anything older than ~3 days here is stale.

**Last update:** 2026-06-16 — **Dissertation write-up prep.** Built the `Final_Writeup/` write-up system + `Final_Writeup/METHODS_LOG.md` (per-notebook reproducibility record) and ran a full verification pass against code/data. **Project renamed NoSegs → STIPPLE** (*Spatial, Transcript-level Inference of Pairwise Proximity in Ligand-receptor Expression*) — name only; on-disk folder, GitHub repo, and HPC `/scratch/.../NoSeggs/...` paths unchanged by design. **KRT8×SCGB3A1 re-run as a genuine 3-strip negative control** (concave hull, n_sim=199 → `results/controls/`, via `scripts/run_controls.sh`; supersedes the 2026-06-04 "dropped" decision). **Doc corrections (README + notebook_audit + decision_log patched):** canonical cleaner is **nb10** (not nb08b); 18 manual cluster IDs recorded; QC = 7.6% noise + 2.0% small + 5.4% manual → **87% kept of 702,873**; strip-assignment + DBSCAN run on PCA-rotated `x_rot_px`, but the **K-function runs on raw coords (rotation-invariant)**; CellChatDB via `liana`; geopandas/anndata/statsmodels/scanpy unused. Environments: HPC Python 3.11.6 (numpy 2.4.4, pandas 3.0.2…), local 3.13.11. **Next:** draft Methods §1 from METHODS_LOG. *(Poster status below is historical — poster delivered for Rome.)*

---

**Last update:** 2026-06-04 (finishing pass) — Poster print-ready pass done. **(1) M6 flowchart rebuilt as native, individually-editable PowerPoint shapes** (drawn in `docs/_build_poster.py` `flowchart()`, no longer a flat PNG): A1-punchy `roundRect` column, **teal `#1A7A7A` marks the two statistical-core stages** (bivariate K→L, permutation envelope), purple prep/output stages, amber option-chips on the tunable steps, green-tick quick-start panel. **(2) "NoSegs" spelling canonical everywhere** incl. both `.pptx` filenames (`git mv` → `NoSegs_Poster_A0.pptx`, `NoSegs_Poster.pptx`); the 3 HPC `scripts/*.sh` deliberately keep "NoSeggs" in their `/scratch/.../NoSeggs/...` paths — renaming would break CREATE re-runs (flag if you ever rename the remote dir). **(3) Results take-home sharpened** — purple banner over Fig F leads with a bold headline + big **5/6** callout + "+2.67 envelope half-widths"; one `[ FILL: confirm headline wording ]` left. **(4) Negative-2 control row dropped** (KRT8×SCGB3A1 was never run — absent from both parquets; was a fabricated placeholder). nb14 re-executed (kernel `noseggs`) so `14_M0_raw_data.png` title reads "NoSegs"; poster rebuilt + COM-rendered, verified clean. **New write-up kit:** `notes/poster_writeup_guide.md` — per-section reading list + key points + what to write + the one reference gap (RSV/lung infection biology, ~5 refs to gather). **Remaining for Harry:** confirm Fig F headline wording, pick 42% vs 43.6% dropout (use consistently), fill Acknowledgements funder/PI, verify QR URL, spacing touch-ups, A4 test print → A0 order. Prior (late): added two nb14 figures and folded them into the A0 layout. **`14_M0_raw_data.png`** (one FOV, 3 panels: DAPI → vendor seg mask → raw transcripts assigned-vs-unassigned, ~40% magenta) leads "The Problem"; **`14_M6_workflow.png`** (8-stage flowchart with amber option-chips + quick-start panel) is the Methods visual anchor. Both 300 dpi, Arial, executed in nb14 (exec 1–9). **nb13 re-executed → B/C/F/G/H regenerated in Arial** (Fig B dpi 300), so the whole poster is now font-consistent. Poster `docs/NoSegs_Poster_A0.pptx` rebuilt with all Arial figs + M0 + M6, COM-rendered, fits clean across all 3 columns — no overflow. **Infra:** notebook execution was failing on a stale global kernelspec — fixed by installing a venv-prefix `noseggs` kernel (ran via `python -m nbconvert --execute --ExecutePreprocessor.kernel_name=noseggs`; the `jupyter` dispatcher swallows output, use `-m nbconvert`). **Remaining for Harry:** fill the `[ FILL ]` prose prompts + fact-checks, touch up spacing in PowerPoint, then A4 test print → A0 order. Prior (evening): A0 poster first assembled (true A0 84.1×118.9 cm; A1 `NoSegs_Poster.pptx` kept as fallback; native editable text/table/pictures, `[ FILL ]` prompts + legends).

**Supervisor:** Dr Dan Nicolau (Nicolau Lab, KCL).

---

## Deliverables in flight

| # | Deliverable | Deadline | State |
|---|---|---|---|
| 1 | NeuroMonster A0 poster — **PRINT** | **Sat 6 June 2026** (fly Rome Sun 7) | All figures ready. Results: Fig F (infection) locked; C/G/H done; D/E banked. **Methodology suite (nb14):** M1 dropout, M2 GMM strips, M3 DBSCAN QC, M4 window iterations, M5 r-scale + T 2D positioning map — all 300 dpi, poster-grade. pptx layout pass not started — use REVISED methodology-forward layout in poster_design.md |
| 2 | MSc dissertation draft → Dan Nicolau | End of June 2026 | 0 / 12,000 words |
| 3 | MSc dissertation submission | 16 July 2026 | — |
| 4 | Full-length presentation | Built in Week 5 (post-supervisor feedback) | — |

---

## What's done

**Pipeline + data (since project start):**
- Full pipeline coded across 19 notebooks (00 → 13).
- Positive controls (KRT8 × KRT18) and negative controls validated at n_sim=99.
- 438-job HPC array on KCL CREATE at n_sim=199 — all jobs OK. Aggregated to `results/lr_panel_results.parquet`.
- DBSCAN noise removal + manual cluster QC.
- Concave hull windows + Shapely edge correction.
- Vendor segmentation-failure figure: `results/figures/01_vendor_segmentation_dropout.png`.
- Strong writeup scaffolding: `notes/intro_methods_reading.md` (120+ refs + key terms), `notes/dissertation_outline.md` (section spine + word budgets), `notes/poster_writeup_guide.md` (per-section reading + key points), `notes/decision_log.md` (per-decision rationale). *(The old `planning_docs/04_background_topics.md` + `05_methods_decisions.md` no longer exist — their content folded into these notes.)*
- Reboot streamline (2026-06-01): repo cleaned, `02_work_to_do.md` archived, `future_work.md` merged, nb09b + old pptx archived, empty `src/` and `tests/` removed, README rewritten with verified novelty positioning.

**Wed 3 June work session:**
- v2 extended panel (81 hand-picked pairs from `extended_panel_rationale.md`) run completed overnight on HPC at n_sim=99. 176 of 243 jobs survived gene + abundance filters (panel constraint: lung-resident epithelial markers absent from CosMx 1000-plex).
- v2 aggregated to `results/lr_panel_results_v2.parquet` (176 rows).
- nb13: pooled v1 + v2 (614 rows), attached pre-registered Group labels (G0 CellChatDB, G1 cytokeratin, G3 epi niche, G4 immune, G5 stromal, G6 negative anchor, G7 infection, G8 paracrine), computed peak SES and r_at_peak per row.
- **Critical interpretive finding:** the bivariate K under label-swap null measures *excess spatial co-association beyond marginal density*, not "do A and B co-occur in the same regions". Pan-tissue markers (cytokeratins) correctly score null; focal spatially-structured pairs score positive; pairs in disjoint sub-niches can score *negative* (anti-association). This is the central methodological framing for the dissertation.
- Fig C generated: group-wise SES violin — clear evidence the method is correctly null on diffuse marginals (G1, G6) and elevated on focal sub-themes (G4, G5, G7, G8).
- Fig D generated then BANKED: SES heatmap (203 pairs × 3 strips). Competent landscape but not poster-grade visual punch.
- Fig E generated then BANKED: volcano (peak SES vs r_at_peak). Discrete r-banding doesn't look polished. Dissertation only.
- **Fig F locked (BIOLOGICAL CENTREPIECE):** per-pair Strip-2 boost for G7 infection pairs, sub-themed. Immune cell communication (5/6 positive), Pure ISG–ISG (3/4 positive), Antiviral × epithelium (1/3 positive, mean −0.50). Framing decision logged in `notes/decision_log.md` Wed PM entry.

## What's next (3 bullets)

1. **Thu 4 June — poster layout: DONE (draft).** `docs/NoSegs_Poster_A0.pptx` built at A0 with the methodology-forward layout (Left = Intro + M1 → M2 → M3; Middle = pipeline + M4 + M5 + T + Fig B + validation table + Fig C; Right = Fig F + Fig G + Fig H; bottom band = conclusions/refs/acks/QR). Next on this: **Harry to fill the `[ FILL ]` prompts** (intro fact-check, results headline, Negative-2 SES, references, funder attribution) and touch up spacing in PowerPoint. Optional: re-run nb13 results figs with Arial pinned for full font consistency.
2. **Fri 5 June — A4 test print** (legibility check at 8× shrink), then A0 order at KCL Print Services / local shop (confirm 24h turnaround).
3. **Sat 6 June — collect A0**, tube/case for travel; fly Rome Sun 7.

*(Done 2026-06-03 evening: Fig G per-strip Louvain network + Fig H L(r)-profile UMAP, both in nb13 §8/§9 → `results/figures/13_G_network.png`, `13_H_umap.png`.)*

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

**Original decision (2026-06-02 PM):** defer until after Rome.
**Revised decision (2026-06-02 evening):** ACCELERATED — submitted Tuesday night to remove queue risk from Wednesday. Completed overnight; aggregated Wed morning. Result feeds both poster (Fig F sub-themes) and dissertation Results chapter.

## Open decisions / blockers awaiting input

- TBD: whether to defer Squidpy benchmark to post-poster (Week 2 onwards). Currently leaning toward defer — Saturday print deadline doesn't leave time.

## Live framing (Moderate–Bold, NeuroMonster-tuned)

**Headline:** Segmentation-free ligand–receptor co-localisation in spatial transcriptomics via bivariate point-pattern statistics.

**30-second pitch (post-Fig F, locked 2026-06-03 PM):** "Most tools for finding cell communication signals in spatial transcriptomics first segment the data into cells — but segmentation fails on noisy clinical samples, dropping 42% of transcripts in ours. I built a pipeline that skips segmentation entirely: each gene becomes a 2D point pattern, and bivariate Ripley's K detects excess spatial co-association beyond what marginal density predicts. Applied to RSV-infected lung, immune cell-communication pairs preferentially co-localise on infected tissue — antigen presentation, T-cell recruitment, antiviral monocyte signalling. And a paradoxical second finding only segmentation-free analysis can surface: antiviral genes that *don't* sit near epithelium on infected tissue, consistent with epithelial loss redistributing antiviral activity into immune cells."

**Do NOT claim** (verified against literature 2026-06-01): "first to use Ripley's K in spatial transcriptomics" (Squidpy has it); "segmentation-free spatial transcriptomics" as a category (FICTURE exists). DO claim: "first segmentation-free L-R co-localisation pipeline using bivariate Ripley's K on transcript point clouds with hull-window edge correction and panel-scale application." See `~/.claude/plans/hello-mr-claude-i-polished-manatee.md` for full landscape.

## Operational plan reference

Current 6-week plan: [`~/.claude/plans/hello-mr-claude-i-polished-manatee.md`](file:///C:/Users/Harry/.claude/plans/hello-mr-claude-i-polished-manatee.md)

Supersedes the now-archived `planning_docs/02_work_to_do.md` (whose 15 May / 31 May milestones slipped).

## Decision history

Live chronological log of methodological decisions and scientific findings:
[`notes/decision_log.md`](decision_log.md). This is what feeds the dissertation
Methods rationale and Discussion limitations sections — keep it updated whenever
we make a methodological call.

## Dissertation outline

Section-by-section outline with key points + figure references + source-doc
pointers: [`notes/dissertation_outline.md`](dissertation_outline.md). Use this
as the *starting point* for the writeup. Pulls together findings from
`decision_log.md`, captions from `poster_design.md`, reading from
`intro_methods_reading.md`, decisions from `planning_docs/05_methods_decisions.md`.

## Handoff context for a new chat session (2026-06-03 evening)

If a new chat picks this up cold, the four files to read in order are:
1. `notes/STATUS.md` (this file) — current state, what's next, framing.
2. `notes/decision_log.md` — every methodological + scientific finding to date,
   newest first.
3. `notes/poster_design.md` — figure list, locked captions, layout plan.
4. `notes/dissertation_outline.md` — bullet outline of the 12,000-word writeup
   with cross-references to the source docs.

After those, read `README.md` for the codebase orientation and
`notes/notebook_audit.md` for per-notebook state.
