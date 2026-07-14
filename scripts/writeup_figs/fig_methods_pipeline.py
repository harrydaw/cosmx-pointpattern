"""Methods fig:pipeline (FIRST DRAFT — Harry to art-direct). Vertical flowchart of the STIPPLE
pipeline, light/IBM, grouped by phase. Pure diagram (no data)."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path
FIG = Path(__file__).resolve().parents[2] / "results" / "figures"
plt.rcParams.update({"font.family": "Arial", "figure.facecolor": "white", "axes.facecolor": "white",
                     "savefig.facecolor": "white"})
IBM = {"blue": "#648FFF", "purple": "#785EF0", "magenta": "#DC267F", "orange": "#FE6100"}

# (label, phase)
steps = [
    ("CosMx transcripts\nS1: 12 FOVs, 702,873 transcripts", "data"),
    ("PCA rotation + per-FOV GMM\n→ three tissue strips", "data"),
    ("DBSCAN + manual QC\n611,150 retained (87%)", "data"),
    ("Concave-hull window\n+ edge correction", "stat"),
    ("Bivariate Ripley's $K$ → $L(r)$\n50 radii, 5–250 px", "stat"),
    ("Label-permutation envelope\n$n_{sim}=199$, seed 42", "infer"),
    ("SES$(r)$; rank by Peak SES", "infer"),
    ("Per-strip co-localisation networks", "out"),
]
PC = {"data": IBM["blue"], "stat": IBM["purple"], "infer": IBM["magenta"], "out": IBM["orange"]}
PH = {"data": "Data preparation", "stat": "Spatial statistic", "infer": "Inference", "out": "Output"}

fig, ax = plt.subplots(figsize=(6.6, 9.2)); ax.set_xlim(0, 10); ax.set_ylim(0, len(steps)*1.32 + 0.4); ax.axis("off")
n = len(steps); y0 = len(steps)*1.32 - 0.2; dy = 1.32; bx, bw, bh = 3.0, 5.2, 0.92
cx = bx + bw/2
centers = []
for i, (lab, ph) in enumerate(steps):
    y = y0 - i*dy
    col = PC[ph]
    ax.add_patch(FancyBboxPatch((bx, y-bh/2), bw, bh, boxstyle="round,pad=0.02,rounding_size=0.12",
                                linewidth=1.6, edgecolor=col, facecolor=col, alpha=0.16, zorder=2))
    ax.add_patch(FancyBboxPatch((bx, y-bh/2), bw, bh, boxstyle="round,pad=0.02,rounding_size=0.12",
                                linewidth=1.6, edgecolor=col, facecolor="none", zorder=3))
    ax.text(cx, y, lab, ha="center", va="center", fontsize=9.3, zorder=4)
    centers.append((y, ph))
# arrows
for i in range(n-1):
    y_top = centers[i][0] - bh/2; y_bot = centers[i+1][0] + bh/2
    ax.add_patch(FancyArrowPatch((cx, y_top), (cx, y_bot), arrowstyle="-|>", mutation_scale=14,
                                 lw=1.4, color="0.35", zorder=1))
# phase brackets/labels on the left
seen = []
for i, (y, ph) in enumerate(centers):
    if ph not in seen:
        seen.append(ph)
        ys = [c[0] for c in centers if c[1] == ph]
        ytop, ybot = max(ys)+bh/2, min(ys)-bh/2
        ax.plot([bx-0.5, bx-0.5], [ybot, ytop], color=PC[ph], lw=3, solid_capstyle="round", zorder=2)
        ax.text(bx-0.75, (ytop+ybot)/2, PH[ph], rotation=90, ha="center", va="center",
                fontsize=9, color=PC[ph], fontweight="bold")
fig.savefig(FIG / "pub_pipeline.png", dpi=300, bbox_inches="tight"); print("saved", FIG / "pub_pipeline.png")
