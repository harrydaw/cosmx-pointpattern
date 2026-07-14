"""Methods fig:mechanics — the estimator's moving parts, one dedicated composite:
A bivariate Ripley's K schematic (count gene-B within r of each gene-A) |
B the two-part radius grid (dense <=50 px, coarse to 250 px) |
C edge correction on the real concave-hull window (disc fraction inside). Light/IBM/Arial."""
import warnings; warnings.filterwarnings("ignore")
import json, io
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Circle, FancyBboxPatch, Polygon as MplPoly
from shapely.geometry import Point

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "results" / "figures"
plt.rcParams.update({"font.family": "Arial", "figure.facecolor": "white", "axes.facecolor": "white",
                     "savefig.facecolor": "white"})
IBM = {"blue": "#648FFF", "purple": "#785EF0", "magenta": "#DC267F", "orange": "#FE6100", "yellow": "#FFB000"}
INK, MUTED = "0.15", "0.45"
def panel(ax, L): ax.text(-0.02, 1.03, L, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom", ha="right")

# nb00 functions for the real edge correction
ns = globals()
nb = json.load(io.open(ROOT / "notebooks" / "00_functions.ipynb", encoding="utf-8"))
for c in nb["cells"]:
    if c["cell_type"] == "code":
        try: exec("".join(c["source"]), ns)
        except Exception: pass

R = np.load(ROOT / "results" / "lr_panel_results.parquet.r_vals.npy")

fig = plt.figure(figsize=(13.5, 4.7))
gs = GridSpec(1, 3, wspace=0.16, left=0.02, right=0.98, top=0.86, bottom=0.06)

# ===== A: bivariate Ripley's K schematic (synthetic, clean) =====
axA = fig.add_subplot(gs[0, 0]); axA.set_aspect("equal"); axA.axis("off")
rng = np.random.default_rng(3)
W0, W1 = 0.5, 9.5
Apts = rng.uniform(1.2, 8.8, (18, 2))
Bpts = rng.uniform(1.2, 8.8, (24, 2))
anchors = Apts[[3, 9, 14]]
Bpts = np.vstack([Bpts, np.repeat(anchors, 6, 0) + rng.normal(0, 0.5, (18, 2))])
rr = 1.5
axA.add_patch(FancyBboxPatch((W0, W0), W1 - W0, W1 - W0, boxstyle="round,pad=0.02,rounding_size=0.3",
                             fill=False, ec="0.55", lw=1.4, ls=(0, (4, 3))))
for fx, fy in anchors:
    axA.add_patch(Circle((fx, fy), rr, fill=True, fc=IBM["blue"], ec="none", alpha=0.10, zorder=1))
    axA.add_patch(Circle((fx, fy), rr, fill=False, ec=IBM["blue"], lw=1.6, zorder=2))
ins = np.zeros(len(Bpts), bool)
for fx, fy in anchors: ins |= np.hypot(Bpts[:, 0] - fx, Bpts[:, 1] - fy) <= rr
axA.scatter(Bpts[~ins, 0], Bpts[~ins, 1], s=45, c=IBM["magenta"], alpha=0.4, lw=0, zorder=3)
axA.scatter(Bpts[ins, 0], Bpts[ins, 1], s=62, c=IBM["magenta"], edgecolors="white", lw=0.8, zorder=4)
axA.scatter(Apts[:, 0], Apts[:, 1], s=42, c=IBM["blue"], alpha=0.5, lw=0, zorder=3)
axA.scatter(anchors[:, 0], anchors[:, 1], s=95, c=IBM["blue"], edgecolors="white", lw=1.2, zorder=5)
fx, fy = anchors[0]
axA.annotate("", xy=(fx + rr, fy), xytext=(fx, fy), arrowprops=dict(arrowstyle="-", color=INK, lw=1.5), zorder=6)
axA.text(fx + rr / 2, fy + 0.32, "r", color=INK, fontsize=15, ha="center", style="italic", zorder=6)
axA.scatter([], [], s=60, c=IBM["blue"], label="Gene A (focal)")
axA.scatter([], [], s=60, c=IBM["magenta"], label="Gene B (counted within r)")
axA.legend(frameon=False, fontsize=8.5, loc="lower center", bbox_to_anchor=(0.5, -0.08), ncol=1, handletextpad=0.3)
axA.set_xlim(-1.0, 10.5); axA.set_ylim(-0.6, 10.5); panel(axA, "A")
axA.set_title("Bivariate Ripley's $K$: count gene B\nwithin $r$ of each gene A", fontsize=10)

# ===== B: the two-part radius grid =====
axB = fig.add_subplot(gs[0, 1]); axB.set_aspect("equal"); axB.axis("off")
for rv in R:
    axB.add_patch(Circle((0, 0), rv, fill=False, ec="0.80", lw=0.5, zorder=1))
axB.add_patch(Circle((0, 0), 50, fill=False, ec=IBM["orange"], lw=2.0, zorder=3))
axB.add_patch(Circle((0, 0), 250, fill=False, ec=INK, lw=1.4, ls="--", zorder=3))
axB.scatter([0], [0], s=70, c=IBM["blue"], edgecolors="white", lw=0.8, zorder=5)
axB.annotate("50 px cut-off\n(30 dense radii within)", (0, 50), (95, 120), fontsize=8.5, color=IBM["orange"],
             ha="left", arrowprops=dict(arrowstyle="->", color=IBM["orange"], lw=1.1))
axB.annotate("250 px max\n(20 coarse radii, 10 px)", (177, 177), (150, 235), fontsize=8.5, color=INK,
             ha="left", arrowprops=dict(arrowstyle="->", color=INK, lw=1.1))
axB.set_xlim(-265, 265); axB.set_ylim(-265, 300); panel(axB, "B")
axB.set_title("Radius grid: 30 dense radii to 50 px,\n20 coarse radii to 250 px", fontsize=10)

# ===== C: edge correction on the real window =====
axC = fig.add_subplot(gs[0, 2]); axC.set_aspect("equal"); axC.axis("off")
XC, YC = "x_global_px_transformed", "y_global_px_transformed"
k = pd.read_parquet(ROOT / "data" / "processed" / "s1_all_strips_cleaned.parquet")
keep = (~k.is_noise) & (~k.is_small_cluster.fillna(False)) & (~k.get("manually_excluded", pd.Series(False, index=k.index)).fillna(False))
sdf = k[keep & (k.strip == "strip_2")]
allxy = sdf[[XC, YC]].values.astype(float)
A = sdf[sdf.target == "SPP1"][[XC, YC]].values.astype(float)
rC = 80.0
hull = get_concave_hull(sdf, x_col=XC, y_col=YC)
hx, hy = np.array(hull.exterior.xy[0]), np.array(hull.exterior.xy[1])
fr = np.array([fraction_inside_hull(Point(p), rC, hull) for p in A])
cand = np.where((fr > 0.45) & (fr < 0.75))[0]
fp = A[cand[np.argmin([hull.exterior.distance(Point(A[i])) for i in cand])]] if len(cand) else A[int(np.argmin(fr))]
frac = fraction_inside_hull(Point(fp), rC, hull)
inside = Point(fp).buffer(rC, resolution=64).intersection(hull)
axC.fill(hx, hy, facecolor="0.96", edgecolor="0.4", lw=1.2, zorder=0)
pad = rC * 1.9
allC = allxy[(np.abs(allxy[:, 0] - fp[0]) < pad) & (np.abs(allxy[:, 1] - fp[1]) < pad)]
axC.scatter(*allC.T, s=5, c="0.80", alpha=0.55, lw=0, zorder=1)
if not inside.is_empty and inside.geom_type == "Polygon":
    ix, iy = inside.exterior.xy
    axC.add_patch(MplPoly(np.column_stack([ix, iy]), closed=True, facecolor=IBM["blue"], alpha=0.35, lw=0, zorder=2))
axC.add_patch(Circle(fp, rC, fill=False, ec=INK, lw=1.6, ls="--", zorder=3))
axC.scatter(*fp, s=75, c=IBM["blue"], edgecolors="white", lw=0.9, zorder=4)
axC.set_xlim(fp[0] - pad, fp[0] + pad); axC.set_ylim(fp[1] - pad, fp[1] + pad); panel(axC, "C")
axC.set_title("Edge correction  $w=1/\\mathrm{frac}$\n(disc %.0f%% inside the window)" % (frac * 100), fontsize=10)

fig.savefig(FIG / "pub_mechanics.png", dpi=300, bbox_inches="tight")
print("saved", FIG / "pub_mechanics.png", "| edge-correction disc frac %.0f%%" % (frac * 100))
