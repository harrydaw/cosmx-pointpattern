"""Methods fig:window — the observation window matters. One strip under three
candidate windows, same extent throughout:
A axis-aligned bounding box | B convex hull | C concave hull (ratio 0.1, used).
The shaded area with no points under it is empty tissue the window wrongly counts
as sampled; the concave hull excludes almost all of it. Light/IBM/Arial."""
import warnings; warnings.filterwarnings("ignore")
import json, io
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle, Polygon as MplPoly
from scipy.spatial import ConvexHull

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "results" / "figures"
plt.rcParams.update({"font.family": "Arial", "figure.facecolor": "white", "axes.facecolor": "white",
                     "savefig.facecolor": "white"})
IBM = {"blue": "#648FFF", "purple": "#785EF0", "magenta": "#DC267F", "orange": "#FE6100", "yellow": "#FFB000"}
INK, MUTED = "0.15", "0.45"
def panel(ax, L): ax.text(-0.02, 1.03, L, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom", ha="right")

# nb00 functions for the real concave hull
ns = globals()
nb = json.load(io.open(ROOT / "notebooks" / "00_functions.ipynb", encoding="utf-8"))
for c in nb["cells"]:
    if c["cell_type"] == "code":
        try: exec("".join(c["source"]), ns)
        except Exception: pass

XC, YC = "x_global_px_transformed", "y_global_px_transformed"
k = pd.read_parquet(ROOT / "data" / "processed" / "s1_all_strips_cleaned.parquet")
keep = (~k.is_noise) & (~k.is_small_cluster.fillna(False)) & \
       (~k.get("manually_excluded", pd.Series(False, index=k.index)).fillna(False))
sdf = k[keep & (k.strip == "strip_2")]
xy = sdf[[XC, YC]].values.astype(float)

# three candidate windows
x0, y0, x1, y1 = xy[:, 0].min(), xy[:, 1].min(), xy[:, 0].max(), xy[:, 1].max()
box_area = (x1 - x0) * (y1 - y0)
ch = ConvexHull(xy); cvx = xy[ch.vertices]
cvx_area = ch.volume  # 2-D ConvexHull.volume is the polygon area
hull = get_concave_hull(sdf, x_col=XC, y_col=YC)
hx, hy = np.array(hull.exterior.xy[0]), np.array(hull.exterior.xy[1])
cav_area = hull.area

# the strip is a long thin curved sliver: rotate to horizontal (swap axes) and
# stack the three windows so each panel is wide-short and fills the space.
rng = np.random.default_rng(0)
sub = xy[rng.choice(len(xy), size=min(35000, len(xy)), replace=False)]
pad = 0.04 * max(x1 - x0, y1 - y0)
mpx = 1e6
def sw(a): return np.column_stack([a[:, 1], a[:, 0]])   # (x,y) -> (y,x) display
sub_d, cvx_d, hull_d = sw(sub), sw(cvx), np.column_stack([hy, hx])
HX, HY = (y0 - pad, y1 + pad), (x0 - pad, x1 + pad)      # display limits

def base(ax, title):
    ax.set_aspect("equal"); ax.axis("off")
    ax.scatter(sub_d[:, 0], sub_d[:, 1], s=2.0, c="0.60", alpha=0.6, lw=0, zorder=2)
    ax.set_xlim(*HX); ax.set_ylim(*HY); ax.set_title(title, fontsize=10.5, pad=6)

fig = plt.figure(figsize=(10.5, 7.6))
gs = GridSpec(3, 1, hspace=0.42, left=0.03, right=0.98, top=0.93, bottom=0.02)

# ---- A bounding box ----
axA = fig.add_subplot(gs[0, 0])
axA.add_patch(Rectangle((y0, x0), y1 - y0, x1 - x0, facecolor=IBM["orange"], alpha=0.16, zorder=1))
axA.add_patch(Rectangle((y0, x0), y1 - y0, x1 - x0, fill=False, ec=IBM["orange"], lw=1.8, zorder=3))
base(axA, "Bounding box   %.1f Mpx$^2$  (100%% of box)" % (box_area / mpx)); panel(axA, "A")

# ---- B convex hull ----
axB = fig.add_subplot(gs[1, 0])
axB.add_patch(MplPoly(cvx_d, closed=True, facecolor=IBM["orange"], alpha=0.16, zorder=1))
axB.add_patch(MplPoly(cvx_d, closed=True, fill=False, ec=IBM["orange"], lw=1.8, zorder=3))
base(axB, "Convex hull   %.1f Mpx$^2$  (%.0f%% of box)" % (cvx_area / mpx, 100 * cvx_area / box_area)); panel(axB, "B")

# ---- C concave hull (used) ----
axC = fig.add_subplot(gs[2, 0])
axC.add_patch(MplPoly(hull_d, closed=True, facecolor=IBM["blue"], alpha=0.22, zorder=1))
axC.add_patch(MplPoly(hull_d, closed=True, fill=False, ec=IBM["blue"], lw=1.8, zorder=3))
base(axC, "Concave hull (ratio 0.1, used)   %.1f Mpx$^2$  (%.0f%% of box, %.0f%% of convex)"
     % (cav_area / mpx, 100 * cav_area / box_area, 100 * cav_area / cvx_area)); panel(axC, "C")

fig.savefig(FIG / "pub_window.png", dpi=300, bbox_inches="tight")
print("saved", FIG / "pub_window.png")
print("areas Mpx2  box=%.1f  convex=%.1f  concave=%.1f  | concave/convex=%.0f%%  concave/box=%.0f%%"
      % (box_area / mpx, cvx_area / mpx, cav_area / mpx, 100 * cav_area / cvx_area, 100 * cav_area / box_area))
