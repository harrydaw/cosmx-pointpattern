"""Intro fig:bottleneck panel B (REAL vendor segmentation). Same real crop of an S1 FOV, two treatments:
LEFT segment+aggregate (real vendor cell masks; off-cell transcripts discarded; -> per-cell matrix);
RIGHT segmentation-free (every transcript kept as a point pattern; radius disc). Light/IBM."""
import warnings, contextlib, io as _io; warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import spatialdata as sd

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "results" / "figures"
plt.rcParams.update({"font.family": "Arial", "figure.facecolor": "white", "axes.facecolor": "white",
                     "savefig.facecolor": "white"})
IBM = {"blue": "#648FFF", "purple": "#785EF0", "magenta": "#DC267F", "orange": "#FE6100"}
INK = "0.15"
FOV = 12; CROP = 320  # full-res px (tight patch, PoC-style)

with contextlib.redirect_stderr(_io.StringIO()):
    sdata = sd.read_zarr(ROOT / "data" / "raw" / "Varsha_1234.zarr")
lab = np.asarray(sdata.labels[f"{FOV}_labels"].data)
pts = sdata.points[f"{FOV}_points"]; pts = pts.compute() if hasattr(pts, "compute") else pts
H, W = lab.shape
# crop centred on the densest transcript region
hb, ye, xe = np.histogram2d(pts.y.values, pts.x.values, bins=30)
iy, ix = np.unravel_index(hb.argmax(), hb.shape)
cy, cx = (ye[iy]+ye[iy+1])/2, (xe[ix]+xe[ix+1])/2
y0, y1 = int(max(0, cy-CROP/2)), int(min(H, cy+CROP/2)); x0, x1 = int(max(0, cx-CROP/2)), int(min(W, cx+CROP/2))
labc = lab[y0:y1, x0:x1]
p = pts[(pts.x >= x0) & (pts.x < x1) & (pts.y >= y0) & (pts.y < y1)].copy()
p["lx"], p["ly"] = p.x - x0, p.y - y0
p["assigned"] = labc[np.clip(p.ly.astype(int), 0, labc.shape[0]-1), np.clip(p.lx.astype(int), 0, labc.shape[1]-1)] > 0
ncell = len(np.unique(labc[labc > 0]))
print(f"FOV {FOV} crop {labc.shape} | transcripts {len(p)} | off-cell {(~p.assigned).mean()*100:.0f}% | cells {ncell}")

# pastel colour per cell (curated, non-orange, to avoid clashing with the off-cell highlight)
pastels = np.array([[0.72,0.80,0.96],[0.76,0.90,0.79],[0.85,0.79,0.94],[0.73,0.89,0.89],
                    [0.90,0.80,0.90],[0.80,0.86,0.73],[0.83,0.84,0.95],[0.90,0.86,0.75]])
nlab = int(labc.max())
cols = np.vstack([[1,1,1], pastels[np.arange(nlab) % len(pastels)]])
seg = cols[labc]

fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 5.6))

# ---------- LEFT: real segmentation + aggregate ----------
axL.imshow(seg, origin="upper", interpolation="nearest")
asg, una = p[p.assigned], p[~p.assigned]
axL.scatter(asg.lx, asg.ly, s=4, c="0.30", alpha=0.7, lw=0, zorder=2)
axL.scatter(una.lx, una.ly, s=7, c=IBM["orange"], alpha=0.75, lw=0, zorder=3)
axL.set_xlim(0, labc.shape[1]); axL.set_ylim(labc.shape[0], 0); axL.axis("off")
axL.set_title("Standard: segment, then aggregate", fontsize=11.5, fontweight="bold", pad=10)
axL.text(0.5, -0.03, f"Grey = in a vendor cell; orange = off-cell, discarded ({(~p.assigned).mean()*100:.0f}%)",
         transform=axL.transAxes, ha="center", va="top", fontsize=9, color=INK)

# ---------- RIGHT: segmentation-free ----------
axR.set_xlim(0, labc.shape[1]); axR.set_ylim(labc.shape[0], 0); axR.set_aspect("equal"); axR.axis("off")
axR.set_title("Segmentation-free: keep the point pattern", fontsize=11.5, fontweight="bold", pad=10)
top2 = p.target.value_counts().head(2).index.tolist()
oth = p[~p.target.isin(top2)]
axR.scatter(oth.lx, oth.ly, s=6, c="0.72", alpha=0.6, lw=0, zorder=1)
for g, col in zip(top2, [IBM["blue"], IBM["magenta"]]):
    m = p.target == g; axR.scatter(p.lx[m], p.ly[m], s=13, c=col, alpha=0.95, lw=0, zorder=3, label=g)
fpt = p[p.target == top2[0]][["lx", "ly"]].values
fp = fpt[len(fpt)//2] if len(fpt) else np.array([labc.shape[1]/2, labc.shape[0]/2])
axR.add_patch(Circle(fp, 70, fill=False, ec=INK, lw=1.5, ls="--", zorder=4))
axR.annotate("r", fp, fp + np.array([42, -42]), color=INK, fontsize=12,
             arrowprops=dict(arrowstyle="<-", color=INK, lw=1.2))
axR.text(0.5, -0.03, "Every transcript kept; spatial relations measured directly",
         transform=axR.transAxes, ha="center", va="top", fontsize=9, color=INK)
axR.legend(frameon=False, fontsize=8.5, loc="upper right", title="Most abundant genes", title_fontsize=8)

fig.subplots_adjust(wspace=0.12, left=0.02, right=0.98, top=0.9, bottom=0.10)
fig.savefig(FIG / "pub_concept.png", dpi=300, bbox_inches="tight"); print("saved", FIG / "pub_concept.png")
