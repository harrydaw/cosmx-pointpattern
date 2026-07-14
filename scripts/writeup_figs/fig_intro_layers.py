"""Intro figure — the layered CosMx output for one FOV: morphology, vendor segmentation,
located transcripts, as an exploded diagonal stack with a zoomed transcript inset.
Light / IBM colourblind-safe / Arial to match the write-up figure style (not the dark deck)."""
import sys, os, warnings, contextlib, io
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, ConnectionPatch
import spatialdata as sd

plt.rcParams.update({"font.family": "Arial", "figure.facecolor": "white", "axes.facecolor": "white",
                     "savefig.facecolor": "white"})
IBM = {"blue": "#648FFF", "purple": "#785EF0", "magenta": "#DC267F", "orange": "#FE6100", "yellow": "#FFB000"}
INK, MUTED = "0.15", "0.45"
ROOT = "."
while not os.path.exists(os.path.join(ROOT, "results")):
    ROOT = os.path.join(ROOT, "..")
FIGDIR = os.path.join(ROOT, "results", "figures")

FOV, SUB = "16", 3
with contextlib.redirect_stderr(io.StringIO()):
    s = sd.read_zarr(os.path.join(ROOT, "data", "raw", "Varsha_1234.zarr"))
img = np.asarray(s.images[f"{FOV}_image"].data[:, ::SUB, ::SUB]).astype(float)
lab = np.asarray(s.labels[f"{FOV}_labels"].data[::SUB, ::SUB])
pts = s.points[f"{FOV}_points"]
pts = pts.compute() if hasattr(pts, "compute") else pts
H, W = lab.shape


def norm(c, lo=1, hi=99.5):
    a, b = np.percentile(c, [lo, hi]); return np.clip((c - a) / (b - a + 1e-9), 0, 1)


# morphology: structure (ch0) warm-grey, nuclei (ch4) blue, on a light base
struct, nuc = norm(img[0]), norm(img[4])
morph = np.ones((H, W, 3))
morph[..., 0] = 1 - 0.75 * struct - 0.25 * nuc
morph[..., 1] = 1 - 0.55 * struct - 0.55 * nuc
morph[..., 2] = 1 - 0.20 * struct - 0.05 * nuc
morph = np.clip(morph, 0, 1)

# segmentation: each cell a soft colour, background white
rng = np.random.default_rng(7)
ncell = int(lab.max())
cols = rng.random((ncell + 1, 3)) * 0.5 + 0.3
cols[0] = [1, 1, 1]
segrgb = cols[lab]

# transcripts: top-8 genes coloured (visible on white), the rest light grey
top = pts.target.value_counts().head(8).index.tolist()
palette = [IBM["blue"], IBM["orange"], "#009E73", IBM["purple"], "#56B4E9", IBM["magenta"], "#8C564B", IBM["yellow"]]
gene_c = dict(zip(top, palette))

fig = plt.figure(figsize=(12, 6.2))
h_frac = 0.52
w_frac = h_frac / (12 / 6.2)
x0, y0 = 0.02, 0.30
dx, dy = 0.135, -0.12
panels = [("Morphology", "5-channel stain", morph, None),
          ("Segmentation", f"{ncell} cells", segrgb, None),
          ("Transcripts", f"{len(pts):,} molecules", None, "pts")]
ax_pts = None
for i, (name, sub, rgb, kind) in enumerate(panels):
    ax = fig.add_axes([x0 + i * dx, y0 + i * dy, w_frac, h_frac])
    if kind == "pts":
        ax_pts = ax
        ax.set_xlim(0, W); ax.set_ylim(H, 0)
        other = pts[~pts.target.isin(top)]
        ax.scatter(other.x / SUB, other.y / SUB, s=1.0, c="0.78", alpha=0.7, linewidths=0, rasterized=True)
        for g in top:
            m = pts.target == g
            ax.scatter(pts.x[m] / SUB, pts.y[m] / SUB, s=3.2, c=gene_c[g], alpha=0.9, linewidths=0, rasterized=True)
    else:
        ax.imshow(rgb, origin="upper")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color("0.35"); sp.set_linewidth(1.3)
    ax.text(0.02, -0.02, name, transform=ax.transAxes, color=INK, fontsize=12.5, fontweight="bold", va="top")
    ax.text(0.02, -0.075, sub, transform=ax.transAxes, color=IBM["magenta"], fontsize=10, fontweight="bold", va="top")

# zoomed transcript inset: centred on the densest patch of tissue (2D histogram peak)
zs = int(min(H, W) * 0.20)
hx, hy = pts.x.values / SUB, pts.y.values / SUB
Hh, xe, ye = np.histogram2d(hx, hy, bins=60, range=[[0, W], [0, H]])
bi, bj = np.unravel_index(np.argmax(Hh), Hh.shape)
cxm, cym = (xe[bi] + xe[bi + 1]) / 2, (ye[bj] + ye[bj + 1]) / 2
zx, zy = int(cxm - zs / 2), int(cym - zs / 2)
ax_pts.add_patch(Rectangle((zx, zy), zs, zs, fill=False, ec=INK, lw=1.4))
axz = fig.add_axes([0.70, 0.30, 0.27, 0.52])
axz.set_xlim(zx, zx + zs); axz.set_ylim(zy + zs, zy)
oth = pts[~pts.target.isin(top)]
axz.scatter(oth.x / SUB, oth.y / SUB, s=8, c="0.80", alpha=0.8, linewidths=0)
for g in top:
    m = pts.target == g
    axz.scatter(pts.x[m] / SUB, pts.y[m] / SUB, s=22, c=gene_c[g], alpha=0.95, linewidths=0, label=g)
axz.set_xticks([]); axz.set_yticks([])
for sp in axz.spines.values():
    sp.set_color(INK); sp.set_linewidth(1.4)
axz.text(0.02, 1.008, "Zoomed transcripts (single molecules)", transform=axz.transAxes, color=INK, fontsize=9.5, fontweight="bold", va="bottom")
axz.legend(frameon=False, fontsize=7.2, loc="lower left", ncol=2, handletextpad=0.2, columnspacing=0.7, markerscale=1.1)
# connector lines from the box to the inset
for cxy, ixy in [((zx + zs, zy), (0, 1)), ((zx + zs, zy + zs), (0, 0))]:
    fig.add_artist(ConnectionPatch(xyA=cxy, coordsA=ax_pts.transData, xyB=ixy, coordsB=axz.transAxes,
                                   color="0.6", lw=0.8, ls=(0, (3, 3))))

# small stat call-outs
fig.text(0.71, 0.20, f"{len(pts):,}", color=IBM["orange"], fontsize=26, fontweight="bold", ha="left", va="bottom")
fig.text(0.71, 0.165, "Transcripts", color=MUTED, fontsize=11, ha="left", va="bottom")
fig.text(0.86, 0.20, f"{ncell}", color=IBM["blue"], fontsize=26, fontweight="bold", ha="left", va="bottom")
fig.text(0.86, 0.165, "Cells", color=MUTED, fontsize=11, ha="left", va="bottom")
fig.text(0.71, 0.12, "Segmentation compresses tens of thousands of located\nmolecules into a few dozen cells.",
         color=MUTED, fontsize=9.5, ha="left", va="top", style="italic")

out = os.path.join(FIGDIR, "pub_intro_layers.png")
plt.savefig(out, dpi=300, bbox_inches="tight")
print("saved", out, "| cells", ncell, "| pts", len(pts), "| top genes", top)
