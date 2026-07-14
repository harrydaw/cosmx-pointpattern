"""Methods fig:dataprep — data preparation on ONE field of view (nb10 style):
A initial data | B DBSCAN noise + manual QC highlighted | C GMM strip assignment
histogram | D final strip-assigned data. Vertical strips, light/IBM/Arial, real S1 data."""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "results" / "figures"
plt.rcParams.update({"font.family": "Arial", "figure.facecolor": "white", "axes.facecolor": "white",
                     "savefig.facecolor": "white"})
IBM = {"blue": "#648FFF", "purple": "#785EF0", "magenta": "#DC267F", "orange": "#FE6100", "yellow": "#FFB000"}
STRIP = {"strip_1": IBM["blue"], "strip_2": IBM["magenta"], "strip_3": IBM["purple"]}
SLAB = {"strip_1": "Strip 1 (control)", "strip_2": "Strip 2 (infected)", "strip_3": "Strip 3 (control)"}
MANRED = "#A11526"
def title(ax, letter, txt): ax.set_title(r"$\bf{%s}$   %s" % (letter, txt), fontsize=11, loc="left", pad=6)

XR, YR = "x_rot_px", "y_rot_px"           # x_rot = across-strip (left-right); y_rot = length (up)
k = pd.read_parquet(ROOT / "data" / "processed" / "s1_all_strips_cleaned.parquet")
k["removed"] = k.is_noise | k.is_small_cluster.fillna(False) | k.manually_excluded.fillna(False)
def PX(df): return df[XR].values
def PY(df): return df[YR].values

FOV = 8
f = k[k.fov == FOV]
fk = f[~f.removed]
x0, x1 = f[XR].min(), f[XR].max(); y0, y1 = f[YR].min(), f[YR].max()
mgx = 0.03 * (x1 - x0); mgy = 0.03 * (y1 - y0)
XL, YL = (x0 - mgx, x1 + mgx), (y0 - mgy, y1 + mgy)
print(f"FOV {FOV}: n={len(f)} kept={len(fk)} ({len(fk)/len(f)*100:.0f}%) strips={f.strip.value_counts().to_dict()}")

# GMM on the clean across-strip coordinate for this FOV
xr = fk[XR].values
gm = GaussianMixture(3, random_state=42).fit(xr.reshape(-1, 1))
order = np.argsort(gm.means_.ravel())
mu = gm.means_.ravel()[order]; sd = np.sqrt(gm.covariances_.ravel())[order]; w = gm.weights_[order]
bnd = [(mu[i] + mu[i + 1]) / 2 for i in range(2)]

fig = plt.figure(figsize=(10.6, 10.0))
gs = GridSpec(2, 2, height_ratios=[1.0, 1.0], hspace=0.16, wspace=0.14,
              left=0.03, right=0.98, top=0.95, bottom=0.05)

def scat(ax):
    ax.set_aspect("equal"); ax.axis("off"); ax.set_xlim(*XL); ax.set_ylim(*YL)

# ===== A: initial data =====
axA = fig.add_subplot(gs[0, 0]); scat(axA)
axA.scatter(PX(f), PY(f), s=3, c="0.6", alpha=0.5, lw=0, rasterized=True)
title(axA, "A", "Initial data")

# ===== B: DBSCAN noise + manual QC highlighted =====
axB = fig.add_subplot(gs[0, 1]); scat(axB)
nz = f.is_noise.fillna(False)
sc = f.is_small_cluster.fillna(False) & ~nz
mn = f.manually_excluded.fillna(False) & ~nz & ~sc
axB.scatter(PX(fk), PY(fk), s=3, c="0.72", alpha=0.6, lw=0, rasterized=True)
axB.scatter(PX(f[nz]), PY(f[nz]), s=3, c=IBM["orange"], alpha=0.5, lw=0, rasterized=True)
axB.scatter(PX(f[sc]), PY(f[sc]), s=7, c=IBM["yellow"], alpha=0.9, lw=0, rasterized=True)
axB.scatter(PX(f[mn]), PY(f[mn]), s=7, c=MANRED, alpha=0.85, lw=0, rasterized=True)
title(axB, "B", "DBSCAN + manual QC")
axB.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor="0.72", markersize=6, label="Retained"),
                    Line2D([0], [0], marker="o", color="w", markerfacecolor=IBM["orange"], markersize=6, label="DBSCAN noise"),
                    Line2D([0], [0], marker="o", color="w", markerfacecolor=IBM["yellow"], markersize=6, label="Small cluster"),
                    Line2D([0], [0], marker="o", color="w", markerfacecolor=MANRED, markersize=6, label="Manual removal")],
           loc="lower center", frameon=False, fontsize=8, ncol=2, bbox_to_anchor=(0.5, -0.02))

# ===== C: GMM strip assignment histogram =====
axC = fig.add_subplot(gs[1, 0])
bins = np.linspace(x0, x1, 70)
axC.hist([fk[fk.strip == s][XR].values for s in STRIP], bins=bins, stacked=True,
         color=[STRIP[s] for s in STRIP], alpha=0.6, edgecolor="none")
xg = np.linspace(x0, x1, 400); bw = bins[1] - bins[0]; N = len(fk)
for m_, s_, w_, st in zip(mu, sd, w, STRIP):
    axC.plot(xg, N * bw * w_ * norm.pdf(xg, m_, s_), color=STRIP[st], lw=2.0)
for b in bnd:
    axC.axvline(b, color="0.25", lw=1.2, ls=(0, (3, 2)))
axC.set_xlabel("Across-strip coordinate (px, PCA-rotated)", fontsize=9.5)
axC.set_ylabel("Transcripts", fontsize=9.5); axC.set_yticks([]); axC.tick_params(labelsize=8)
for sp in ("top", "right"): axC.spines[sp].set_visible(False)
title(axC, "C", "GMM strip assignment")

# ===== D: final strip-assigned data =====
axD = fig.add_subplot(gs[1, 1]); scat(axD)
for st, c in STRIP.items():
    m = fk.strip.values == st
    axD.scatter(PX(fk)[m], PY(fk)[m], s=3, c=c, alpha=0.6, lw=0, rasterized=True)
title(axD, "D", "Strip-assigned data")
axD.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=STRIP[s], markersize=6, label=SLAB[s])
                    for s in STRIP], loc="lower center", frameon=False, fontsize=8, ncol=1, bbox_to_anchor=(0.5, -0.02))

fig.savefig(FIG / "pub_dataprep.png", dpi=300, bbox_inches="tight")
print("saved", FIG / "pub_dataprep.png", "| GMM means", np.round(mu).tolist(), "| boundaries", np.round(bnd).tolist())
