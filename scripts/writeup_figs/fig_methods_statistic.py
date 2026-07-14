"""Methods fig:statistic — the bivariate statistic illustrated on SPP1 x CD44, strip 2 (real data).
A: dense tissue patch (all transcripts grey = structure) with SPP1/CD44 highlighted | B: pair-counting
within r | C: edge-correction disc fraction | D: observed L(r) vs envelope | E: SES(r). No step-legend panel."""
import warnings; warnings.filterwarnings("ignore")
import json, io
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Circle, Polygon as MplPoly, Rectangle
from shapely.geometry import Point

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "results" / "figures"
plt.rcParams.update({"font.family": "Arial", "figure.facecolor": "white", "axes.facecolor": "white",
                     "savefig.facecolor": "white", "axes.linewidth": 0.8, "figure.dpi": 110,
                     "axes.spines.top": False, "axes.spines.right": False})
IBM = {"blue": "#648FFF", "purple": "#785EF0", "magenta": "#DC267F", "orange": "#FE6100", "yellow": "#FFB000"}
BAND, INK, MUTED = "0.86", "0.15", "0.45"
def panel(ax, L): ax.text(-0.02, 1.04, L, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom", ha="right")

ns = globals()
nb = json.load(io.open(ROOT / "notebooks" / "00_functions.ipynb", encoding="utf-8"))
for c in nb["cells"]:
    if c["cell_type"] == "code":
        try: exec("".join(c["source"]), ns)
        except Exception: pass
assert "get_concave_hull" in ns and "fraction_inside_hull" in ns

R = np.load(ROOT / "results" / "lr_panel_results.parquet.r_vals.npy")
GA, GB, ST = "SPP1", "CD44", "strip_2"
XC, YC = "x_global_px_transformed", "y_global_px_transformed"
k = pd.read_parquet(ROOT / "data" / "processed" / "s1_all_strips_cleaned.parquet")
keep = (~k.is_noise) & (~k.is_small_cluster.fillna(False)) & (~k.get("manually_excluded", pd.Series(False, index=k.index)).fillna(False))
sdf = k[keep & (k.strip == ST)]
allxy = sdf[[XC, YC]].values.astype(float)
A = sdf[sdf.target == GA][[XC, YC]].values.astype(float)
B = sdf[sdf.target == GB][[XC, YC]].values.astype(float)
P = np.vstack([A, B])
# patch (S x S) centred to contain the most SPP1/CD44 while showing dense tissue structure
S = 1000.0; hlf = S / 2
cnt_box = [((np.abs(P[:,0]-p[0]) < hlf) & (np.abs(P[:,1]-p[1]) < hlf)).sum() for p in P]
cen = P[int(np.argmax(cnt_box))]
X0, X1, Y0, Y1 = cen[0]-hlf, cen[0]+hlf, cen[1]-hlf, cen[1]+hlf
inpatch = lambda XY: (XY[:,0] > X0) & (XY[:,0] < X1) & (XY[:,1] > Y0) & (XY[:,1] < Y1)
pAll, pA, pB = allxy[inpatch(allxy)], A[inpatch(A)], B[inpatch(B)]
print(f"patch: {len(pAll)} transcripts | {GA} {len(pA)} | {GB} {len(pB)}")

# focal SPP1 (in patch) with most CD44 within rB -> panel B zoom
rB, halfB = 70.0, 175.0
cnt = np.array([(np.hypot(pB[:,0]-p[0], pB[:,1]-p[1]) <= rB).sum() for p in pA]) if len(pA) else np.array([0])
focus = pA[int(np.argmax(cnt))] if len(pA) else cen

fig = plt.figure(figsize=(9.5, 8.6))
gs = GridSpec(2, 2, height_ratios=[1.08, 0.92], hspace=0.32, wspace=0.24,
              left=0.07, right=0.98, top=0.93, bottom=0.10)

# ===== A: dense tissue patch, SPP1/CD44 highlighted =====
axA = fig.add_subplot(gs[0, 0])
axA.scatter(*pAll.T, s=2.2, c="0.68", alpha=0.5, lw=0, zorder=1)          # all transcripts = structure
axA.scatter(*pA.T, s=20, c=IBM["blue"], alpha=0.95, lw=0, zorder=3, label=f"{GA} (n={len(pA)})")
axA.scatter(*pB.T, s=20, c=IBM["magenta"], alpha=0.95, lw=0, zorder=3, label=f"{GB} (n={len(pB)})")
axA.add_patch(Rectangle((focus[0]-halfB, focus[1]-halfB), 2*halfB, 2*halfB, fill=False, ec=INK, lw=1.2, ls="--", zorder=4))
axA.set_xlim(X0, X1); axA.set_ylim(Y0, Y1); axA.set_aspect("equal"); axA.axis("off"); panel(axA, "A")
axA.legend(frameon=False, fontsize=8.5, loc="upper center", ncol=2, markerscale=1.6, bbox_to_anchor=(0.5, -0.02))
axA.set_title("Tissue patch (infected strip)\nGrey = all transcripts; box = panel B", fontsize=9.5)

# ===== B: pair-counting within r =====
axB = fig.add_subplot(gs[0, 1])
selA = pA[(np.abs(pA[:,0]-focus[0]) < halfB) & (np.abs(pA[:,1]-focus[1]) < halfB)]
selB = pB[(np.abs(pB[:,0]-focus[0]) < halfB) & (np.abs(pB[:,1]-focus[1]) < halfB)]
within = selB[np.hypot(selB[:,0]-focus[0], selB[:,1]-focus[1]) <= rB]
pAllB = pAll[(np.abs(pAll[:,0]-focus[0]) < halfB) & (np.abs(pAll[:,1]-focus[1]) < halfB)]
axB.scatter(*pAllB.T, s=6, c="0.82", alpha=0.6, lw=0, zorder=0)   # tissue context
axB.scatter(*selA.T, s=16, c=IBM["blue"], alpha=0.45, lw=0)
axB.scatter(*selB.T, s=16, c=IBM["magenta"], alpha=0.4, lw=0)
axB.add_patch(Circle(focus, rB, fill=False, ec=INK, lw=1.6, ls="--"))
axB.scatter(*within.T, s=30, c=IBM["magenta"], edgecolors="white", lw=0.5, zorder=4)
axB.scatter(*focus, s=80, c=IBM["blue"], edgecolors="white", lw=0.8, zorder=5)
axB.annotate("r", focus, focus + np.array([rB*0.5, rB*0.5]), color=INK, fontsize=12,
             arrowprops=dict(arrowstyle="<-", color=INK, lw=1.2))
axB.set_xlim(focus[0]-halfB, focus[0]+halfB); axB.set_ylim(focus[1]-halfB, focus[1]+halfB)
axB.set_aspect("equal"); axB.axis("off"); panel(axB, "B")
axB.set_title(f"Count {GB} within r of each {GA}\n({len(within)} within r = {rB:.0f} px here)", fontsize=9.5)

# ===== D, E from the screen output =====
d = pd.concat([pd.read_parquet(ROOT/"results"/"lr_panel_results.parquet"),
               pd.read_parquet(ROOT/"results"/"lr_panel_results_v2.parquet")], ignore_index=True)
row = d[(d.ligand==GA)&(d.receptor==GB)&(d.strip==ST)].iloc[0]
lo, hi, ob = np.array(row.l_lo,float), np.array(row.l_hi,float), np.array(row.l_obs,float)
mid = (hi+lo)/2; Hh = np.where((hi-lo)/2==0, np.nan, (hi-lo)/2); ses = (ob-mid)/Hh
axD = fig.add_subplot(gs[1, 0])
axD.fill_between(R, lo, hi, color=BAND, label="Permutation envelope")
axD.plot(R, ob, color=IBM["blue"], lw=2.0, label="Observed L(r)")
exc = ob > hi
if exc.any(): axD.scatter(R[exc], ob[exc], s=14, c=INK, zorder=5, label="Above envelope")
axD.set_xlabel("r (px)", fontsize=10); axD.set_ylabel("L(r)", fontsize=10); axD.tick_params(labelsize=8)
axD.legend(frameon=False, fontsize=8, loc="upper left"); panel(axD, "C")
axD.set_title(f"{GA} × {GB}, Strip 2 (infected)", fontsize=9.5)

axE = fig.add_subplot(gs[1, 1])
axE.axhspan(-1, 1, color="0.92", zorder=0)
axE.axhline(0, color="0.55", ls="--", lw=0.9); axE.axhline(1, color="0.7", ls=":", lw=0.8)
axE.plot(R, ses, color=IBM["magenta"], lw=2.0)
axE.scatter(R[np.nanargmax(ses)], np.nanmax(ses), s=40, c=INK, zorder=5)
axE.annotate(f"Peak SES {np.nanmax(ses):.2f}", (R[np.nanargmax(ses)], np.nanmax(ses)),
             textcoords="offset points", xytext=(-8, 8), fontsize=8.5, ha="right", color=INK)
axE.set_xlabel("r (px)", fontsize=10); axE.set_ylabel("SES(r)", fontsize=10); axE.tick_params(labelsize=8)
panel(axE, "D"); axE.set_title("Standardised effect size", fontsize=9.5)

fig.savefig(FIG / "pub_statistic.png", dpi=300, bbox_inches="tight")
print("saved", FIG / "pub_statistic.png")
