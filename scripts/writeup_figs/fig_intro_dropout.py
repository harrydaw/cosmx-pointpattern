"""Intro fig:bottleneck panel A: vendor segmentation dropout, light/IBM style, no baked title.
Overall assigned/unassigned + per-FOV consistency. Numbers from the production parquet."""
import warnings; warnings.filterwarnings("ignore")
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "results" / "figures"
plt.rcParams.update({"font.family": "Arial", "figure.facecolor": "white", "axes.facecolor": "white",
                     "savefig.facecolor": "white", "axes.linewidth": 0.8, "figure.dpi": 110,
                     "axes.spines.top": False, "axes.spines.right": False})
IBM = {"blue": "#648FFF", "purple": "#785EF0", "magenta": "#DC267F", "orange": "#FE6100", "yellow": "#FFB000"}
def panel(ax, L): ax.text(-0.02, 1.05, L, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom", ha="right")

df = pd.read_parquet(ROOT / "data" / "processed" / "s1_all_strips_cleaned.parquet", columns=["cell_ID", "fov"])
N = len(df); n_un = int((df.cell_ID == 0).sum()); n_as = N - n_un
per = df.groupby("fov").cell_ID.apply(lambda s: (s == 0).mean() * 100)
mean_un = n_un / N * 100
print(f"N={N} assigned={n_as} ({n_as/N*100:.1f}%) unassigned={n_un} ({mean_un:.1f}%) | FOVs={df.fov.nunique()}")

fig, (a0, a1) = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw={"width_ratios": [1, 1.5], "wspace": 0.3})
# overall
bars = a0.bar(["assigned\n(cell ID > 0)", "unassigned\n(cell ID = 0)"], [n_as, n_un],
              color=[IBM["blue"], IBM["orange"]], width=0.62)
for r, n in zip(bars, [n_as, n_un]):
    a0.text(r.get_x() + r.get_width()/2, n + N*0.012, f"{n:,}\n({n/N*100:.1f}%)", ha="center", va="bottom", fontsize=9)
a0.set_ylabel("Transcripts", fontsize=10); a0.set_ylim(0, N*0.62); a0.tick_params(labelsize=9)
a0.set_yticks([0, 1e5, 2e5, 3e5, 4e5]); a0.set_yticklabels(["0", "100k", "200k", "300k", "400k"])
# per-FOV
a1.bar(per.index.astype(str), per.values, color=IBM["orange"], width=0.7)
a1.axhline(mean_un, color="0.25", ls="--", lw=1.1, label=f"Mean {mean_un:.1f}%")
a1.set_xlabel("Field of view", fontsize=10); a1.set_ylabel("% unassigned", fontsize=10)
a1.set_ylim(0, 60); a1.tick_params(labelsize=9); a1.legend(frameon=False, fontsize=9, loc="upper right")
plt.tight_layout(); fig.savefig(FIG / "pub_dropout.png", dpi=300, bbox_inches="tight")
print("saved", FIG / "pub_dropout.png")
