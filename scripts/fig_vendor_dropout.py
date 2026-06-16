"""Figure: vendor segmentation dropout (cell_ID == 0) on the full dataset.

Produces results/figures/01_vendor_segmentation_dropout.png — a two-panel
defence of the segmentation-free choice:
  (a) overall bar chart: assigned vs unassigned across the whole dataset
  (b) per-FOV bar chart: unassigned fraction per FOV, to show 42% is not
      driven by a handful of bad FOVs

Run from repo root:
    .venv/Scripts/python.exe scripts/fig_vendor_dropout.py
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO   = Path(__file__).resolve().parents[1]
DATA   = REPO / 'data' / 'processed' / 's1_all_strips_cleaned.parquet'
OUT    = REPO / 'results' / 'figures' / '01_vendor_segmentation_dropout.png'
OUT.parent.mkdir(parents=True, exist_ok=True)

# cell_ID is the vendor's segmentation assignment (0 = unassigned).
# Our QC flags (is_noise / is_small_cluster / manually_excluded) are
# downstream of and orthogonal to this column, so it's valid to read
# vendor-segmentation quality off the cleaned parquet.
df = pd.read_parquet(DATA, columns=['cell_ID', 'fov'])

n_total      = len(df)
n_unassigned = int((df['cell_ID'] == 0).sum())
n_assigned   = n_total - n_unassigned
pct_unassigned = 100 * n_unassigned / n_total

print(f'Total transcripts:        {n_total:,}')
print(f'Unassigned (cell_ID=0):   {n_unassigned:,}  ({pct_unassigned:.1f}%)')
print(f'Assigned (cell_ID > 0):   {n_assigned:,}    ({100 - pct_unassigned:.1f}%)')

# Per-FOV breakdown
per_fov = (
    df.assign(unassigned=df['cell_ID'].eq(0))
      .groupby('fov')
      .agg(n_total=('unassigned', 'size'),
           n_unassigned=('unassigned', 'sum'))
)
per_fov['pct_unassigned'] = 100 * per_fov['n_unassigned'] / per_fov['n_total']
per_fov = per_fov.sort_index()
print('\nPer-FOV unassigned %:')
print(per_fov[['n_total', 'n_unassigned', 'pct_unassigned']].to_string())

# ── Figure ────────────────────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5),
                               gridspec_kw={'width_ratios': [1, 2]})

# (a) Overall stacked bar
counts  = [n_assigned, n_unassigned]
labels  = ['Assigned\n(cell_ID > 0)', 'Unassigned\n(cell_ID = 0)']
colors  = ['#4C72B0', '#C44E52']
bars    = ax1.bar(labels, counts, color=colors, width=0.6, edgecolor='black', linewidth=0.5)

for bar, n in zip(bars, counts):
    pct = 100 * n / n_total
    ax1.text(bar.get_x() + bar.get_width() / 2,
             bar.get_height() + n_total * 0.01,
             f'{n:,}\n({pct:.1f}%)',
             ha='center', va='bottom', fontsize=10)

ax1.set_ylabel('Transcript count')
ax1.set_title(f'Vendor segmentation drops {pct_unassigned:.0f}% of transcripts\n'
              f'(total N = {n_total:,})', fontsize=11)
ax1.set_ylim(0, max(counts) * 1.18)
ax1.grid(axis='y', alpha=0.3)
ax1.set_axisbelow(True)

# (b) Per-FOV unassigned %
fov_labels = [str(f) for f in per_fov.index]
ax2.bar(fov_labels, per_fov['pct_unassigned'],
        color='#C44E52', edgecolor='black', linewidth=0.5, width=0.7)
ax2.axhline(pct_unassigned, color='black', lw=1, ls='--',
            label=f'Dataset mean ({pct_unassigned:.1f}%)')
ax2.set_xlabel('FOV')
ax2.set_ylabel('% unassigned (cell_ID = 0)')
ax2.set_title('Per-FOV breakdown — dropout is consistent across FOVs',
              fontsize=11)
ax2.set_ylim(0, max(per_fov['pct_unassigned'].max() * 1.15, 60))
ax2.legend(loc='upper right', fontsize=9)
ax2.grid(axis='y', alpha=0.3)
ax2.set_axisbelow(True)

plt.tight_layout()
plt.savefig(OUT, dpi=160, bbox_inches='tight')
print(f'\nSaved -> {OUT.relative_to(REPO)}')
