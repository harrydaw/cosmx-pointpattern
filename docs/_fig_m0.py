"""Prototype M0 raw-data 3-panel figure. Iterate here, then port to nb14."""
import warnings, logging
warnings.filterwarnings('ignore'); logging.disable(logging.CRITICAL)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import spatialdata as sd
from pathlib import Path

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 11, 'axes.titlesize': 14, 'axes.titleweight': 'bold',
    'savefig.dpi': 300,
})
C_ASSIGNED = '#648FFF'
C_UNASSIGNED = '#DC267F'
FIGS = Path('results/figures')
FOV = 4

sdata = sd.read_zarr('data/raw/Varsha_1234.zarr')
img = np.asarray(sdata.images[f'{FOV}_image'][0])          # DAPI channel
lab = np.asarray(sdata.labels[f'{FOV}_labels'])            # segmentation mask
pts = sdata.points[f'{FOV}_points'].compute()
H, W = lab.shape
unassigned = pts['cell_ID'] == 0
frac_un = unassigned.mean()
n_cells = int(lab.max())

# DAPI contrast stretch
lo, hi = np.percentile(img, [1, 99.5])
disp = np.clip((img - lo) / (hi - lo + 1e-9), 0, 1)

# random LUT for segmentation labels (0 = background -> dark)
rng = np.random.default_rng(0)
lut = rng.random((n_cells + 1, 3)) * 0.6 + 0.35
lut[0] = (0.07, 0.07, 0.10)
seg_rgb = lut[np.clip(lab, 0, n_cells)]

fig, ax = plt.subplots(1, 3, figsize=(18, 6.6))

ax[0].imshow(disp, cmap='gray', interpolation='nearest')
ax[0].set_title('Raw morphology image\n(DAPI nuclear stain)')

ax[1].imshow(seg_rgb, interpolation='nearest')
ax[1].set_title(f'Vendor cell segmentation\n({n_cells} cells — the layer we discard)')

# transcripts: assigned vs unassigned
ax[2].scatter(pts.loc[~unassigned, 'x'], pts.loc[~unassigned, 'y'],
              s=0.6, c=C_ASSIGNED, marker='.', linewidths=0, rasterized=True,
              label=f'In a cell ({(~unassigned).mean()*100:.1f}%)')
ax[2].scatter(pts.loc[unassigned, 'x'], pts.loc[unassigned, 'y'],
              s=0.6, c=C_UNASSIGNED, marker='.', linewidths=0, rasterized=True,
              label=f'Unassigned, cell_ID=0 ({frac_un*100:.1f}%)')
ax[2].set_xlim(0, W); ax[2].set_ylim(0, H); ax[2].invert_yaxis()
ax[2].set_aspect('equal')
ax[2].set_facecolor('white')
ax[2].set_title('Raw transcripts\n(the layer NoSegs carries through)')
leg = ax[2].legend(loc='lower right', fontsize=10, markerscale=18,
                   framealpha=0.95)

for a in ax:
    a.set_xticks([]); a.set_yticks([])
    for s in a.spines.values():
        s.set_visible(False)

fig.suptitle(
    f'One field of view (FOV {FOV}): segmentation drops {frac_un*100:.0f}% of '
    f'transcripts — NoSegs keeps every point',
    fontsize=16, fontweight='bold', y=1.02)
fig.tight_layout()
out = FIGS / '14_M0_raw_data.png'
fig.savefig(out, dpi=300, bbox_inches='tight')
print('saved', out, '| frac_un', round(frac_un, 3), '| n_cells', n_cells)
