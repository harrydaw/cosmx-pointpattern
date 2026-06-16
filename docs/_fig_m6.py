"""Prototype M6 workflow flowchart + quick-start. Iterate, then port to nb14."""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'savefig.dpi': 300,
})
FIGS = Path('results/figures')

PURPLE = '#500778'
PURPLE_L = '#EDE3F3'
STAGE = '#6A1B9A'
OPT_BG = '#FFF4D6'
OPT_ED = '#E0A800'
IN_BG = '#1A9850'
OUT_BG = '#2166AC'
TXT = '#2C2C2C'

fig, ax = plt.subplots(figsize=(11, 15))
ax.set_xlim(0, 100); ax.set_ylim(0, 140)
ax.axis('off')

def box(x, y, w, h, text, fc, ec, tc='white', fs=14, bold=True, rad=0.08):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0.2,rounding_size={rad*40}',
                       linewidth=1.6, edgecolor=ec, facecolor=fc, mutation_aspect=1)
    ax.add_patch(p)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', color=tc,
            fontsize=fs, fontweight='bold' if bold else 'normal', wrap=True)

def opt(x, y, w, h, text):
    p = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.15,rounding_size=2',
                       linewidth=1.2, edgecolor=OPT_ED, facecolor=OPT_BG)
    ax.add_patch(p)
    ax.text(x + w/2, y + h/2, text, ha='center', va='center', color='#6B5200',
            fontsize=10.5, fontstyle='italic')

def arrow(y0, y1, x=33):
    ax.add_patch(FancyArrowPatch((x, y0), (x, y1), arrowstyle='-|>',
                 mutation_scale=22, linewidth=2.2, color=PURPLE))

ax.text(50, 137, 'The NoSegs workflow', ha='center', fontsize=20,
        fontweight='bold', color=PURPLE)
ax.text(50, 133, 'from raw transcripts to signalling networks — no segmentation',
        ha='center', fontsize=12, fontstyle='italic', color='#6A1B9A')
ax.text(78.5, 129.5, 'amber = options / choices you can tune',
        ha='center', fontsize=10, fontstyle='italic', color='#6B5200')

# stage column geometry
bx, bw, bh = 8, 50, 8.5
ox, ow, oh = 62, 33, 7.0
stages = [
    ('INPUT:  raw transcripts  (x, y, gene, FOV)', IN_BG, IN_BG, None),
    ('1.  PCA align + per-FOV GMM\nsplit into tissue strips', STAGE, PURPLE, '3-component\nGaussian mixture'),
    ('2.  DBSCAN noise removal', STAGE, PURPLE, 'adaptive ε =\np97 1-NN, [20,30] px'),
    ('3.  Observation window', STAGE, PURPLE, 'rectangle / convex /\nconcave hull'),
    ('4.  Bivariate Ripley K → L(r)', STAGE, PURPLE, 'multi-scale r-grid\n(cell + tissue)'),
    ('5.  Label-swap permutation\nenvelope', STAGE, PURPLE, 'n_sim = 99 / 199\n+ edge correction'),
    ('6.  Effect-size ranking (SES)', STAGE, PURPLE, 'binary test  or\ncontinuous SES'),
    ('OUTPUTS:  networks · UMAP · heatmaps', OUT_BG, OUT_BG, None),
]
y = 122
gap = 15.2
for i, (label, fc, ec, option) in enumerate(stages):
    box(bx, y, bw, bh, label, fc, ec, fs=13.5)
    if option:
        opt(ox, y + (bh-oh)/2, ow, oh, option)
        ax.add_patch(FancyArrowPatch((bx+bw, y+bh/2), (ox, y+oh/2+ (bh-oh)/2),
                     arrowstyle='-', linewidth=1.0, color=OPT_ED, linestyle=(0,(3,2))))
    if i < len(stages) - 1:
        arrow(y - 0.5, y - gap + bh + 0.5, x=bx + bw/2)
    y -= gap

# Quick-start panel (bottom)
qx, qy, qw, qh = 4, 0.5, 92, 13.0
p = FancyBboxPatch((qx, qy), qw, qh, boxstyle='round,pad=0.3,rounding_size=3',
                   linewidth=1.8, edgecolor=PURPLE, facecolor=PURPLE_L)
ax.add_patch(p)
ax.text(qx + 3, qy + qh - 1.9, 'Quick start — use NoSegs when:', fontsize=14,
        fontweight='bold', color=PURPLE, va='center')
uses = [
    'segmentation fails or is unavailable on your sample',
    'you want ligand–receptor co-localisation straight from transcripts',
    'imaging-based platforms (CosMx, MERFISH, Xenium)',
    '≥ ~50 transcripts per gene per region of interest',
]
for j, u in enumerate(uses):
    cy = qy + qh - 4.4 - 2.2 * j
    ax.text(qx + 3, cy, '✓', fontsize=12, color=IN_BG, va='center',
            fontweight='bold', fontfamily='DejaVu Sans')
    ax.text(qx + 6, cy, u, fontsize=11.5, color=TXT, va='center')

fig.savefig(FIGS / '14_M6_workflow.png', dpi=300, bbox_inches='tight')
print('saved 14_M6_workflow.png')
