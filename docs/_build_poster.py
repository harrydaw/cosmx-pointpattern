"""Build the A0 NeuroMonster poster from the methodology-forward layout.

Produces docs/NoSeggs_Poster_A0.pptx — fully editable native PowerPoint
(text boxes, table, picture objects). Figures inserted as pictures; all
text/guides/legends are editable. Run with the repo venv python.
"""
from pathlib import Path
from PIL import Image
from pptx import Presentation
from pptx.util import Cm, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from pptx.oxml.ns import qn

REPO = Path(__file__).resolve().parent.parent
FIGS = REPO / "results" / "figures"
ASSETS = REPO / "docs" / "_poster_assets"
OUT = REPO / "docs" / "NoSeggs_Poster_A0.pptx"

# ---- palette ----
PURPLE   = RGBColor(0x50, 0x07, 0x78)   # KCL purple, header/footer/bottom band
PANEL    = RGBColor(0xF3, 0xEE, 0xF7)   # light column panel
HEADING  = RGBColor(0x50, 0x07, 0x78)
SUBHEAD  = RGBColor(0x6A, 0x1B, 0x9A)
BODY     = RGBColor(0x2C, 0x2C, 0x2C)
META     = RGBColor(0xC9, 0xA8, 0xD4)
WHITE    = RGBColor(0xFF, 0xFF, 0xFF)
GUIDE    = RGBColor(0x9A, 0x7A, 0xAE)   # grey-purple for fill-in guides
LEGEND   = RGBColor(0x3A, 0x3A, 0x3A)
RULE     = RGBColor(0xC9, 0xA8, 0xD4)
F = "Arial"

prs = Presentation()
prs.slide_width  = Cm(84.1)
prs.slide_height = Cm(118.9)
slide = prs.slides.add_slide(prs.slide_layouts[6])  # blank
shapes = slide.shapes


def _no_line(sp):
    sp.line.fill.background()


def rect(x, y, w, h, fill, rounded=False, line=None):
    shp = MSO_SHAPE.ROUNDED_RECTANGLE if rounded else MSO_SHAPE.RECTANGLE
    sp = shapes.add_shape(shp, Cm(x), Cm(y), Cm(w), Cm(h))
    sp.fill.solid()
    sp.fill.fore_color.rgb = fill
    if line is None:
        _no_line(sp)
    else:
        sp.line.color.rgb = line
        sp.line.width = Pt(1)
    sp.shadow.inherit = False
    if rounded:
        try:
            sp.adjustments[0] = 0.03
        except Exception:
            pass
    return sp


def textbox(x, y, w, h, paras, anchor=MSO_ANCHOR.TOP, align=PP_ALIGN.LEFT,
            wrap=True):
    """paras: list of dicts with keys:
       text, size, bold, italic, color, bullet(bool), align, space_after,
       runs (optional list of (text, dict-overrides)), line_spacing
    """
    tb = shapes.add_textbox(Cm(x), Cm(y), Cm(w), Cm(h))
    tf = tb.text_frame
    tf.word_wrap = wrap
    tf.vertical_anchor = anchor
    tf.margin_left = Cm(0.1); tf.margin_right = Cm(0.1)
    tf.margin_top = Cm(0.05); tf.margin_bottom = Cm(0.05)
    for i, pdef in enumerate(paras):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = pdef.get("align", align)
        if pdef.get("space_after") is not None:
            p.space_after = Pt(pdef["space_after"])
        if pdef.get("space_before") is not None:
            p.space_before = Pt(pdef["space_before"])
        if pdef.get("line_spacing"):
            p.line_spacing = pdef["line_spacing"]
        runs = pdef.get("runs")
        if runs is None:
            runs = [(pdef.get("text", ""), {})]
        bullet = pdef.get("bullet", False)
        for txt, ov in runs:
            r = p.add_run()
            r.text = txt
            r.font.name = ov.get("font", F)
            r.font.size = Pt(ov.get("size", pdef.get("size", 16)))
            r.font.bold = ov.get("bold", pdef.get("bold", False))
            r.font.italic = ov.get("italic", pdef.get("italic", False))
            r.font.color.rgb = ov.get("color", pdef.get("color", BODY))
        if bullet:
            _set_bullet(p)
    return tb


def _set_bullet(p):
    pPr = p._pPr
    if pPr is None:
        pPr = p._p.get_or_add_pPr()
    pPr.set("indent", "-228600")
    pPr.set("marL", "228600")
    buFont = pPr.makeelement(qn("a:buFont"), {"typeface": "Arial"})
    buChar = pPr.makeelement(qn("a:buChar"), {"char": "•"})
    pPr.append(buFont)
    pPr.append(buChar)


def figure(path, x_col, colw, y, target_w=None, center=True):
    """Place picture; width target_w (cm) or colw. Returns bottom y (cm)."""
    im = Image.open(path)
    ar = im.size[0] / im.size[1]
    w = target_w if target_w else colw
    h = w / ar
    x = x_col + (colw - w) / 2 if center else x_col
    shapes.add_picture(str(path), Cm(x), Cm(y), width=Cm(w))
    return y + h


def heading(text, x, w, y, size=30):
    textbox(x, y, w, 1.5, [{"text": text, "size": size, "bold": True,
                             "color": HEADING}])
    # underline rule
    rl = rect(x, y + 1.45, w, 0.07, PURPLE)
    return y + 1.7


def subhead(text, x, w, y, size=21):
    textbox(x, y, w, 1.1, [{"text": text, "size": size, "bold": True,
                             "color": SUBHEAD}])
    return y + 1.15


def legend(figtag, captext, x, w, y, h):
    textbox(x, y, w, h, [{
        "runs": [(figtag + " ", {"bold": True, "size": 14, "color": HEADING}),
                 (captext, {"size": 14, "color": LEGEND})],
        "line_spacing": 1.02,
    }])
    return y + h


# ================= HEADER =================
rect(0, 0, 84.1, 12.4, PURPLE)
textbox(2.0, 1.0, 70, 5.2, [{
    "text": "NoSeggs: Segmentation-Free Ligand–Receptor Co-Localisation in Spatial Transcriptomics",
    "size": 58, "bold": True, "color": WHITE, "line_spacing": 1.0}])
textbox(2.0, 6.5, 72, 1.7, [{
    "text": "Bivariate point-pattern statistics on raw transcripts — no cells, no clusters, no segmentation.",
    "size": 26, "italic": True, "color": META}])
textbox(2.0, 8.5, 76, 1.3, [{
    "text": "Harry Woodward   ·   Supervisor: Dr Dan Nicolau, Nicolau Lab   ·   MSc Bioinformatics, King's College London   ·   2026",
    "size": 22, "color": WHITE}])
textbox(2.0, 10.0, 76, 1.5, [{
    "text": "Platform: NanoString CosMx SMI    ·    Tissue: severe asthmatic, RSV-infected human lung    ·    ~2.66 M transcripts    ·    1,000-gene panel",
    "size": 18, "color": META}])
# KCL logo top-right
_logo = Image.open(ASSETS / "pic_48.png")
_lar = _logo.size[0] / _logo.size[1]
_lw = 8.6
shapes.add_picture(str(ASSETS / "pic_48.png"), Cm(84.1 - _lw - 1.6),
                   Cm(1.1), width=Cm(_lw))

# ================= COLUMN PANELS =================
COL_TOP, COL_BOT = 13.2, 100.6
col1, col2, col3 = 2.0, 29.2, 56.4
colw = 25.7
for cx in (col1, col2, col3):
    rect(cx, COL_TOP, colw, COL_BOT - COL_TOP, PANEL, rounded=True)

GUI = {"color": GUIDE, "italic": True, "size": 13}  # guide-run override

# ================= LEFT COLUMN =================
x, w = col1 + 0.6, colw - 1.2
figw = colw - 1.2
y = COL_TOP + 0.5

y = heading("Introduction", x, w, y)
textbox(x, y, w, 13.0, [
    {"text": "Spatial transcriptomics maps gene expression at single-molecule resolution in intact tissue. Almost every downstream analysis — cell-type calling, ligand–receptor inference, niche detection — requires a fragile first step: assigning each transcript to a cell.",
     "size": 16, "space_after": 8, "line_spacing": 1.04},
    {"text": "On the CosMx sample analysed here, vendor segmentation leaves ~42% of transcripts unassigned and most subcellular labels missing. Third-party tools (Cellpose, MOSAIK) cannot recover the signal from poor DAPI staining, so standard downstream tools (CellChat, LIANA, Bento) become unusable.",
     "size": 16, "space_after": 8, "line_spacing": 1.04},
    {"text": "NoSeggs takes a different route: treat each gene as a 2D point pattern in tissue space, and use bivariate Ripley's K-functions with label-permutation envelopes to detect ligand–receptor co-localisation directly. No cells. No clusters. Just transcripts.",
     "size": 16, "space_after": 6, "line_spacing": 1.04},
    {"runs": [("[ FILL: confirm 42% dropout wording matches M1 = 43.6%; add one sentence on why this dataset matters clinically. ]", GUI)]},
])
y += 13.4

y = heading("The Problem — Segmentation Fails", x, w, y, size=26)
y = figure(FIGS / "14_M0_raw_data.png", col1 + 0.6, figw, y)
y += 0.25
y = legend("Fig M0.",
           "One field of view, three layers. Left: raw DAPI morphology. Middle: the vendor cell-segmentation mask (148 cells) — the layer NoSeggs discards. Right: the raw transcripts, with ~40% never assigned to any cell (magenta). NoSeggs carries every point through.",
           x, w, y, 3.6)
y += 0.3
y = figure(FIGS / "14_M1_segmentation_dropout.png", col1 + 0.6, figw, y)
y += 0.25
y = legend("Fig M1.",
           "Across all FOVs, vendor segmentation never assigns 43.6% of transcripts to a cell, so every cell-level tool silently discards them — consistently across the whole dataset.",
           x, w, y, 2.9)
y += 0.3

y = heading("From Transcripts to Point Patterns", x, w, y, size=26)
y = figure(FIGS / "14_M2_gmm_strips.png", col1 + 0.6, figw, y)
y += 0.25
y = legend("Fig M2.",
           "After PCA aligns the tissue axis, a per-FOV 3-component Gaussian mixture separates the three tissue strips; boundaries are set per FOV, with manual overrides where the mixture is ambiguous.",
           x, w, y, 3.0)
y += 0.4
y = subhead("Adaptive noise removal", x, w, y)
y = figure(FIGS / "14_M3_dbscan_qc.png", col1 + 0.6, figw, y)
y += 0.25
y = legend("Fig M3.",
           "DBSCAN with an adaptively-set ε (97th-percentile 1-NN, clipped 20–30 px) flags diffuse background and stray clusters; 87% of transcripts survive QC.",
           x, w, y, 2.9)

# ================= MIDDLE COLUMN =================
x, w = col2 + 0.6, colw - 1.2
figw_big = 20.0          # compressed width for tall figures
y = COL_TOP + 0.5

y = heading("Methods — The NoSeggs Pipeline", x, w, y)
y = figure(FIGS / "14_M6_workflow.png", col2 + 0.6, colw - 1.2, y, target_w=16.5)
y += 0.2
y = legend("Fig M6.",
           "The end-to-end workflow and a quick-start guide. Eight stages from raw transcripts to signalling networks; amber chips mark the tunable options at each step (GMM components, adaptive ε, window type, multi-scale r-grid, n_sim, binary vs continuous SES).",
           x, w, y, 3.2)
y += 0.25
textbox(x, y, w, 3.2, [
    {"runs": [("Significance → effect size.  ", {"bold": True, "size": 15, "color": HEADING}),
              ("Under the production concave-hull window even the positive control sits within the permutation envelope — the window is the dominant lever for inferential power. We rank pairs by continuous standardised effect size (SES), the continuous form of the global rank envelope test (Myllymäki 2017), not binary significance.",
               {"size": 15})],
     "space_after": 5, "line_spacing": 1.03},
])
y += 3.4

y = subhead("Observation window", x, w, y)
y = figure(FIGS / "14_M4_window_iterations.png", col2 + 0.6, colw - 1.2, y, target_w=figw_big)
y += 0.2
y = legend("Fig M4.",
           "The observation window is the dominant lever for inferential power. A concave hull tracks the true tissue extent, shrinking the area ~20× vs the bounding box and tightening the permutation null.",
           x, w, y, 2.7)
y += 0.3

y = subhead("Multi-scale r-grid & edge correction", x, w, y)
y = figure(FIGS / "14_M5_rscale.png", col2 + 0.6, colw - 1.2, y)
y += 0.2
y = legend("Fig M5.",
           "A deliberate multi-scale r-grid resolves both cell-scale and tissue-scale co-association; Shapely edge correction down-weights discs that spill outside the tissue window.",
           x, w, y, 2.6)
y += 0.35

y = heading("Where NoSeggs Sits", x, w, y, size=26)
y = figure(FIGS / "14_T_tool_positioning.png", col2 + 0.6, colw - 1.2, y, target_w=figw_big)
y += 0.2
y = legend("Fig T.",
           "No published tool screens L–R co-localisation on raw transcripts without segmentation. CellChat/LIANA, Squidpy, Bento and Baysor require segmentation or cell-typing; FICTURE is segmentation-free but answers a different question; spatstat is a generic point-pattern library. NoSeggs satisfies all four capabilities.",
           x, w, y, 3.4)
y += 0.35

y = heading("Validation & Controls", x, w, y, size=26)
y = figure(FIGS / "13_pos_vs_neg_sanity.png", col2 + 0.6, colw - 1.2, y, target_w=figw_big)
y += 0.2
y = legend("Fig B.",
           "Both the positive control (KRT8×KRT18) and the housekeeping×specific negative (MALAT1×KRT18) fall inside the label-swap envelope on the infected strip, yet their effect sizes separate by 2.67 envelope half-widths. We therefore rank by continuous peak SES, not binary pass/fail.",
           x, w, y, 3.4)
y += 0.35

# validation table (native editable table)
tbl_rows = [
    ("Control", "Pair", "Expected", "Observed"),
    ("Positive", "KRT8 × KRT18", "Above null", "SES = −2.18 (top)"),
    ("Negative 1", "MALAT1 × KRT18", "Far below", "SES = −4.85"),
    ("Negative 2", "KRT8 × SCGB3A1", "Below", "[ FILL ]"),
]
tb_h = 5.6
gtbl = shapes.add_table(len(tbl_rows), 4, Cm(x), Cm(y), Cm(w), Cm(tb_h)).table
gtbl.columns[0].width = Cm(5.2)
gtbl.columns[1].width = Cm(7.6)
gtbl.columns[2].width = Cm(5.4)
gtbl.columns[3].width = Cm(w - 18.2)
for ri, row in enumerate(tbl_rows):
    for ci, val in enumerate(row):
        cell = gtbl.cell(ri, ci)
        cell.margin_left = Cm(0.2); cell.margin_right = Cm(0.1)
        cell.margin_top = Cm(0.05); cell.margin_bottom = Cm(0.05)
        cell.vertical_anchor = MSO_ANCHOR.MIDDLE
        tf = cell.text_frame
        tf.word_wrap = True
        p = tf.paragraphs[0]
        r = p.add_run(); r.text = val
        r.font.name = F
        r.font.size = Pt(14 if ri else 14)
        r.font.bold = (ri == 0)
        if ri == 0:
            r.font.color.rgb = WHITE
            cell.fill.solid(); cell.fill.fore_color.rgb = PURPLE
        else:
            r.font.color.rgb = BODY
            cell.fill.solid()
            cell.fill.fore_color.rgb = WHITE if ri % 2 else PANEL
y += tb_h + 0.4

y = figure(FIGS / "13_C_group_ses_distribution.png", col2 + 0.6, colw - 1.2, y, target_w=figw_big)
y += 0.2
y = legend("Fig C.",
           "Method is correctly null for diffuse markers (cytokeratins, housekeeping anchors) and elevated for focal pairs (paracrine, infection, stromal, immune). Bivariate K under this null detects excess spatial co-association beyond marginal density.",
           x, w, y, 3.0)

# ================= RIGHT COLUMN =================
x, w = col3 + 0.6, colw - 1.2
figw = colw - 1.2
y = COL_TOP + 0.5

y = heading("Results — Infection Biology", x, w, y)
textbox(x, y, w, 4.3, [
    {"text": "On RSV-infected lung, immune cell-communication pairs preferentially co-localise on infected tissue — antigen presentation, T-cell recruitment, antiviral monocyte signalling — recovered with no segmentation step.",
     "size": 16, "space_after": 5, "line_spacing": 1.04},
    {"runs": [("[ FILL: 1-line headline you want the reader to leave with. ]", GUI)]},
])
y += 4.5

y = subhead("Infection signalling boost", x, w, y)
y = figure(FIGS / "13_F_infection_signature.png", col3 + 0.6, figw, y)
y += 0.2
y = legend("Fig F.",
           "Per-pair Strip-2 boost for the surviving infection-response (G7) pairs. Immune cell-communication pairs are positive (5/6): antigen presentation onto T cells (HLA-DRA×CD3D, +0.74), chemokine T-cell recruitment (CXCL10×CD3D, +0.59), macrophage antigen presentation (HLA-DRA×CD68, +0.32). A paradoxical second finding only a segmentation-free method can surface: antiviral×epithelial pairs (ISG15×KRT8, MX1×KRT5) move off infected tissue, consistent with epithelial loss redistributing antiviral activity into immune cells.",
           x, w, y, 7.2)
y += 0.35

y = subhead("Co-localisation networks", x, w, y)
y = figure(FIGS / "13_G_network.png", col3 + 0.6, figw, y)
y += 0.2
y = legend("Fig G.",
           "Top-50 spatial co-association edges per strip with Louvain community detection. Communities are spatial signalling modules — interferon response, chemokine signalling, stromal ECM — recovered without any cell segmentation. Strip 2 (infected) carries more and larger modules than the controls.",
           x, w, y, 4.4)
y += 0.35

y = subhead("Spatial signatures (UMAP)", x, w, y)
y = figure(FIGS / "13_H_umap.png", col3 + 0.6, figw, y, target_w=18.5)
y += 0.2
y = legend("Fig H.",
           "Each (pair, strip) treated as a 50-dim L(r) curve-shape vector, embedded with UMAP. Stromal, immune and infection-response pairs occupy distinct regions — the method recovers biological grouping as soft visual structure from raw transcripts alone.",
           x, w, y, 4.4)

# ================= BOTTOM BAND =================
BB_TOP = 101.0
rect(0, BB_TOP, 84.1, 15.6, PURPLE)
rect(0, BB_TOP, 84.1, 0.12, META)
# Conclusions + next steps
textbox(2.0, BB_TOP + 0.5, 31, 0.9,
        [{"text": "Conclusions & Next Steps", "size": 22, "bold": True, "color": WHITE}])
textbox(2.0, BB_TOP + 1.6, 31, 9.0, [
    {"text": "Bivariate Ripley's K on raw transcript point clouds detects L–R co-localisation with no cell segmentation step.",
     "size": 14.5, "color": WHITE, "bullet": True, "space_after": 4, "line_spacing": 1.02},
    {"text": "Discriminates positive from negative controls by effect size (separation = +2.67 envelope half-widths under the production window).",
     "size": 14.5, "color": WHITE, "bullet": True, "space_after": 4, "line_spacing": 1.02},
    {"text": "Fills a vacant niche: Bento, Squidpy, CellChat, LIANA all require segmentation — NoSeggs does not.",
     "size": 14.5, "color": WHITE, "bullet": True, "space_after": 9, "line_spacing": 1.02},
    {"runs": [("Next: ", {"bold": True, "size": 14.5, "color": META}),
              ("extended cytokeratin/housekeeping anchor panel; marker-gene pseudo-segmentation for cell-type stratification; benchmark vs Squidpy gr.ripley on cell-level reductions.",
               {"size": 14.5, "color": WHITE})],
     "line_spacing": 1.02},
])
# Key references
textbox(35.5, BB_TOP + 0.5, 30, 0.9,
        [{"text": "Key References", "size": 22, "bold": True, "color": WHITE}])
refs = [
    "Ripley BD (1976) J. Appl. Prob. 13:255–266 — the K-function.",
    "Baddeley A, Rubak E, Turner R (2015) Spatial Point Patterns. CRC Press.",
    "Mah CK et al. (2024) Genome Biology 25:82 — Bento.",
    "Palla G et al. (2022) Nature Methods 19:171–178 — Squidpy.",
    "Myllymäki M et al. (2017) J. R. Stat. Soc. B 79:381–404 — rank envelope test.",
]
textbox(35.5, BB_TOP + 1.6, 30, 9.0,
        [{"runs": [(f"{i+1}. ", {"bold": True, "size": 13, "color": META}),
                   (r, {"size": 13, "color": WHITE})],
          "space_after": 4, "line_spacing": 1.02}
         for i, r in enumerate(refs)] +
        [{"runs": [("[ FILL: confirm citations / add a 6th if room. ]", GUI)]}])
# Acknowledgements
textbox(67.0, BB_TOP + 0.5, 13.5, 0.9,
        [{"text": "Acknowledgements", "size": 22, "bold": True, "color": WHITE}])
textbox(67.0, BB_TOP + 1.6, 13.0, 7.5, [
    {"text": "Dataset and HPC compute courtesy of King's College London (CREATE HPC, msc_appbio). RSV samples via [TBD]. CellChatDB via Jin et al. 2021. Supervised by Dr Dan Nicolau, Nicolau Lab.",
     "size": 13, "color": WHITE, "line_spacing": 1.04, "space_after": 6},
    {"runs": [("[ FILL: dataset PI / funder attribution. ]", GUI)]},
])
# QR code bottom-right
shapes.add_picture(str(ASSETS / "pic_50.png"), Cm(76.6), Cm(BB_TOP + 9.6),
                   width=Cm(5.4))

# ================= FOOTER BAR =================
rect(0, 116.9, 84.1, 2.0, RGBColor(0x3A, 0x05, 0x57))
textbox(2.0, 117.15, 80, 1.4,
        [{"text": "github.com/harrydaw/cosmx-pointpattern    ·    harrywoodward99@gmail.com    ·    MSc Bioinformatics — King's College London — 2026",
          "size": 14, "color": META, "align": PP_ALIGN.CENTER}])

prs.save(str(OUT))
print("saved", OUT)
print("slide cm:", round(prs.slide_width/360000, 2), round(prs.slide_height/360000, 2))
