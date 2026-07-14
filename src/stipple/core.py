"""STIPPLE core library: segmentation-free ligand-receptor co-localisation.

Bivariate, edge-corrected Ripley's K/L analysis on raw transcript point patterns,
with a label-permutation envelope and a continuous standardised effect size (SES).

Extracted from notebooks/00_functions.ipynb as an importable module. The public
entry points are re-exported in ``stipple/__init__.py``; see the README quick start.
"""

import json
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from shapely.geometry import MultiPoint, Point, Polygon as ShapelyPolygon
import shapely
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
import matplotlib.cm as cm


def add_rotated_coords(df, tissue_angle_deg=None,
                       x_col='x_global_px', y_col='y_global_px'):
    """
    Add x_rot_px / y_rot_px columns aligned so the cross-strip axis runs along x_rot.

    By default (tissue_angle_deg=None) the rotation angle is auto-detected via PCA:
    the first principal component (maximum-variance axis, i.e. along-strip direction)
    is rotated to be vertical, so strips separate cleanly along the rotated x-axis.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain x_col and y_col.
    tissue_angle_deg : float or None
        Anticlockwise rotation in degrees to apply.
        None (default) -> auto-detect using PCA on all transcript coordinates.
    x_col, y_col : str
        Source coordinate columns. Defaults: x_global_px, y_global_px.

    Returns
    -------
    pd.DataFrame
        Copy of df with x_rot_px and y_rot_px added.
        Original coordinate columns are untouched.

    Diagnostic plot
    ---------------
    Left panel  - original coordinates with black (original) and red (rotated) axes
    Right panel - rotated coordinates; strips should now separate along x_rot_px
    """
    from sklearn.decomposition import PCA

    coords   = df[[x_col, y_col]].values.astype(np.float64)
    centroid = coords.mean(axis=0)
    coords_c = coords - centroid

    if tissue_angle_deg is None:
        pca = PCA(n_components=2)
        pca.fit(coords_c)
        # PC1 = long axis (along-strip direction)
        # Rotate so PC1 becomes vertical -> strips separate along rotated x
        pc1       = pca.components_[0]
        angle_pc1 = np.degrees(np.arctan2(pc1[1], pc1[0]))
        tissue_angle_deg = 90.0 - angle_pc1
        tissue_angle_deg = ((tissue_angle_deg + 180) % 360) - 180  # normalise to [-180, 180]
        print(f"PCA auto-detected rotation: {tissue_angle_deg:.1f} deg anticlockwise")
    else:
        print(f"Using manual rotation: {tissue_angle_deg:.1f} deg anticlockwise")

    theta = np.radians(tissue_angle_deg)
    R     = np.array([[np.cos(theta), -np.sin(theta)],
                      [np.sin(theta),  np.cos(theta)]])
    rotated = coords_c @ R.T + centroid

    result = df.copy()
    result['x_rot_px'] = rotated[:, 0]
    result['y_rot_px'] = rotated[:, 1]

    # ── Diagnostic plot ────────────────────────────────────────────────────
    rng      = np.random.default_rng(0)
    n_sample = min(50_000, len(coords))
    idx      = rng.choice(len(coords), n_sample, replace=False)

    x_range  = coords[:, 0].max() - coords[:, 0].min()
    y_range  = coords[:, 1].max() - coords[:, 1].min()
    ax_len   = max(x_range, y_range) * 0.12

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

    # Left: original coords + axes
    ax1.scatter(coords[idx, 0], coords[idx, 1],
                s=0.3, alpha=0.25, c='steelblue', rasterized=True)
    cx, cy = centroid
    # Original axes (black)
    ax1.annotate('', xy=(cx + ax_len, cy), xytext=(cx, cy),
                 arrowprops=dict(arrowstyle='->', color='black', lw=2))
    ax1.annotate('', xy=(cx, cy + ax_len), xytext=(cx, cy),
                 arrowprops=dict(arrowstyle='->', color='black', lw=2))
    ax1.text(cx + ax_len * 1.05, cy, 'x', ha='left', va='center', fontsize=9)
    ax1.text(cx, cy + ax_len * 1.05, 'y', ha='center', va='bottom', fontsize=9)
    # Rotated axes (red)
    ax1.annotate('', xy=(cx + ax_len * np.cos(theta),
                         cy + ax_len * np.sin(theta)), xytext=(cx, cy),
                 arrowprops=dict(arrowstyle='->', color='red', lw=2))
    ax1.annotate('', xy=(cx + ax_len * np.cos(theta + np.pi/2),
                         cy + ax_len * np.sin(theta + np.pi/2)), xytext=(cx, cy),
                 arrowprops=dict(arrowstyle='->', color='red', lw=2))
    ax1.text(cx + ax_len * np.cos(theta) * 1.15,
             cy + ax_len * np.sin(theta) * 1.15,
             'x_rot', ha='center', va='center', fontsize=9, color='red')
    ax1.set_title(
        f'Original coordinates\n'
        f'Black = original axes | Red = rotated axes ({tissue_angle_deg:.1f}\u00b0 ACW)',
        fontsize=10)
    ax1.set_aspect('equal')
    ax1.set_xlabel(x_col, fontsize=8)
    ax1.set_ylabel(y_col, fontsize=8)

    # Right: rotated coords
    ax2.scatter(rotated[idx, 0], rotated[idx, 1],
                s=0.3, alpha=0.25, c='steelblue', rasterized=True)
    ax2.set_title(
        'Rotated coordinates (x_rot_px, y_rot_px)\n'
        'Strips should now separate cleanly along x_rot_px',
        fontsize=10)
    ax2.set_aspect('equal')
    ax2.set_xlabel('x_rot_px', fontsize=8)
    ax2.set_ylabel('y_rot_px', fontsize=8)

    plt.tight_layout()
    plt.show()

    return result


def reassign_strips_gmm(df, x_col='x_rot_px', strip_col='strip',
                        fov_n_strips=None, random_state=42,
                        clean_only=True):
    """
    Re-run per-FOV GMM strip assignment on a rotated coordinate column.

    By default the GMM is *fitted* only on clean transcripts (those with
    is_noise=False and manually_excluded=False, when those columns exist),
    so diffuse DBSCAN noise and auto-excluded small clusters do not
    distort the component means. Strip labels are *predicted* for every
    row so downstream filters remain straightforward.

    Parameters
    ----------
    df : pd.DataFrame  — must contain 'fov' and x_col
    x_col : str        — column to fit GMM on (default 'x_rot_px')
    strip_col : str    — column to write labels into (default 'strip')
    fov_n_strips : dict or None
        {fov_id: tuple_of_strip_names}
        Override the strips present in specific FOVs. Number of GMM
        components = len(tuple). Components are ordered by ascending
        mean and mapped to the supplied strip names in that order.
        Examples:
          {4: ('strip_1', 'strip_2')}        # left two strips only
          {12: ('strip_1', 'strip_3')}       # strip_2 absent
          {7: ('strip_2',)}                  # single strip
        FOVs not listed default to ('strip_1', 'strip_2', 'strip_3').
    random_state : int
    clean_only : bool

    Returns
    -------
    result : pd.DataFrame
        Copy of df with strip_col reassigned.
    boundaries : dict
        {fov_id: [b1, b2, ...]}
        Positional boundary values between adjacent strips in
        fov_n_strips[fov]. Length = len(strips) - 1. Empty list for
        single-strip FOVs. Each value is a float (midpoint of adjacent
        GMM component means) or None if the fit was skipped.
    """
    from sklearn.mixture import GaussianMixture

    if fov_n_strips is None:
        fov_n_strips = {}

    result     = df.copy()
    fovs       = sorted(df['fov'].unique())
    boundaries = {}

    clean_mask = pd.Series(True, index=result.index)
    if clean_only:
        if 'is_noise' in result.columns:
            clean_mask &= ~result['is_noise']
        if 'manually_excluded' in result.columns:
            clean_mask &= ~result['manually_excluded']

    n_fit = clean_mask.sum()
    print(f"Reassigning strips per-FOV on '{x_col}' "
          f"(fitting on {n_fit:,} clean / {len(result):,} total rows):")

    for fov_id in fovs:
        fov_mask    = result['fov'] == fov_id
        fit_mask    = fov_mask & clean_mask
        strip_names = fov_n_strips.get(fov_id, ('strip_1', 'strip_2', 'strip_3'))
        n_comp      = len(strip_names)

        if fit_mask.sum() < n_comp:
            print(f"  FOV {fov_id:2d}: only {fit_mask.sum()} clean points "
                  f"(< {n_comp}) — skipping GMM, leaving labels untouched")
            boundaries[fov_id] = [None] * (n_comp - 1)
            continue

        if n_comp == 1:
            result.loc[fov_mask, strip_col] = strip_names[0]
            boundaries[fov_id] = []
            print(f"  FOV {fov_id:2d}: single strip → all assigned to {strip_names[0]}")
            continue

        X_fit = result.loc[fit_mask, x_col].values.reshape(-1, 1).astype(np.float64)
        gmm   = GaussianMixture(n_components=n_comp, random_state=random_state)
        gmm.fit(X_fit)

        order     = np.argsort(gmm.means_[:, 0])
        label_map = {old: strip_names[new] for new, old in enumerate(order)}

        X_all       = result.loc[fov_mask, x_col].values.reshape(-1, 1).astype(np.float64)
        raw_labels  = gmm.predict(X_all)
        new_strip   = np.array([label_map[l] for l in raw_labels])
        result.loc[fov_mask, strip_col] = new_strip

        means  = gmm.means_[order, 0]
        # Positional boundaries: midpoint between adjacent component means
        bvals  = [float((means[i] + means[i + 1]) / 2) for i in range(n_comp - 1)]
        boundaries[fov_id] = bvals

        clean_strip = result.loc[fit_mask, strip_col]
        counts      = {s: int((clean_strip == s).sum()) for s in strip_names}
        mean_str    = '[' + ', '.join(f'{m:.0f}' for m in means) + ']'
        b_str       = '[' + ', '.join(f'{b:.0f}' for b in bvals) + ']'
        print(f"  FOV {fov_id:2d}: strips={strip_names}  means={mean_str}  "
              f"boundaries={b_str}  "
              + '  '.join(f'{s}={counts[s]:,}' for s in strip_names))

    return result, boundaries


def get_coords(df: pd.DataFrame,
               gene: str,
               x_col: str = 'x_global_px',
               y_col: str = 'y_global_px') -> np.ndarray:
    """Extract (N, 2) coordinate array for a given gene. Returns (0, 2) if absent."""
    coords = df.loc[df['target'] == gene, [x_col, y_col]].values
    return coords.astype(np.float64)


def k_to_l(k_vals: np.ndarray, r_vals: np.ndarray) -> np.ndarray:
    """Variance-stabilising transform: L(r) = sqrt(K(r)/pi) - r. Under CSR, L(r) = 0."""
    return np.sqrt(np.maximum(k_vals, 0) / np.pi) - r_vals


def get_convex_hull(df: pd.DataFrame,
                    x_col: str = 'x_global_px',
                    y_col: str = 'y_global_px'):
    """
    Compute the convex hull of all transcript coordinates in a strip DataFrame.

    Returns
    -------
    shapely.geometry.Polygon

    Raises
    ------
    ValueError
        If fewer than 3 points are present or all points are collinear.
    """
    coords = df[[x_col, y_col]].values.astype(np.float64)
    if len(coords) < 3:
        raise ValueError(
            f'Cannot compute convex hull: strip contains {len(coords)} point(s), '
            'minimum 3 required.'
        )
    hull = MultiPoint(coords).convex_hull
    if hull.geom_type != 'Polygon':
        raise ValueError(
            f'Convex hull is degenerate ({hull.geom_type}): all {len(coords)} '
            'points are collinear.'
        )
    return hull


def get_concave_hull(df: pd.DataFrame,
                     x_col: str = 'x_global_px',
                     y_col: str = 'y_global_px',
                     concave_ratio: float = 0.1):
    """
    Compute a concave hull of transcript coordinates using Shapely 2.x.

    Parameters
    ----------
    df : pd.DataFrame
        Strip- or (FOV, strip)-level transcript DataFrame.
    x_col, y_col : str
        Coordinate column names.
    concave_ratio : float in [0, 1]
        Controls hull tightness. 0 = as tight as possible (hugs data edge),
        1 = identical to convex hull. Default 0.1.

    Returns
    -------
    shapely.geometry.Polygon

    Notes
    -----
    Uses shapely.concave_hull() (Shapely >= 2.0 required).
    allow_holes=False prevents the hull from cutting holes through tissue.
    Falls back to convex hull if the concave result is degenerate
    (e.g. MultiPolygon or empty for very sparse data).
    """
    coords = df[[x_col, y_col]].values.astype(np.float64)
    if len(coords) < 3:
        raise ValueError(
            f'Cannot compute hull: strip contains {len(coords)} point(s), '
            'minimum 3 required.'
        )
    mp = MultiPoint(coords)
    hull = shapely.concave_hull(mp, ratio=concave_ratio, allow_holes=False)
    # Fall back to convex hull if concave result is not a simple Polygon
    if hull.geom_type != 'Polygon' or hull.is_empty:
        hull = mp.convex_hull
    if hull.geom_type != 'Polygon':
        raise ValueError(
            f'Hull is degenerate ({hull.geom_type}): all points may be collinear.'
        )
    return hull


def fraction_inside_hull(point, r, hull, resolution: int = 64):
    """
    Fraction of a disc of radius r centred at point that lies inside hull.

    Used as the per-point edge-correction weight in the polygon path of
    bivariate_k(). Approximated via Shapely polygon intersection.

    Parameters
    ----------
    point : tuple of float — (x, y)
    r : float
    hull : shapely.geometry.Polygon
    resolution : int
        Segments per quarter-circle for the disc approximation. Default 64.

    Returns
    -------
    float in (0, 1]
    """
    disc = Point(point).buffer(r, resolution=resolution)
    if hull.contains(disc):
        return 1.0
    intersection = hull.intersection(disc)
    return intersection.area / disc.area


def sample_in_polygon(polygon, n: int, rng) -> np.ndarray:
    """Sample n points uniformly inside a convex polygon via rejection sampling."""
    minx, miny, maxx, maxy = polygon.bounds
    pts = []
    while len(pts) < n:
        candidates = np.column_stack([
            rng.uniform(minx, maxx, n * 2),
            rng.uniform(miny, maxy, n * 2),
        ])
        mask = shapely.contains_xy(polygon, candidates[:, 0], candidates[:, 1])
        pts.extend(candidates[mask].tolist())
    return np.array(pts[:n])


def load_custom_window(path: str):
    """
    Load a user-drawn polygon from a JSON file saved by draw_custom_window() in 09b.

    Parameters
    ----------
    path : str
        Path to JSON file with key 'vertices': list of [x, y] pairs.

    Returns
    -------
    shapely.geometry.Polygon
    """
    with open(path) as f:
        data = json.load(f)
    poly = ShapelyPolygon(data['vertices'])
    if not poly.is_valid:
        raise ValueError(f'Loaded polygon from {path} is not valid: {poly.geom_type}')
    return poly


def get_window(df: pd.DataFrame,
               window_type: str = 'concave',
               custom_path: str = None,
               concave_ratio: float = 0.1,
               x_col: str = 'x_global_px',
               y_col: str = 'y_global_px'):
    """
    Return the observation window for a strip or (FOV, strip) DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
    window_type : str
        'rect'    — (x_min, x_max, y_min, y_max) tuple
        'hull'    — convex hull Shapely Polygon
        'concave' — concave hull Shapely Polygon (default)
        'custom'  — polygon loaded from custom_path JSON (drawn in 09b)
    custom_path : str or None
        Required when window_type='custom'.
    concave_ratio : float
        Passed to get_concave_hull(). Ignored for other window types.
        0 = tightest fit, 1 = convex hull. Default 0.1.
    x_col, y_col : str

    Returns
    -------
    tuple or shapely.geometry.Polygon
    """
    if window_type == 'rect':
        return (
            df[x_col].min(), df[x_col].max(),
            df[y_col].min(), df[y_col].max(),
        )
    elif window_type == 'hull':
        return get_convex_hull(df, x_col, y_col)
    elif window_type == 'concave':
        return get_concave_hull(df, x_col, y_col, concave_ratio=concave_ratio)
    elif window_type == 'custom':
        if custom_path is None:
            raise ValueError("custom_path must be provided when window_type='custom'")
        return load_custom_window(custom_path)
    else:
        raise ValueError(
            f"Unknown window_type: {window_type!r}. "
            "Use 'rect', 'hull', 'concave', or 'custom'."
        )


def bivariate_k(coords_a: np.ndarray,
                coords_b: np.ndarray,
                r_vals: np.ndarray,
                window,
                resolution: int = 64,
                max_n: int | None = None,
                seed: int | None = None) -> np.ndarray:
    """
    Bivariate (cross-type) Ripley's K-function K_AB(r).

    Dispatches on window type:
      tuple           -> rectangular isotropic arc-fraction edge correction
      Shapely Polygon -> Shapely disc-intersection edge correction

    Under CSR, K_AB(r) = pi * r^2.

    Parameters
    ----------
    coords_a : np.ndarray, shape (n_a, 2)
    coords_b : np.ndarray, shape (n_b, 2)
    r_vals : np.ndarray, shape (n_r,)
    window : tuple or shapely.geometry.Polygon
    resolution : int
        Disc approximation resolution for the polygon path. Ignored for tuple windows.
    max_n : int or None
        If set, subsample coords_a and coords_b independently to at most
        max_n points before computing K. Caps the per-call cost on very
        abundant genes (e.g. MALAT1). None (default) = no subsampling.
    seed : int or None
        RNG seed for the subsample. Only used when max_n triggers a
        subsample. For a like-for-like observed vs envelope comparison,
        pass the same max_n and seed to compute_envelope.

    Returns
    -------
    np.ndarray, shape (n_r,)
    """
    if max_n is not None:
        _rng = np.random.default_rng(seed)
        if len(coords_a) > max_n:
            coords_a = coords_a[_rng.choice(len(coords_a), max_n, replace=False)]
        if len(coords_b) > max_n:
            coords_b = coords_b[_rng.choice(len(coords_b), max_n, replace=False)]

    n_a = len(coords_a)
    n_b = len(coords_b)

    if isinstance(window, tuple):
        x_min, x_max, y_min, y_max = window
        area = (x_max - x_min) * (y_max - y_min)
    else:
        area = window.area
    lambda_b = n_b / area

    tree_b = cKDTree(coords_b)
    r_max = r_vals[-1]
    neighbours = tree_b.query_ball_point(coords_a, r_max)

    k_vals = np.zeros(len(r_vals))

    for idx_a in range(n_a):
        x_i, y_i = coords_a[idx_a]
        b_indices = neighbours[idx_a]
        if len(b_indices) == 0:
            continue

        b_pts = coords_b[b_indices]
        dists = np.sqrt(((b_pts - coords_a[idx_a]) ** 2).sum(axis=1))

        if isinstance(window, tuple):
            dx_left  = x_i - x_min
            dx_right = x_max - x_i
            dy_bot   = y_i - y_min
            dy_top   = y_max - y_i

            for d_ij in dists:
                if d_ij == 0:
                    continue
                angle_outside = 0.0
                for dist_to_edge in [dx_left, dx_right, dy_bot, dy_top]:
                    if dist_to_edge < d_ij:
                        angle_outside += 2.0 * np.arccos(
                            np.clip(dist_to_edge / d_ij, -1.0, 1.0)
                        )
                frac = max((2.0 * np.pi - angle_outside) / (2.0 * np.pi), 0.01)
                weight = 1.0 / frac
                for k, r in enumerate(r_vals):
                    if d_ij <= r:
                        k_vals[k] += weight

        else:
            for d_ij in dists:
                if d_ij == 0:
                    continue
                frac = fraction_inside_hull(
                    (x_i, y_i), d_ij, window, resolution=resolution
                )
                frac = max(frac, 0.01)
                weight = 1.0 / frac
                for k, r in enumerate(r_vals):
                    if d_ij <= r:
                        k_vals[k] += weight

    k_vals /= (n_a * lambda_b)
    return k_vals


def compute_envelope(coords_a: np.ndarray,
                     coords_b: np.ndarray,
                     r_vals: np.ndarray,
                     window,
                     n_sim: int = 99,
                     seed: int = 42,
                     resolution: int = 64,
                     max_n: int | None = None) -> tuple:
    """
    Monte Carlo permutation envelope for the bivariate L-function.

    Null model: pool all transcript locations and permute gene labels.
    Window (tuple or Polygon) is held fixed.

    Parameters
    ----------
    window : tuple or shapely.geometry.Polygon
    n_sim : int  - 99 gives pointwise alpha ~0.02.
    resolution : int - disc resolution for polygon path.
    max_n : int or None
        If set, subsample coords_a and coords_b (each capped to max_n)
        once at the start using `seed`, then build the permutation pool
        from the subsamples. For a like-for-like comparison, the observed
        call to bivariate_k must use the same max_n and seed.

    Returns
    -------
    l_lo, l_hi : np.ndarray, shape (n_r,)
    """
    rng = np.random.default_rng(seed)

    if max_n is not None:
        if len(coords_a) > max_n:
            coords_a = coords_a[rng.choice(len(coords_a), max_n, replace=False)]
        if len(coords_b) > max_n:
            coords_b = coords_b[rng.choice(len(coords_b), max_n, replace=False)]

    n_a = len(coords_a)
    pooled = np.vstack([coords_a, coords_b])
    sim_l = np.zeros((n_sim, len(r_vals)))

    for s in range(n_sim):
        idx = rng.permutation(len(pooled))
        sim_a = pooled[idx[:n_a]]
        sim_b = pooled[idx[n_a:]]
        k_sim = bivariate_k(sim_a, sim_b, r_vals, window, resolution=resolution)
        sim_l[s] = k_to_l(k_sim, r_vals)

    return sim_l.min(axis=0), sim_l.max(axis=0)


def plot_diagnostics(coords_a: np.ndarray,
                     coords_b: np.ndarray,
                     window,
                     gene_a: str = 'Gene A',
                     gene_b: str = 'Gene B',
                     title: str = '',
                     ax=None) -> None:
    """
    Scatter plot of two transcript patterns with window overlay.
    Handles both tuple (rectangular) and Shapely Polygon windows.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(coords_a[:, 0], coords_a[:, 1],
               s=8, alpha=0.6, color='steelblue',
               label=f'{gene_a} (n={len(coords_a)})')
    ax.scatter(coords_b[:, 0], coords_b[:, 1],
               s=8, alpha=0.6, color='tomato',
               label=f'{gene_b} (n={len(coords_b)})')

    if isinstance(window, tuple):
        x_min, x_max, y_min, y_max = window
        rect = plt.Rectangle(
            (x_min, y_min), x_max - x_min, y_max - y_min,
            linewidth=1.5, edgecolor='black', facecolor='none',
            linestyle='--', label='Window (rect)'
        )
        ax.add_patch(rect)
    else:
        hx, hy = window.exterior.xy
        ax.plot(hx, hy, color='black', lw=1.5, linestyle='--', label='Window')

    ax.set_xlabel('x_global_px')
    ax.set_ylabel('y_global_px')
    ax.set_title(title or f'{gene_a} × {gene_b}')
    ax.legend(fontsize=7, markerscale=2, loc='lower right', framealpha=0.9)
    ax.set_aspect('equal')


def run_pair_analysis(strip_df: pd.DataFrame,
                      gene_a: str,
                      gene_b: str,
                      r_vals: np.ndarray,
                      window_type: str = 'hull',
                      custom_path: str = None,
                      concave_ratio: float = 0.1,
                      n_sim: int = 99,
                      seed: int = 42,
                      resolution: int = 64,
                      diagnostics: bool = False,
                      strip_name: str = '') -> dict:
    """
    Full bivariate co-localisation analysis for one gene pair in one strip.

    Single entry point for all window types.

    Parameters
    ----------
    strip_df : pd.DataFrame
        Transcript DataFrame for a single strip (all genes).
    gene_a, gene_b : str
    r_vals : np.ndarray
    window_type : str
        'rect', 'hull', 'concave', or 'custom'. Default 'hull'.
    custom_path : str or None
        Required when window_type='custom'.
    concave_ratio : float
        Passed to get_window() when window_type='concave'. Ignored otherwise.
        0 = tightest fit, 1 = convex hull. Default 0.1.
    n_sim, seed, resolution : int
    diagnostics : bool
    strip_name : str

    Returns
    -------
    dict with keys:
        'coords_a', 'coords_b', 'window', 'k_obs', 'l_obs', 'l_lo', 'l_hi',
        'n_a', 'n_b'
    """
    coords_a = get_coords(strip_df, gene_a)
    coords_b = get_coords(strip_df, gene_b)
    window   = get_window(strip_df,
                          window_type=window_type,
                          custom_path=custom_path,
                          concave_ratio=concave_ratio)

    n_a, n_b = len(coords_a), len(coords_b)

    if diagnostics:
        fig, ax = plt.subplots(figsize=(8, 6))
        plot_diagnostics(coords_a, coords_b, window,
                         gene_a=gene_a, gene_b=gene_b,
                         title=f'{gene_a} × {gene_b} | {strip_name}',
                         ax=ax)
        plt.show()

    print(f'  {gene_a}: {n_a} transcripts')
    print(f'  {gene_b}: {n_b} transcripts')
    if isinstance(window, tuple):
        x_min, x_max, y_min, y_max = window
        narrow = min(x_max - x_min, y_max - y_min)
        print(f'  Window: rect, narrowest dim = {narrow:.0f} px')
    else:
        print(f'  Window: {window_type}, area = {window.area:.0f} px²')
    print(f'  R_MAX = {r_vals[-1]:.0f} px')

    t0 = time.time()
    k_obs = bivariate_k(coords_a, coords_b, r_vals, window, resolution=resolution)
    l_obs = k_to_l(k_obs, r_vals)
    print(f'  Observed L(r) in {time.time() - t0:.1f}s')

    t0 = time.time()
    l_lo, l_hi = compute_envelope(
        coords_a, coords_b, r_vals, window,
        n_sim=n_sim, seed=seed, resolution=resolution
    )
    print(f'  Envelope in {time.time() - t0:.1f}s')

    return {
        'coords_a': coords_a,
        'coords_b': coords_b,
        'window':   window,
        'k_obs':    k_obs,
        'l_obs':    l_obs,
        'l_lo':     l_lo,
        'l_hi':     l_hi,
        'n_a':      n_a,
        'n_b':      n_b,
    }


def plot_strip_assignments(df, fovs=None,
                           x_col='x_global_px', y_col='y_global_px',
                           strip_col='strip', boundaries=None,
                           fov_n_strips=None, clean_only=True):
    """
    Diagnostic plot of strip assignments per FOV.

    For each FOV, shows two panels:
      Left  — transcript scatter coloured by strip
      Right — x-coordinate histogram per strip with boundary line(s)

    ``boundaries[fov_id]`` is a positional list aligned with
    ``fov_n_strips[fov_id]`` (or the default 3-strip tuple). Entry ``i``
    is the boundary between ``strip_names[i]`` and ``strip_names[i+1]``
    and may be:
      * ``None``                 — no line drawn (falls back to midpoint
                                   of adjacent strip means if both are present)
      * ``number``               — vertical dashed line at x = number
      * ``[(x1, y1), (x2, y2)]`` — diagonal line drawn on scatter panel;
                                   summarised as a vertical at midpoint x
                                   on the histogram panel

    By default, noise and manually-excluded transcripts are filtered out
    before plotting.

    Returns
    -------
    pd.DataFrame  — boundary summary per FOV
    """
    strip_colours = {'strip_1': '#e41a1c', 'strip_2': '#377eb8', 'strip_3': '#4daf4a'}
    strip_order   = ['strip_1', 'strip_2', 'strip_3']

    if fov_n_strips is None:
        fov_n_strips = {}

    plot_df = df
    if clean_only:
        mask = pd.Series(True, index=plot_df.index)
        if 'is_noise' in plot_df.columns:
            mask &= ~plot_df['is_noise']
        if 'manually_excluded' in plot_df.columns:
            mask &= ~plot_df['manually_excluded']
        plot_df = plot_df[mask]

    if fovs is None:
        fovs = sorted(plot_df['fov'].unique())

    n_fovs = len(fovs)
    fig, axes = plt.subplots(n_fovs, 2,
                             figsize=(14, 4.5 * n_fovs),
                             squeeze=False)
    boundary_records = []

    def _describe(b):
        if b is None:                    return None
        if isinstance(b, (int, float)):  return round(float(b))
        (x1, y1), (x2, y2) = b
        return f'({x1:.0f},{y1:.0f})->({x2:.0f},{y2:.0f})'

    for row_idx, fov_id in enumerate(fovs):
        fov_df = plot_df[plot_df['fov'] == fov_id]
        ax_sc  = axes[row_idx, 0]
        ax_hi  = axes[row_idx, 1]

        strip_names = fov_n_strips.get(fov_id, ('strip_1', 'strip_2', 'strip_3'))
        bvals       = list((boundaries or {}).get(fov_id, []))

        # Scatter
        for strip in strip_order:
            sub = fov_df[fov_df[strip_col] == strip]
            ax_sc.scatter(sub[x_col], sub[y_col],
                          s=0.5, alpha=0.4, c=strip_colours[strip],
                          label=strip, rasterized=True)
        ax_sc.set_title(f'FOV {fov_id}  strips={strip_names}', fontsize=10)
        ax_sc.set_xlabel(x_col, fontsize=8)
        ax_sc.set_ylabel(y_col, fontsize=8)
        ax_sc.set_aspect('equal')
        ax_sc.legend(markerscale=8, fontsize=7, loc='upper left')

        # Histogram
        x_means = {}
        for strip in strip_order:
            sub = fov_df[fov_df[strip_col] == strip][x_col]
            if len(sub) == 0:
                continue
            ax_hi.hist(sub, bins=80, density=True, alpha=0.5,
                       color=strip_colours[strip], label=strip)
            x_means[strip] = sub.mean()

        y_min, y_max = fov_df[y_col].min(), fov_df[y_col].max()

        rec = {'fov': fov_id, 'strips': strip_names}

        # Iterate through adjacent pairs and matched boundaries
        for i in range(len(strip_names) - 1):
            s1, s2 = strip_names[i], strip_names[i + 1]
            bkey   = f'b_{s1[-1]}{s2[-1]}'   # e.g. 'b_12', 'b_23', 'b_13'
            gkey   = f'gap_{s1[-1]}{s2[-1]}'

            b = bvals[i] if i < len(bvals) else None

            if (s1 not in x_means or s2 not in x_means) and b is None:
                continue
            if b is None:
                b = (x_means[s1] + x_means[s2]) / 2

            if s1 in x_means and s2 in x_means:
                gap = x_means[s2] - x_means[s1]
                rec[gkey] = round(gap)
                colour = 'black' if gap > 0 else 'red'
            else:
                colour = 'black'

            if isinstance(b, (int, float)):
                bx = float(b)
                rec[bkey] = round(bx)
                ax_hi.axvline(bx, color=colour, lw=1.5, ls='--',
                              label=f'{bkey} ({bx:.0f} px)')
                ax_sc.axvline(bx, color=colour, lw=1, ls='--', alpha=0.7)
            else:
                (x1, y1), (x2, y2) = b
                ax_sc.plot([x1, x2], [y1, y2],
                           color=colour, lw=1.5, ls='--', alpha=0.9,
                           label=f'{bkey} line')
                if y2 != y1:
                    t_lo = (y_min - y1) / (y2 - y1)
                    t_hi = (y_max - y1) / (y2 - y1)
                    xs = [x1 + t_lo * (x2 - x1), x1 + t_hi * (x2 - x1)]
                    ys = [y_min, y_max]
                    ax_sc.plot(xs, ys, color=colour, lw=0.8, ls=':', alpha=0.5)
                midx = (x1 + x2) / 2
                rec[bkey] = _describe(b)
                ax_hi.axvline(midx, color=colour, lw=1.5, ls='--',
                              label=f'{bkey} (line, midx={midx:.0f})')

        boundary_records.append(rec)
        ax_hi.set_xlabel(f'{x_col}  (px)', fontsize=8)
        ax_hi.set_ylabel('density', fontsize=8)
        ax_hi.set_title(f'FOV {fov_id} — x distribution', fontsize=10)
        ax_hi.legend(fontsize=6)

    plt.tight_layout()
    plt.show()

    summary = pd.DataFrame(boundary_records)
    print('\nStrip boundary summary (negative gap = strip means reversed / collapsed):')
    print(summary.to_string(index=False))
    return summary


def apply_strip_overrides(df, overrides, fov_n_strips=None,
                          gmm_boundaries=None,
                          x_col='x_rot_px', y_col='y_rot_px',
                          strip_col='strip'):
    """
    Override strip assignments for specific FOVs using positional
    boundary lists aligned with fov_n_strips.

    For each FOV in ``overrides``, ``overrides[fov]`` is a list of
    boundary values, one per adjacent strip pair in
    ``fov_n_strips[fov]``. List length must equal
    ``len(strip_names) - 1``.

    Each boundary value can be:
      * ``None``                 — fall back to ``gmm_boundaries[fov][i]``
      * ``number``               — vertical cutoff at x = number
      * ``[(x1, y1), (x2, y2)]`` — diagonal line through two points;
                                   boundary x is linearly interpolated in y.

    Strip labels for all transcripts in the FOV are rewritten by counting
    how many boundaries each point lies to the right of. Boundaries are
    assumed monotone in x (ascending); if the user supplies them out of
    order the resulting slot assignment is undefined.

    Parameters
    ----------
    df : pd.DataFrame
    overrides : dict
        ``{fov_id: [b1, b2, ...]}``
    fov_n_strips : dict or None
        Same dict passed to ``reassign_strips_gmm``. FOVs absent default
        to ``('strip_1', 'strip_2', 'strip_3')``.
    gmm_boundaries : dict or None
        Output of ``reassign_strips_gmm`` (positional lists). Used to
        resolve ``None`` entries in ``overrides``. Required if any
        override entry is ``None``.
    x_col, y_col, strip_col : str

    Examples
    --------
    FOV_N_STRIPS = {12: ('strip_1', 'strip_3')}

    STRIP_OVERRIDES = {
        5:  [11200, 12400],                           # 3-strip, two verticals
        11: [None, [(12700, 0), (13200, -3800)]],     # 3-strip, 1|2 kept as GMM
        12: [12900],                                  # 2-strip (1,3), one boundary
    }
    df = apply_strip_overrides(df, STRIP_OVERRIDES,
                               fov_n_strips=FOV_N_STRIPS,
                               gmm_boundaries=gmm_boundaries)
    """
    if not overrides:
        return df
    if fov_n_strips is None:
        fov_n_strips = {}
    if gmm_boundaries is None:
        gmm_boundaries = {}

    def _boundary_x_at(b, y):
        if isinstance(b, (int, float)):
            return np.full_like(y, float(b), dtype=np.float64)
        (x1, y1), (x2, y2) = b
        if y2 == y1:
            return np.full_like(y, (x1 + x2) / 2.0, dtype=np.float64)
        return x1 + (y - y1) * (x2 - x1) / (y2 - y1)

    def _describe(b):
        if b is None:
            return 'None'
        if isinstance(b, (int, float)):
            return f'vertical x={b}'
        (x1, y1), (x2, y2) = b
        return f'line ({x1:.0f},{y1:.0f}) → ({x2:.0f},{y2:.0f})'

    result = df.copy()

    for fov_id, bvals_user in overrides.items():
        strip_names = fov_n_strips.get(fov_id, ('strip_1', 'strip_2', 'strip_3'))
        n_expected  = len(strip_names) - 1

        if len(bvals_user) != n_expected:
            raise ValueError(
                f'FOV {fov_id}: expected {n_expected} boundary value(s) for strips '
                f'{strip_names}, got {len(bvals_user)}'
            )
        if n_expected == 0:
            print(f'  FOV {fov_id}: single-strip FOV — no boundaries to apply')
            continue

        # Resolve None entries from gmm_boundaries
        bvals = list(bvals_user)
        for i, v in enumerate(bvals):
            if v is None:
                fallback = gmm_boundaries.get(fov_id, [])
                if i < len(fallback) and fallback[i] is not None:
                    bvals[i] = fallback[i]
                else:
                    raise ValueError(
                        f'FOV {fov_id}: boundary {i} is None and no GMM fallback '
                        f'available (pass gmm_boundaries from reassign_strips_gmm)'
                    )

        mask = result['fov'] == fov_id
        n    = mask.sum()
        if n == 0:
            print(f'  [warn] FOV {fov_id} not found — skipping')
            continue

        x = result.loc[mask, x_col].values.astype(np.float64)
        y = result.loc[mask, y_col].values.astype(np.float64)

        # slot = count of boundaries point is at-or-right-of
        slot = np.zeros(n, dtype=int)
        for b in bvals:
            bx = _boundary_x_at(b, y)
            slot += (x >= bx).astype(int)

        new_labels = np.array([strip_names[s] for s in slot])

        old_counts = result.loc[mask, strip_col].value_counts().to_dict()
        result.loc[mask, strip_col] = new_labels
        new_counts = result.loc[mask, strip_col].value_counts().to_dict()

        desc = '  |  '.join(_describe(b) for b in bvals)
        print(f'  FOV {fov_id} strips={strip_names}:  {desc}')
        print(f'    Before: {old_counts}')
        print(f'    After:  {new_counts}')

    return result


def plot_fov_overview(df, cluster_table, min_cluster_size, label_threshold=0,
                      x_col='x_rot_px', y_col='y_rot_px'):
    """
    Pre-strip cluster overview: one panel per FOV showing DBSCAN cluster quality.

    Used after DBSCAN but before strip assignment.

    Colour scheme:
      steelblue — kept cluster (n >= min_cluster_size)
      orange    — small cluster (n < min_cluster_size), auto-excluded
      lightgrey — DBSCAN noise

    All kept clusters are annotated with `id=X  n=N` by default
    (label_threshold=0), so every candidate for manual exclusion is
    identifiable. Raise `label_threshold` to hide labels on small kept
    clusters if the plot gets too cluttered.

    Parameters
    ----------
    df : pd.DataFrame  — output of dbscan_noise_filter()
    cluster_table : pd.DataFrame  — from build_cluster_table (no strip required)
    min_cluster_size : int
    label_threshold : int
        Minimum n_transcripts for a kept cluster to be annotated.
        Default 0 = label every kept cluster.
    x_col, y_col : str  — coordinate columns (default: x_rot_px / y_rot_px)
    """
    import math
    fovs   = sorted(df['fov'].unique())
    n_fovs = len(fovs)
    n_cols = min(4, n_fovs)
    n_rows = math.ceil(n_fovs / n_cols)

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(6 * n_cols, 5 * n_rows),
                             squeeze=False)
    axes = axes.ravel()

    for i, fov_id in enumerate(fovs):
        ax        = axes[i]
        fov_df    = df[df['fov'] == fov_id]
        fov_clust = cluster_table[cluster_table['fov'] == fov_id]

        # Noise
        noise_pts = fov_df[fov_df['is_noise']]
        if len(noise_pts):
            ax.scatter(noise_pts[x_col], noise_pts[y_col],
                       s=0.1, alpha=0.15, c='lightgrey', rasterized=True)

        for _, row in fov_clust.iterrows():
            mask = (~fov_df['is_noise']) & (fov_df['dbscan_label'] == row['dbscan_label'])
            pts  = fov_df[mask]
            if len(pts) == 0:
                continue
            size     = row['n_transcripts']
            is_small = size < min_cluster_size
            colour   = 'orange' if is_small else 'steelblue'
            ax.scatter(pts[x_col], pts[y_col],
                       s=0.3, alpha=0.4, c=colour, rasterized=True)
            # Annotate every kept cluster at/above label_threshold
            if (not is_small) and size >= label_threshold:
                ax.annotate(f'id={int(row.global_id)}\nn={size}',
                            (pts[x_col].mean(), pts[y_col].mean()),
                            fontsize=6, ha='center', va='center',
                            color='black',
                            bbox=dict(boxstyle='round,pad=0.15',
                                      facecolor='white', edgecolor='none',
                                      alpha=0.7))

        ax.set_title(f'FOV {fov_id}', fontsize=9)
        ax.set_aspect('equal')
        ax.set_xlabel(x_col, fontsize=7)
        ax.set_ylabel(y_col, fontsize=7)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle(
        f'Pre-strip cluster overview  |  '
        f'blue=kept  orange=small(<{min_cluster_size})  grey=noise',
        fontsize=11
    )
    plt.tight_layout()
    plt.show()

    ct = cluster_table.copy()
    ct['small?'] = ct['n_transcripts'] < min_cluster_size
    print(ct.to_string(index=False))
    n_small = ct['small?'].sum()
    print(f'\n{n_small} small clusters will be auto-excluded')


def plot_strip_tissue_overview(df, x_col='x_rot_px', y_col='y_rot_px',
                               strip_col='strip'):
    """
    Combined scatter of all clean transcripts across all FOVs, coloured by strip.

    Call this after strip assignment (reassign_strips_gmm + apply_strip_overrides)
    as a final sanity check before window assignment.

    Noise and manually_excluded transcripts are greyed out if those columns exist.

    Parameters
    ----------
    df : pd.DataFrame  — must have strip_col; optionally is_noise / manually_excluded
    x_col, y_col : str — coordinate columns (default: x_rot_px / y_rot_px)
    strip_col : str
    """
    strip_colours = {'strip_1': '#e41a1c', 'strip_2': '#377eb8', 'strip_3': '#4daf4a'}

    fig, ax = plt.subplots(figsize=(12, 10))

    # Grey background: noise / excluded
    mask_bad = pd.Series(False, index=df.index)
    if 'is_noise' in df.columns:
        mask_bad |= df['is_noise']
    if 'manually_excluded' in df.columns:
        mask_bad |= df['manually_excluded']

    bad_pts = df[mask_bad]
    if len(bad_pts):
        rng = np.random.default_rng(1)
        s   = min(30_000, len(bad_pts))
        idx = rng.choice(len(bad_pts), s, replace=False)
        ax.scatter(bad_pts.iloc[idx][x_col], bad_pts.iloc[idx][y_col],
                   s=0.1, alpha=0.1, c='lightgrey', rasterized=True)

    # Clean transcripts coloured by strip
    df_clean = df[~mask_bad]
    for strip in ['strip_1', 'strip_2', 'strip_3']:
        sub = df_clean[df_clean[strip_col] == strip]
        if len(sub) == 0:
            continue
        rng = np.random.default_rng(0)
        s   = min(50_000, len(sub))
        idx = rng.choice(len(sub), s, replace=False)
        ax.scatter(sub.iloc[idx][x_col], sub.iloc[idx][y_col],
                   s=0.5, alpha=0.4, c=strip_colours[strip],
                   label=f'{strip}  (n={len(sub):,})',
                   rasterized=True)

    ax.set_title('All FOVs — final strip assignments (clean transcripts)',
                 fontsize=12)
    ax.set_xlabel(x_col, fontsize=9)
    ax.set_ylabel(y_col, fontsize=9)
    ax.set_aspect('equal')
    ax.legend(markerscale=10, fontsize=10, loc='upper left')
    plt.tight_layout()
    plt.show()


def dbscan_noise_filter(
    df,
    x_col='x_global_px_transformed',
    y_col='y_global_px_transformed',
    group_cols=None,
    percentile=97,
    eps_floor=20,
    eps_ceiling=30,
    min_samples=5,
    min_cluster_size=150,
    diagnostics=False,
):
    """
    Flag spatial noise transcripts using adaptive DBSCAN.
    
    Derives eps per group from the 1-NN distance distribution at a given
    percentile, clipped to [eps_floor, eps_ceiling]. Points labelled as
    DBSCAN noise (cluster = -1) are flagged but not removed. Surviving
    clusters smaller than min_cluster_size transcripts are additionally
    flagged as is_small_cluster.
    
    Parameters
    ----------
    df : pd.DataFrame
        Transcript dataframe with coordinate columns.
    x_col, y_col : str
        Column names for spatial coordinates.
    group_cols : list of str or None
        Columns defining independent groups for DBSCAN (e.g. ['fov', 'strip']).
        If None, runs on the entire dataframe as one group.
    percentile : float
        Percentile of 1-NN distance distribution used to derive eps.
    eps_floor, eps_ceiling : float
        Bounds for the adaptive eps value.
    min_samples : int
        DBSCAN min_samples parameter.
    min_cluster_size : int
        Clusters with fewer than this many transcripts are flagged as
        is_small_cluster. Default 500. Set to 0 to disable.
    diagnostics : bool
        If True, prints a summary table and plots each group.
        
    Returns
    -------
    pd.DataFrame
        Copy of input with added columns:
        - 'dbscan_label': cluster label (-1 = noise)
        - 'is_noise': bool flag (DBSCAN -1 labels)
        - 'is_small_cluster': bool flag (surviving clusters below min_cluster_size)
        - 'eps_used': eps value applied to that group
        - 'eps_raw': raw adaptive eps before clipping
    """
    result = df.copy()
    result['dbscan_label'] = -1
    result['is_noise'] = True
    result['is_small_cluster'] = False
    result['eps_used'] = np.nan
    result['eps_raw'] = np.nan
    
    # Define groups
    if group_cols is None:
        group_cols = ['_dummy']
        result['_dummy'] = 'all'
    
    groups = result.groupby(group_cols)
    n_groups = len(groups)
    
    # Storage for summary
    summary_records = []
    
    # Set up diagnostic plot if needed
    if diagnostics:
        n_cols = min(3, n_groups)
        n_rows = int(np.ceil(n_groups / n_cols))
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
        if n_groups == 1:
            axes = np.array([axes])
        axes = axes.ravel()
    
    for i, (group_key, group_df) in enumerate(groups):
        idx = group_df.index
        coords = group_df[[x_col, y_col]].values
        
        # Adaptive eps
        nn = NearestNeighbors(n_neighbors=2).fit(coords)
        distances, _ = nn.kneighbors(coords)
        nn_dist = distances[:, 1]
        raw_eps = np.percentile(nn_dist, percentile)
        eps = np.clip(raw_eps, eps_floor, eps_ceiling)
        
        # DBSCAN
        db = DBSCAN(eps=eps, min_samples=min_samples).fit(coords)
        labels = db.labels_
        noise_mask = labels == -1
        
        # Small-cluster flag: surviving clusters below min_cluster_size
        small_cluster_mask = np.zeros(len(labels), dtype=bool)
        if min_cluster_size > 0:
            unique_labels, counts = np.unique(labels[labels >= 0], return_counts=True)
            small_labels = unique_labels[counts < min_cluster_size]
            small_cluster_mask = np.isin(labels, small_labels)
        
        # Store results
        result.loc[idx, 'dbscan_label'] = labels
        result.loc[idx, 'is_noise'] = noise_mask
        result.loc[idx, 'is_small_cluster'] = small_cluster_mask
        result.loc[idx, 'eps_used'] = eps
        result.loc[idx, 'eps_raw'] = raw_eps
        
        # Summary
        n_noise = noise_mask.sum()
        n_small = small_cluster_mask.sum()
        pct = 100 * n_noise / len(labels)
        label = group_key if isinstance(group_key, str) else ', '.join(str(g) for g in group_key)
        summary_records.append({
            'group': label,
            'n_total': len(labels),
            'n_noise': n_noise,
            'n_small_cluster': n_small,
            'pct_noise': round(pct, 1),
            'eps_raw': round(raw_eps, 1),
            'eps_used': round(eps, 1),
        })
        
        # Diagnostic plot
        if diagnostics:
            kept_mask = ~noise_mask & ~small_cluster_mask
            ax = axes[i]
            ax.scatter(coords[kept_mask, 0], coords[kept_mask, 1],
                       s=0.3, alpha=0.3, c='steelblue', label='kept')
            ax.scatter(coords[noise_mask, 0], coords[noise_mask, 1],
                       s=2, alpha=0.8, c='red', label='noise')
            ax.scatter(coords[small_cluster_mask, 0], coords[small_cluster_mask, 1],
                       s=2, alpha=0.8, c='orange', label=f'small (<{min_cluster_size})')
            ax.set_title(f'{label}\neps={eps:.1f} (raw={raw_eps:.1f}) | {n_noise} noise ({pct:.1f}%) | {n_small} small')
            ax.set_aspect('equal')
            ax.legend(markerscale=5, fontsize=7)
    
    # Clean up dummy column
    if '_dummy' in result.columns:
        result.drop(columns='_dummy', inplace=True)
    
    if diagnostics:
        # Hide unused axes
        for j in range(i + 1, len(axes)):
            axes[j].set_visible(False)
        plt.tight_layout()
        plt.show()
        
        print(pd.DataFrame(summary_records).to_string(index=False))
    
    return result


def build_cluster_table(df):
    """
    Build a summary table of all surviving DBSCAN clusters with globally unique IDs.

    Works both before strip assignment (groups by fov + dbscan_label) and after
    (groups by fov + strip + dbscan_label). Strip column presence is auto-detected.

    Parameters
    ----------
    df : pd.DataFrame
        Output of dbscan_noise_filter().

    Returns
    -------
    pd.DataFrame
        Columns: global_id, fov, [strip,] dbscan_label, n_transcripts
    """
    non_noise = df[~df['is_noise']]
    has_strip = ('strip' in non_noise.columns and
                 non_noise['strip'].notna().all() and
                 (non_noise['strip'] != '').all())

    records = []
    gid = 0
    if has_strip:
        for (fov, strip, label), group in non_noise.groupby(['fov', 'strip', 'dbscan_label']):
            records.append({'global_id': gid, 'fov': fov, 'strip': strip,
                            'dbscan_label': int(label), 'n_transcripts': len(group)})
            gid += 1
    else:
        for (fov, label), group in non_noise.groupby(['fov', 'dbscan_label']):
            records.append({'global_id': gid, 'fov': fov,
                            'dbscan_label': int(label), 'n_transcripts': len(group)})
            gid += 1
    return pd.DataFrame(records)


def plot_cluster_overview(df, cluster_table, min_cluster_size=150,
                          label_threshold=5000, exclude_fov_strips=None):
    """
    Show all 3 strips in global coordinates, coloured by cluster size category.

    Colour scheme:
      steelblue  — normal (large) cluster, kept
      orange     — small cluster (below min_cluster_size), auto-excluded
      red        — FOV-strip excluded via exclude_fov_strips
      lightgrey  — DBSCAN noise (already flagged)

    Annotations (id=X n=N) are only drawn for clusters with n >= label_threshold
    to keep the plot readable when tissue has hundreds of sub-clusters.

    Parameters
    ----------
    df : pd.DataFrame  — output of dbscan_noise_filter()
    cluster_table : pd.DataFrame  — output of build_cluster_table()
    min_cluster_size : int
    label_threshold : int
        Only annotate cluster centroids if n_transcripts >= this. Default 5000.
    exclude_fov_strips : list of (fov, strip) tuples or None
        FOV-strip pairs highlighted in red (will be excluded).
    """
    if exclude_fov_strips is None:
        exclude_fov_strips = []
    excl_set = set(map(tuple, exclude_fov_strips))

    strips = ['strip_1', 'strip_2', 'strip_3']
    fig, axes = plt.subplots(1, 3, figsize=(24, 8))
    cmap = cm.get_cmap('tab20', 20)

    for ax, strip_name in zip(axes, strips):
        strip_df = df[df['strip'] == strip_name]
        strip_clusters = cluster_table[cluster_table['strip'] == strip_name]

        # DBSCAN noise in grey
        noise_pts = strip_df[strip_df['is_noise']]
        if len(noise_pts):
            ax.scatter(noise_pts['x_global_px'], noise_pts['y_global_px'],
                       s=0.1, alpha=0.15, c='lightgrey', rasterized=True)

        for _, row in strip_clusters.iterrows():
            gid_val  = int(row['global_id'])
            fov      = row['fov']
            label    = row['dbscan_label']
            size     = row['n_transcripts']
            is_small = size < min_cluster_size
            is_excl  = (fov, strip_name) in excl_set

            mask = (
                (~strip_df['is_noise']) &
                (strip_df['fov'] == fov) &
                (strip_df['dbscan_label'] == label)
            )
            pts = strip_df[mask]
            if len(pts) == 0:
                continue

            if is_excl:
                colour, s_size, alpha = 'red', 0.5, 0.6
            elif is_small:
                colour, s_size, alpha = 'orange', 3, 0.9
            else:
                colour, s_size, alpha = cmap(gid_val % 20), 0.3, 0.4

            ax.scatter(pts['x_global_px'], pts['y_global_px'],
                       s=s_size, alpha=alpha, c=[colour], rasterized=True)

            # Only annotate large clusters — keeps plot readable
            if size >= label_threshold:
                cx = pts['x_global_px'].mean()
                cy = pts['y_global_px'].mean()
                ax.annotate(
                    f'id={gid_val}\nn={size}',
                    (cx, cy), fontsize=5.5, ha='center', va='center',
                    color='darkred' if is_excl else 'black',
                    fontweight='bold',
                )

        ax.set_title(strip_name, fontsize=12)
        ax.set_aspect('equal')
        ax.set_xlabel('x_global_px')
        ax.set_ylabel('y_global_px')

    legend_items = [
        plt.scatter([], [], s=10, c='steelblue',  alpha=0.7, label='kept'),
        plt.scatter([], [], s=10, c='orange',     alpha=0.9, label=f'small (<{min_cluster_size})'),
        plt.scatter([], [], s=10, c='red',        alpha=0.7, label='FOV-strip excluded'),
        plt.scatter([], [], s=10, c='lightgrey',  alpha=0.5, label='DBSCAN noise'),
    ]
    fig.legend(handles=legend_items, loc='lower center', ncol=4, fontsize=9,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(
        f'Cluster overview  |  annotating clusters with n ≥ {label_threshold}',
        fontsize=13,
    )
    plt.tight_layout()
    plt.show()

    # Summary table
    ct = cluster_table.copy()
    ct['small?'] = ct['n_transcripts'] < min_cluster_size
    ct['fov_excl?'] = ct.apply(
        lambda r: (r['fov'], r['strip']) in excl_set, axis=1
    )
    print(ct.to_string(index=False))
    n_small = (ct['n_transcripts'] < min_cluster_size).sum()
    n_excl  = ct['fov_excl?'].sum()
    print(f'\n{n_small} small clusters | {n_excl} FOV-strip excluded cluster(s)')


def apply_cleanup(df, cluster_table, min_cluster_size, also_exclude=None,
                  exclude_fov_strips=None):
    """
    Add a manually_excluded column to df.

    Works with or without a strip column in cluster_table (auto-detected).

    A transcript is manually_excluded if:
      - Its cluster has fewer than min_cluster_size transcripts, OR
      - Its cluster global_id is in also_exclude, OR
      - Its (fov, strip) pair is in exclude_fov_strips (only when strip assigned).

    Parameters
    ----------
    df : pd.DataFrame
    cluster_table : pd.DataFrame  (from build_cluster_table)
    min_cluster_size : int
    also_exclude : list of int  (global_id values) or None
    exclude_fov_strips : list of (fov, strip) tuples or None

    Returns
    -------
    pd.DataFrame with updated is_small_cluster and manually_excluded columns.
    """
    if also_exclude is None:
        also_exclude = []
    if exclude_fov_strips is None:
        exclude_fov_strips = []

    has_strip = 'strip' in cluster_table.columns
    excl_set  = set(map(tuple, exclude_fov_strips))

    result = df.copy()
    result['is_small_cluster']  = False
    result['manually_excluded'] = False

    for _, row in cluster_table.iterrows():
        gid_val = int(row['global_id'])
        fov     = row['fov']
        label   = row['dbscan_label']
        size    = row['n_transcripts']

        mask = (
            ~result['is_noise'] &
            (result['fov'] == fov) &
            (result['dbscan_label'] == label)
        )
        if has_strip:
            mask &= (result['strip'] == row['strip'])

        if size < min_cluster_size:
            result.loc[mask, 'is_small_cluster']  = True
            result.loc[mask, 'manually_excluded'] = True

        if gid_val in also_exclude:
            result.loc[mask, 'manually_excluded'] = True

    # FOV-strip exclusions (only meaningful after strip assignment)
    if excl_set:
        if not has_strip:
            print('  [warn] exclude_fov_strips ignored — strip not yet assigned')
        else:
            fov_strip_mask = result.apply(
                lambda r: (r['fov'], r['strip']) in excl_set, axis=1
            )
            result.loc[fov_strip_mask, 'manually_excluded'] = True
            print(f'  FOV-strip exclusions: {fov_strip_mask.sum():,} transcripts '
                  f'across {len(excl_set)} (fov, strip) pair(s)')

    return result


def load_lr_index(path: str) -> dict:
    """
    Load the LR panel index JSON.

    Parameters
    ----------
    path : str
        Path to lr_panel_index.json.

    Returns
    -------
    dict
        Keyed by 'LIGAND|RECEPTOR'. Each value is a dict mapping strip name
        to {'ligand_n': int, 'receptor_n': int, 'viable': bool}.
    """
    with open(path) as f:
        return json.load(f)


def query_lr_panel(
    index: dict,
    genes_present,
    strip: str = None,
    min_n: int = 50,
) -> list:
    """
    Return viable (ligand, receptor) pairs given a set of observed genes.

    Parameters
    ----------
    index : dict
        Output of load_lr_index().
    genes_present : set or list of str
        Gene names observed in the target tissue.
    strip : str or None
        If given, only consider viability in that strip ('strip_1', 'strip_2',
        'strip_3'). If None, a pair is included if viable in any strip.
    min_n : int
        Minimum transcript count to call a gene viable in a strip.
        Uses the stored n values; set to 0 to use only the pre-computed
        viable flags without reapplying a threshold.

    Returns
    -------
    list of (ligand, receptor) tuples
    """
    genes_present = set(genes_present)
    strips_to_check = [strip] if strip else ['strip_1', 'strip_2', 'strip_3']
    results = []

    for key, strip_data in index.items():
        lig, rec = key.split('|')
        if lig not in genes_present or rec not in genes_present:
            continue
        for s in strips_to_check:
            sd = strip_data.get(s, {})
            if min_n > 0:
                viable = (sd.get('ligand_n', 0) >= min_n and
                          sd.get('receptor_n', 0) >= min_n)
            else:
                viable = sd.get('viable', False)
            if viable:
                results.append((lig, rec))
                break

    return results
