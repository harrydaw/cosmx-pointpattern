"""STIPPLE: segmentation-free ligand-receptor co-localisation on raw transcripts.

Public entry points, grouped by stage of the pipeline. The single most useful
function for a new dataset is :func:`run_pair_analysis`, which takes a strip's
transcript table and two gene names and returns the observed L(r), its
permutation envelope and the standardised effect size (SES).

See the README quick start for a minimal end-to-end example.
"""
from .core import (
    # data access
    get_coords, load_lr_index, query_lr_panel,
    # observation windows and edge correction
    get_window, get_convex_hull, get_concave_hull, load_custom_window,
    fraction_inside_hull, sample_in_polygon,
    # bivariate K / L / permutation envelope / SES
    bivariate_k, k_to_l, compute_envelope, run_pair_analysis,
    # preprocessing (dataset-specific: rotation, strip assignment, QC)
    add_rotated_coords, reassign_strips_gmm, dbscan_noise_filter,
    build_cluster_table, apply_cleanup, apply_strip_overrides,
    # plotting / diagnostics
    plot_diagnostics, plot_strip_assignments, plot_fov_overview,
    plot_strip_tissue_overview, plot_cluster_overview,
)

__version__ = "1.0.0"

__all__ = [
    "run_pair_analysis", "bivariate_k", "k_to_l", "compute_envelope",
    "get_window", "get_convex_hull", "get_concave_hull", "load_custom_window",
    "fraction_inside_hull", "sample_in_polygon", "get_coords",
    "load_lr_index", "query_lr_panel",
    "add_rotated_coords", "reassign_strips_gmm", "dbscan_noise_filter",
    "build_cluster_table", "apply_cleanup", "apply_strip_overrides",
    "plot_diagnostics", "plot_strip_assignments", "plot_fov_overview",
    "plot_strip_tissue_overview", "plot_cluster_overview",
]
