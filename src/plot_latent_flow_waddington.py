#!/usr/bin/env python

import argparse
import re
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.patches import Patch
from scipy.ndimage import gaussian_filter
import numpy as np
import pandas as pd
import seaborn as sns

def truncate_colormap(cmap_name, minval=0.10, maxval=0.90, n=256):
    cmap = mpl.colormaps[cmap_name]
    colors = cmap(np.linspace(minval, maxval, n))
    return mpl.colors.LinearSegmentedColormap.from_list(
        f"{cmap_name}_trunc",
        colors,
        N=n,
    )


ROOT_COLOR = "dodgerblue"
ROOT_SIZE = 60

LANDSCAPE_TRAJ_COLOR = "0.15"
LANDSCAPE_TRAJ_LW = 0.7

SHOW_TREE_LINES_ON_LANDSCAPES = False
FACTOR_NODE_SIZE = 6
FACTOR_TIP_SIZE = 10

WADDINGTON_CMAP = truncate_colormap("Greys", minval=0.15, maxval=0.75)
WADDINGTON_ALPHA_SURF = 0.55
WADDINGTON_SHADE_ALPHA = 0.12
WADDINGTON_HILL_SMOOTH_SIGMA = 1.5
# Waddington outputs use the within-depth density of reconstructed tree edges.
# Each unique edge is counted once, and every tree-depth slice is normalized
# independently so coalesced ancestral branches do not dominate later lineages.
WADDINGTON_TRAJECTORY_DENSITY_BINS = 90
WADDINGTON_TRAJECTORY_DENSITY_SMOOTH_SIGMA = 1.0
WADDINGTON_TRAJECTORY_SAMPLES_PER_BIN = 2
WADDINGTON_TRAJECTORY_SAMPLES_MAX = 100
# In the Waddington-style view, high-density trajectory regions are drawn as
# downward deflections of the terrain lines. Because tree depth is plotted with
# root at the top, adding to y makes the line dip downward visually.
WADDINGTON_VALLEY_STRENGTH = 0.12
WADDINGTON_LINE_EVERY = 2
WADDINGTON_LINE_COLOR = "0.25"
WADDINGTON_LINE_ALPHA = 0.28
WADDINGTON_LINE_LW = 0.35

# Waddington-specific terminal-state density ledge.
# This replaces the separate below-axis histogram in the Waddington outputs.
# High terminal density is drawn as a deeper downward cut at the front edge.
WADDINGTON_DENSITY_LEDGE = True
WADDINGTON_DENSITY_LEDGE_HEIGHT = 0.15
WADDINGTON_DENSITY_LEDGE_BINS = 75
WADDINGTON_DENSITY_LEDGE_SMOOTH_SIGMA = 1.0
WADDINGTON_DENSITY_LEDGE_FACE_COLOR = "0.93"
WADDINGTON_DENSITY_LEDGE_EDGE_COLOR = "0.35"
WADDINGTON_DENSITY_LEDGE_LINE_COLOR = "0.45"
WADDINGTON_DENSITY_LEDGE_ALPHA = 0.75
WADDINGTON_DENSITY_LEDGE_LINE_ALPHA = 0.28
WADDINGTON_DENSITY_LEDGE_LW = 0.35
WADDINGTON_DENSITY_LEDGE_N_HORIZONTAL_LINES = 6
WADDINGTON_DENSITY_LEDGE_VALLEY_STRENGTH = 0.70

WADDINGTON_LABEL_HIST_BINS = 65
WADDINGTON_LABEL_HIST_HEIGHT = 0.84
WADDINGTON_LABEL_HIST_Y_OFFSET = -1.18
WADDINGTON_LABEL_HIST_ALPHA = 0.78

parser = argparse.ArgumentParser()
parser.add_argument("latent_tree_nodes_tsv")
parser.add_argument("out_prefix")
parser.add_argument("--fit-summary", required=True)
parser.add_argument("--active-variance-threshold", type=float, default=1e-5)
parser.add_argument("--min_delta_quantile", type=float, default=0.0)
parser.add_argument("--waddington-density-bins", type=int, default=WADDINGTON_TRAJECTORY_DENSITY_BINS)
parser.add_argument("--waddington-density-smooth-sigma", type=float, default=WADDINGTON_TRAJECTORY_DENSITY_SMOOTH_SIGMA)
parser.add_argument("--waddington-samples-per-bin", type=float, default=WADDINGTON_TRAJECTORY_SAMPLES_PER_BIN)
parser.add_argument("--waddington-samples-max", type=int, default=WADDINGTON_TRAJECTORY_SAMPLES_MAX)
parser.add_argument("--waddington-relief-smooth-sigma", type=float, default=WADDINGTON_HILL_SMOOTH_SIGMA)
parser.add_argument("--waddington-valley-strength", type=float, default=WADDINGTON_VALLEY_STRENGTH)
parser.add_argument("--waddington-line-every", type=int, default=WADDINGTON_LINE_EVERY)
parser.add_argument("--waddington-line-alpha", type=float, default=WADDINGTON_LINE_ALPHA)
parser.add_argument("--waddington-line-width", type=float, default=WADDINGTON_LINE_LW)
parser.add_argument("--waddington-surface-alpha", type=float, default=WADDINGTON_ALPHA_SURF)
parser.add_argument("--waddington-shade-alpha", type=float, default=WADDINGTON_SHADE_ALPHA)
parser.add_argument("--waddington-ledge-height", type=float, default=WADDINGTON_DENSITY_LEDGE_HEIGHT)
parser.add_argument("--waddington-ledge-bins", type=int, default=WADDINGTON_DENSITY_LEDGE_BINS)
parser.add_argument("--waddington-ledge-smooth-sigma", type=float, default=WADDINGTON_DENSITY_LEDGE_SMOOTH_SIGMA)
parser.add_argument("--waddington-ledge-valley-strength", type=float, default=WADDINGTON_DENSITY_LEDGE_VALLEY_STRENGTH)
parser.add_argument("--waddington-label-hist-bins", type=int, default=WADDINGTON_LABEL_HIST_BINS)
parser.add_argument("--waddington-label-hist-height", type=float, default=WADDINGTON_LABEL_HIST_HEIGHT)
parser.add_argument("--waddington-label-hist-y-offset", type=float, default=WADDINGTON_LABEL_HIST_Y_OFFSET)
parser.add_argument(
    "--barcode-discrete-label-tsv",
    help=(
        "Optional TSV with barcode and tissue columns. Adds one aligned terminal-cell "
        "histogram per tissue below each Waddington-style panel."
    ),
)
args = parser.parse_args()

# The dedicated tuned-Waddington entry point exposes these as command-line
# controls. Defaults exactly reproduce the settings used by this full script.
WADDINGTON_TRAJECTORY_DENSITY_BINS = args.waddington_density_bins
WADDINGTON_TRAJECTORY_DENSITY_SMOOTH_SIGMA = args.waddington_density_smooth_sigma
WADDINGTON_TRAJECTORY_SAMPLES_PER_BIN = args.waddington_samples_per_bin
WADDINGTON_TRAJECTORY_SAMPLES_MAX = args.waddington_samples_max
WADDINGTON_HILL_SMOOTH_SIGMA = args.waddington_relief_smooth_sigma
WADDINGTON_VALLEY_STRENGTH = args.waddington_valley_strength
WADDINGTON_LINE_EVERY = args.waddington_line_every
WADDINGTON_LINE_ALPHA = args.waddington_line_alpha
WADDINGTON_LINE_LW = args.waddington_line_width
WADDINGTON_ALPHA_SURF = args.waddington_surface_alpha
WADDINGTON_SHADE_ALPHA = args.waddington_shade_alpha
WADDINGTON_DENSITY_LEDGE_HEIGHT = args.waddington_ledge_height
WADDINGTON_DENSITY_LEDGE_BINS = args.waddington_ledge_bins
WADDINGTON_DENSITY_LEDGE_SMOOTH_SIGMA = args.waddington_ledge_smooth_sigma
WADDINGTON_DENSITY_LEDGE_VALLEY_STRENGTH = args.waddington_ledge_valley_strength
WADDINGTON_LABEL_HIST_BINS = args.waddington_label_hist_bins
WADDINGTON_LABEL_HIST_HEIGHT = args.waddington_label_hist_height
WADDINGTON_LABEL_HIST_Y_OFFSET = args.waddington_label_hist_y_offset


sns.set_theme(
    context="paper",
    style="white",
    font_scale=1.15,
    rc={
        "axes.edgecolor": "0.15",
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
    },
)

mpl.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


df = pd.read_csv(args.latent_tree_nodes_tsv, sep="\t")

factor_cols = [c for c in df.columns if c.startswith("factor_")]
if len(factor_cols) == 0:
    raise ValueError("No factor_* columns found.")

summary = pd.read_csv(args.fit_summary, sep="\t", index_col=0)
active_factor_cols = []
for parameter, value in summary["value"].items():
    match = re.fullmatch(r"sigma2_latent_LF([0-9]+)", parameter)
    if match and float(value) > args.active_variance_threshold:
        active_factor_cols.append(f"factor_{int(match.group(1))}")
missing = [factor for factor in active_factor_cols if factor not in factor_cols]
if missing:
    raise ValueError("Active factors absent from latent-flow TSV: " + ", ".join(missing))
if not active_factor_cols:
    raise ValueError("No active factors found in fit summary.")
factor_cols = active_factor_cols

df["node_id"] = df["node_id"].astype(int)
df["parent_id"] = df["parent_id"].astype(int)
df["is_tip"] = df["is_tip"].astype(int)

barcode_tissue = None
if args.barcode_discrete_label_tsv:
    barcode_tissue = pd.read_csv(args.barcode_discrete_label_tsv, sep="\t", dtype=str)
    required_columns = {"barcode", "tissue"}
    missing_columns = required_columns - set(barcode_tissue.columns)
    if missing_columns:
        raise ValueError(
            "Barcode tissue TSV is missing required column(s): "
            + ", ".join(sorted(missing_columns))
        )

    barcode_tissue = barcode_tissue[["barcode", "tissue"]].dropna().copy()
    conflicting = (
        barcode_tissue.groupby("barcode", sort=False)["tissue"].nunique() > 1
    )
    if conflicting.any():
        examples = ", ".join(conflicting[conflicting].index[:5])
        raise ValueError(f"Barcodes have conflicting tissue labels: {examples}")

    barcode_tissue = barcode_tissue.drop_duplicates("barcode")
    tissue_lookup = barcode_tissue.set_index("barcode")["tissue"]
    tip_names = df.loc[df["is_tip"] == 1, "node_name"].astype(str)
    matched_tissues = tip_names.map(tissue_lookup)

    if matched_tissues.notna().sum() == 0:
        raise ValueError("No terminal node names matched the barcode tissue TSV.")

    unmatched_count = int(matched_tissues.isna().sum())
    if unmatched_count:
        print(
            f"Warning: {unmatched_count} terminal cells have no tissue label and "
            "will be omitted from tissue histograms.",
            file=sys.stderr,
        )

    df["tissue"] = df["node_name"].astype(str).map(tissue_lookup)

root = df.loc[df["parent_id"] < 0]
if len(root) != 1:
    raise ValueError("Expected exactly one root.")

root_f = root.iloc[0][factor_cols].to_numpy(float)
F = df[factor_cols].to_numpy(float)
df["latent_distance_from_root"] = np.linalg.norm(F - root_f[None, :], axis=1)
root_row = df.loc[df["parent_id"] < 0].iloc[0]

# Tree-trajectory overlays need each child's reconstructed values alongside
# those of its parent. This is derived entirely from the existing node table.
parent_df = df[
    ["node_id", "pc1", "pc2", "latent_distance_from_root"] + factor_cols
].copy()
parent_df = parent_df.rename(
    columns={
        "node_id": "parent_id",
        "pc1": "parent_pc1",
        "pc2": "parent_pc2",
        "latent_distance_from_root": "parent_latent_distance_from_root",
        **{c: f"parent_{c}" for c in factor_cols},
    }
)
edges = df[df["parent_id"] >= 0].merge(parent_df, on="parent_id", how="left")

# Preserve the existing optional edge-change filter used for trajectory lines.
child_f = edges[factor_cols].to_numpy(float)
parent_f = edges[[f"parent_{c}" for c in factor_cols]].to_numpy(float)
edges["latent_delta_norm"] = np.linalg.norm(child_f - parent_f, axis=1)
min_delta = edges["latent_delta_norm"].quantile(args.min_delta_quantile)
edges = edges[edges["latent_delta_norm"] >= min_delta]


def format_label(label):
    return label.replace("_", " ").capitalize()


def _ensure_parent_time_col(time_col):
    col = f"parent_{time_col}"
    if col not in edges.columns:
        lookup = df.set_index("node_id")[time_col].to_dict()
        edges[col] = edges["parent_id"].map(lookup)


def _scaled_relief_from_density(density_grid):
    Z = np.nan_to_num(density_grid, nan=0.0)
    Z = gaussian_filter(Z, sigma=WADDINGTON_HILL_SMOOTH_SIGMA)

    zmax = float(np.nanpercentile(Z, 98))
    if not np.isfinite(zmax) or zmax <= 0:
        zmax = 1.0

    Z = Z / zmax
    Z = np.clip(Z, 0.0, 1.0)
    return Z


def _draw_waddington_relief_landscape(ax, xedges, yedges, density_grid, norm, cmap):
    """Draw a 2D Waddington-like landscape without changing panel geometry.

    Waddington panels encode locally smoothed root-to-tip trajectory density.
    Each terminal cell contributes one complete route through the tree, so
    branches shared by many descendant cells become darker, deeper valleys.
    """
    xcenters = 0.5 * (xedges[:-1] + xedges[1:])
    ycenters = 0.5 * (yedges[:-1] + yedges[1:])
    relief = _scaled_relief_from_density(density_grid)

    # Subtle grey background gives the density field a terrain-like texture.
    # Use the same lower-origin, full-edge extent as the density image so the
    # shaded density and curved terrain lines remain visually aligned.
    ax.imshow(
        relief,
        extent=(
            float(xedges[0]),
            float(xedges[-1]),
            float(yedges[0]),
            float(yedges[-1]),
        ),
        origin="lower",
        cmap="Greys",
        alpha=WADDINGTON_SHADE_ALPHA,
        aspect="auto",
        interpolation="bilinear",
        zorder=0,
    )

    im = ax.imshow(
        density_grid,
        extent=(
            float(xedges[0]),
            float(xedges[-1]),
            float(yedges[0]),
            float(yedges[-1]),
        ),
        origin="lower",
        cmap=cmap,
        norm=norm,
        alpha=WADDINGTON_ALPHA_SURF,
        aspect="auto",
        interpolation="bilinear",
        rasterized=True,
        zorder=1,
    )

    finite_relief = relief[np.isfinite(relief)]
    if len(finite_relief) > 0 and np.nanmax(finite_relief) > 0:
        for idx in range(0, len(ycenters), WADDINGTON_LINE_EVERY):
            y = ycenters[idx]
            z = relief[idx, :].copy()
            y_offsets = WADDINGTON_VALLEY_STRENGTH * z

            for _ in range(3):
                y_sample = y + y_offsets
                row_pos = np.interp(y_sample, ycenters, np.arange(len(ycenters)))
                z_sample = np.empty_like(z)

                for j in range(len(xcenters)):
                    r = row_pos[j]
                    r0 = int(np.floor(r))
                    r1 = min(r0 + 1, relief.shape[0] - 1)
                    r0 = max(r0, 0)
                    frac = r - r0
                    z_sample[j] = (1.0 - frac) * relief[r0, j] + frac * relief[r1, j]

                y_offsets = WADDINGTON_VALLEY_STRENGTH * z_sample

            ax.plot(
                xcenters,
                y + y_offsets,
                color=WADDINGTON_LINE_COLOR,
                lw=WADDINGTON_LINE_LW,
                alpha=WADDINGTON_LINE_ALPHA,
                zorder=2,
            )

    return im


def _waddington_ledge_height(yedges):
    if not WADDINGTON_DENSITY_LEDGE:
        return 0.0
    y_range = float(yedges[-1]) - float(yedges[0])
    return WADDINGTON_DENSITY_LEDGE_HEIGHT * y_range


def _draw_waddington_density_ledge(ax, component_col, xedges, yedges):
    if not WADDINGTON_DENSITY_LEDGE:
        return

    tips = df[df["is_tip"] == 1]
    vals = tips[component_col].to_numpy(float)
    vals = vals[np.isfinite(vals)]

    if len(vals) == 0:
        return

    xmin = float(xedges[0])
    xmax = float(xedges[-1])
    y_front = float(yedges[-1])
    ledge_h = _waddington_ledge_height(yedges)
    if ledge_h <= 0:
        return

    counts, bin_edges = np.histogram(
        vals,
        bins=WADDINGTON_DENSITY_LEDGE_BINS,
        range=(xmin, xmax),
    )
    density = counts.astype(float)
    density = gaussian_filter(density, sigma=WADDINGTON_DENSITY_LEDGE_SMOOTH_SIGMA)

    if np.nanmax(density) > 0:
        density = density / np.nanmax(density)
    else:
        density[:] = 0.0

    xcenters = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    xcurve = np.r_[xmin, xcenters, xmax]
    density_curve = np.r_[0.0, density, 0.0]

    y_bottom = y_front + ledge_h
    valley_depth = (
        WADDINGTON_DENSITY_LEDGE_VALLEY_STRENGTH
        * ledge_h
        * density_curve
    )
    y_valley = np.minimum(y_front + valley_depth, y_bottom)

    # Fill the full ledge as a quiet solid face. This makes the non-valley
    # part read as the front face of the Waddington surface, not as a histogram.
    ax.fill_between(
        [xmin, xmax],
        [y_front, y_front],
        [y_bottom, y_bottom],
        color=WADDINGTON_DENSITY_LEDGE_FACE_COLOR,
        alpha=WADDINGTON_DENSITY_LEDGE_ALPHA,
        linewidth=0,
        zorder=17,
        clip_on=True,
    )

    # Lightly open the valley area only. The horizontal lines below are then
    # restricted to this density-defined cut, so they stay in the valley rather
    # than running across the solid bottom ledge.
    ax.fill_between(
        xcurve,
        y_front,
        y_valley,
        color="white",
        alpha=0.28,
        linewidth=0,
        zorder=18,
        clip_on=True,
    )

    # Use the same approximate spacing as the main Waddington terrain lines,
    # but draw only scaled copies of the density-valley outline. This keeps the
    # spacing visually continuous while keeping the lines on the valley side.
    dy_grid = float(np.median(np.diff(yedges)))
    line_step = max(dy_grid * WADDINGTON_LINE_EVERY, ledge_h / 18.0)
    baselines = np.arange(y_front + 0.55 * line_step, y_bottom, line_step)

    min_visible_depth = 0.025 * ledge_h
    inside_valley = valley_depth > min_visible_depth

    for baseline in baselines:
        level = (baseline - y_front) / ledge_h

        # A scaled copy of the valley outline. Larger levels sit deeper in the
        # ledge; masking removes lines where the density valley is too shallow.
        y_line = y_front + level * valley_depth
        y_line = np.minimum(y_line, y_bottom)

        # Only show this line where the valley is actually deep enough to
        # contain this level. This preserves a solid, unlined bottom outside
        # the terminal-density valleys.
        visible = inside_valley & (valley_depth >= level * min(ledge_h, np.nanmax(valley_depth)))
        y_line[~visible] = np.nan

        ax.plot(
            xcurve,
            y_line,
            color=WADDINGTON_DENSITY_LEDGE_LINE_COLOR,
            lw=WADDINGTON_DENSITY_LEDGE_LW,
            alpha=WADDINGTON_DENSITY_LEDGE_LINE_ALPHA,
            zorder=22,
            clip_on=True,
        )

    # Soft front seam and bottom ledge boundary. The front seam is intentionally
    # faint so the ledge feels continuous with the tree-depth landscape above.
    ax.plot(
        [xmin, xmax],
        [y_front, y_front],
        color=WADDINGTON_DENSITY_LEDGE_EDGE_COLOR,
        lw=0.40,
        alpha=0.38,
        zorder=23,
        clip_on=True,
    )

    ax.plot(
        [xmin, xmax],
        [y_bottom, y_bottom],
        color=WADDINGTON_DENSITY_LEDGE_EDGE_COLOR,
        lw=0.55,
        alpha=0.55,
        zorder=21,
        clip_on=True,
    )

    y_valley_edge = y_valley.copy()
    y_valley_edge[~inside_valley] = np.nan
    ax.plot(
        xcurve,
        y_valley_edge,
        color=WADDINGTON_DENSITY_LEDGE_EDGE_COLOR,
        lw=0.60,
        alpha=0.72,
        zorder=24,
        clip_on=True,
    )


def _apply_waddington_y_limits(ax, yedges):
    y_min = float(yedges[0])
    y_front = float(yedges[-1])
    ledge_h = _waddington_ledge_height(yedges)

    # Extend only the Waddington-style plots downward so the density ledge is
    # part of the main axis, while keeping tick labels restricted to tree depth.
    ax.set_ylim(y_front + ledge_h, y_min)

    ticks = [t for t in ax.get_yticks() if y_min <= t <= y_front]
    ax.set_yticks(ticks)


def _draw_waddington_label_histograms(ax, component_col, xedges):
    if barcode_tissue is None:
        return

    tips = df[(df["is_tip"] == 1) & df["tissue"].notna()]
    if len(tips) == 0:
        return

    labels = sorted(tips["tissue"].unique())
    colors = sns.color_palette("tab10", n_colors=len(labels))
    bin_edges = np.linspace(
        float(xedges[0]),
        float(xedges[-1]),
        WADDINGTON_LABEL_HIST_BINS + 1,
    )

    row_height = WADDINGTON_LABEL_HIST_HEIGHT / len(labels)
    legend_handles = []
    x_min = float(xedges[0])
    x_max = float(xedges[-1])
    x_ticks = [tick for tick in ax.get_xticks() if x_min <= tick <= x_max]
    x_labels = ax.xaxis.get_major_formatter().format_ticks(x_ticks)

    for row, (label, color) in enumerate(zip(labels, colors)):
        vals = tips.loc[tips["tissue"] == label, component_col].to_numpy(float)
        vals = vals[np.isfinite(vals)]
        counts, _ = np.histogram(vals, bins=bin_edges)
        max_count = max(int(counts.max()), 1)

        hist_ax = ax.inset_axes(
            [
                0.0,
                WADDINGTON_LABEL_HIST_Y_OFFSET
                + (len(labels) - row - 1) * row_height,
                1.0,
                row_height * 0.66,
            ],
            transform=ax.transAxes,
        )
        hist_ax.stairs(
            counts,
            bin_edges,
            baseline=0,
            fill=True,
            color=color,
            alpha=WADDINGTON_LABEL_HIST_ALPHA,
            linewidth=0.35,
        )
        hist_ax.set_xlim(x_min, x_max)
        hist_ax.set_ylim(0, max_count)
        hist_ax.set_yticks([0, max_count])
        hist_ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%d"))
        hist_ax.tick_params(
            axis="y",
            labelsize=4.5,
            left=True,
            length=3.0,
            width=0.8,
            direction="out",
            pad=1,
        )
        if row == len(labels) - 1:
            hist_ax.set_xticks(x_ticks, labels=x_labels)
            hist_ax.tick_params(
                axis="x",
                labelsize=6,
                bottom=True,
                length=3.0,
                width=0.8,
                direction="out",
                pad=2,
            )
        else:
            hist_ax.tick_params(axis="x", bottom=False, labelbottom=False)
        hist_ax.set_facecolor("none")
        hist_ax.spines["top"].set_visible(False)
        hist_ax.spines["right"].set_visible(False)
        hist_ax.spines["bottom"].set_color("0.55")
        hist_ax.spines["bottom"].set_linewidth(0.8)
        hist_ax.spines["left"].set_color("0.65")
        hist_ax.spines["left"].set_linewidth(0.8)
        legend_handles.append(Patch(facecolor=color, alpha=WADDINGTON_LABEL_HIST_ALPHA, label=label))

    ax.legend(
        handles=legend_handles,
        loc="center left",
        bbox_to_anchor=(
            1.04,
            WADDINGTON_LABEL_HIST_Y_OFFSET
            + 0.5 * WADDINGTON_LABEL_HIST_HEIGHT,
        ),
        ncol=1,
        frameon=False,
        fontsize=6,
        handlelength=1.0,
        handleheight=0.7,
        labelspacing=0.45,
        borderaxespad=0,
    )


def _trajectory_density_grid(component_col, time_col):
    """Estimate the within-depth density of reconstructed tree trajectories.

    Every unique parent-child edge is sampled exactly once. After smoothing,
    each horizontal time-bin is scaled by its own maximum. Thus intensity shows
    relative trajectory density within a tree-depth slice rather than raw route
    counts, which would overemphasize early branches shared by many descendants.
    """
    _ensure_parent_time_col(time_col)

    node_lookup = df.set_index("node_id")

    x_all = df[component_col].to_numpy(float)
    y_all = df[time_col].to_numpy(float)

    xmin = float(np.nanmin(x_all))
    xmax = float(np.nanmax(x_all))
    ymin = float(np.nanmin(y_all))
    ymax = float(np.nanmax(y_all))

    xpad = 0.10 * (xmax - xmin) if xmax > xmin else 0.5
    ypad = 0.10 * (ymax - ymin) if ymax > ymin else 0.5

    xedges = np.linspace(
        xmin - xpad,
        xmax + xpad,
        WADDINGTON_TRAJECTORY_DENSITY_BINS + 1,
    )
    yedges = np.linspace(
        ymin - ypad,
        ymax + ypad,
        WADDINGTON_TRAJECTORY_DENSITY_BINS + 1,
    )

    xbin = float(xedges[1] - xedges[0])
    ybin = float(yedges[1] - yedges[0])

    xs = []
    ys = []

    for child_id, child_row in node_lookup.iterrows():
        parent_id = int(child_row["parent_id"])
        if parent_id < 0 or parent_id not in node_lookup.index:
            continue

        xp = float(node_lookup.loc[parent_id, component_col])
        xc = float(child_row[component_col])
        yp = float(node_lookup.loc[parent_id, time_col])
        yc = float(child_row[time_col])

        if not (
            np.isfinite(xp) and np.isfinite(xc)
            and np.isfinite(yp) and np.isfinite(yc)
        ):
            continue

        dx_bins = abs(xc - xp) / max(xbin, 1e-12)
        dy_bins = abs(yc - yp) / max(ybin, 1e-12)

        n_samples = int(
            np.ceil(max(dx_bins, dy_bins) * WADDINGTON_TRAJECTORY_SAMPLES_PER_BIN)
        ) + 1
        n_samples = max(2, min(WADDINGTON_TRAJECTORY_SAMPLES_MAX, n_samples))

        t = np.linspace(0.0, 1.0, n_samples)
        xs.append(xp + t * (xc - xp))
        ys.append(yp + t * (yc - yp))

    if len(xs) == 0:
        density_grid = np.full(
            (WADDINGTON_TRAJECTORY_DENSITY_BINS, WADDINGTON_TRAJECTORY_DENSITY_BINS),
            np.nan,
        ).T
        return xedges, yedges, density_grid

    xs = np.concatenate(xs)
    ys = np.concatenate(ys)

    density, _, _ = np.histogram2d(xs, ys, bins=[xedges, yedges])
    density = gaussian_filter(
        density.astype(float),
        sigma=WADDINGTON_TRAJECTORY_DENSITY_SMOOTH_SIGMA,
    )

    # Normalize separately at every tree-depth bin. Using a row maximum rather
    # than a global maximum keeps later, more dispersed lineages visible while
    # retaining their relative density pattern along the component axis.
    row_max = np.max(density, axis=0, keepdims=True)
    density = np.divide(
        density,
        row_max,
        out=np.zeros_like(density),
        where=row_max > 0,
    )

    density_grid = density.T
    density_grid[density_grid <= 1e-12] = np.nan

    return xedges, yedges, density_grid


def _draw_component_tree_trajectories(ax, component_col, time_col):
    if SHOW_TREE_LINES_ON_LANDSCAPES:
        parent_component_col = f"parent_{component_col}"

        for _, e in edges.iterrows():
            ax.plot(
                [e[parent_component_col], e[component_col]],
                [e[f"parent_{time_col}"], e[time_col]],
                color=LANDSCAPE_TRAJ_COLOR,
                lw=LANDSCAPE_TRAJ_LW,
                alpha=0.45,
                solid_capstyle="round",
                zorder=3,
            )

        internal = df[df["is_tip"] == 0]
        tips = df[df["is_tip"] == 1]

        ax.scatter(
            internal[component_col], internal[time_col],
            color="0.20", s=FACTOR_NODE_SIZE, alpha=0.18,
            linewidths=0, rasterized=True, zorder=4,
        )
        ax.scatter(
            tips[component_col], tips[time_col],
            color="0.05", s=FACTOR_TIP_SIZE, alpha=0.30,
            linewidths=0, rasterized=True, zorder=5,
        )
    ax.scatter(
        root_row[component_col], root_row[time_col],
        s=ROOT_SIZE, color=ROOT_COLOR,
        edgecolors="white", linewidths=0.6, zorder=10,
    )


def _draw_factor_tree_trajectories(ax, factor_col, time_col):
    if SHOW_TREE_LINES_ON_LANDSCAPES:
        parent_factor_col = f"parent_{factor_col}"

        for _, e in edges.iterrows():
            ax.plot(
                [e[parent_factor_col], e[factor_col]],
                [e[f"parent_{time_col}"], e[time_col]],
                color=LANDSCAPE_TRAJ_COLOR,
                lw=LANDSCAPE_TRAJ_LW,
                alpha=0.45,
                solid_capstyle="round",
                zorder=3,
            )

        internal = df[df["is_tip"] == 0]
        tips = df[df["is_tip"] == 1]

        ax.scatter(
            internal[factor_col], internal[time_col],
            color="0.20", s=FACTOR_NODE_SIZE, alpha=0.18,
            linewidths=0, rasterized=True, zorder=4,
        )
        ax.scatter(
            tips[factor_col], tips[time_col],
            color="0.05", s=FACTOR_TIP_SIZE, alpha=0.30,
            linewidths=0, rasterized=True, zorder=5,
        )
    ax.scatter(
        root_row[factor_col], root_row[time_col],
        s=ROOT_SIZE, color=ROOT_COLOR,
        edgecolors="white", linewidths=0.6, zorder=10,
    )


def plot_waddington_landscapes():
    """Plot overall latent distance followed by active factors in one row."""
    time_col = "tree_depth"
    _ensure_parent_time_col(time_col)

    panels = [(
        "latent_distance_from_root",
        "Euclidean distance from root",
        "overall",
    )]
    panels.extend((factor, format_label(factor), "factor") for factor in factor_cols)
    n_panels = len(panels)

    fig, axes = plt.subplots(
        1, n_panels,
        figsize=(3.55 * n_panels, 3.35),
        squeeze=False,
    )
    fig.subplots_adjust(
        left=0.045,
        right=0.985,
        bottom=0.30,
        top=0.86,
        wspace=0.52,
    )
    axes = axes.ravel()

    for ax, (component_col, xlabel, panel_type) in zip(axes, panels):
        xedges, yedges, density_grid = _trajectory_density_grid(component_col, time_col)
        norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
        im = _draw_waddington_relief_landscape(
            ax, xedges, yedges, density_grid, norm, WADDINGTON_CMAP
        )
        cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.030)
        cbar.set_label("Relative trajectory density", fontsize=8)
        cbar.ax.tick_params(labelsize=6)
        cbar.outline.set_visible(False)

        if panel_type == "overall":
            _draw_component_tree_trajectories(ax, component_col, time_col)
        else:
            _draw_factor_tree_trajectories(ax, component_col, time_col)
        ax.set_xlim(float(xedges[0]), float(xedges[-1]))
        _apply_waddington_y_limits(ax, yedges)
        _draw_waddington_density_ledge(ax, component_col, xedges, yedges)
        ax.set_xlabel(xlabel, fontsize=9, labelpad=5)
        ax.set_ylabel(format_label(time_col), fontsize=9, labelpad=6)
        ax.tick_params(labelsize=8)
        ax.tick_params(
            axis="y",
            which="major",
            length=3.5,
            width=0.8,
            direction="out",
            left=True,
        )
        ax.tick_params(
            axis="x",
            which="major",
            bottom=True,
            length=3.5,
            width=0.8,
            direction="out",
        )
        sns.despine(ax=ax)
        _draw_waddington_label_histograms(ax, component_col, xedges)

    fig.savefig(f"{args.out_prefix}.waddington.pdf", bbox_inches="tight")
    plt.close(fig)


plot_waddington_landscapes()
