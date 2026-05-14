#!/usr/bin/env python

import argparse
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
import numpy as np
import pandas as pd
import seaborn as sns
import umap


UMAP_N_NEIGHBORS = 15
UMAP_MIN_DIST = 0.3
UMAP_RANDOM_STATE = 1

CMAP = sns.color_palette("YlOrBr", as_cmap=True)
TREE_CMAP = mpl.colormaps["RdBu_r"]

ARROW_GAP_FRAC = 0.14
ARROW_ALPHA = 0.25
ARROW_COLOR = "black"
ARROW_LW = 0.5

ROOT_COLOR = "dodgerblue"
ROOT_SIZE = 60

VECTOR_GRID_N = 40
VECTOR_BANDWIDTH_FRAC = 0.08

FIGSIZE = (5.6, 4.7)
AX_RECT = [0.12, 0.14, 0.68, 0.78]

CIRCULAR_TREE_COLS = 5

LANDSCAPE_GRID_N = 80
LANDSCAPE_CMAP = "RdBu_r"
LANDSCAPE_ALPHA_SURF = 0.85
LANDSCAPE_TRAJ_COLOR = "0.15"
LANDSCAPE_TRAJ_LW = 0.7
LANDSCAPE_TRAJ_ALPHA = 0.55

FACTOR_LANDSCAPE_COLS = 5
FACTOR_DENSITY_BINS = 80
FACTOR_DENSITY_SMOOTH_SIGMA = 2.0
FACTOR_DENSITY_CMAP = "YlOrBr"
FACTOR_DENSITY_ALPHA = 0.92
FACTOR_NODE_SIZE = 7
FACTOR_TIP_SIZE = 12


parser = argparse.ArgumentParser()
parser.add_argument("latent_tree_nodes_tsv")
parser.add_argument("out_prefix")
parser.add_argument("--min_delta_quantile", type=float, default=0.0)
args = parser.parse_args()


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

df["node_id"] = df["node_id"].astype(int)
df["parent_id"] = df["parent_id"].astype(int)
df["is_tip"] = df["is_tip"].astype(int)

root = df.loc[df["parent_id"] < 0]
if len(root) != 1:
    raise ValueError("Expected exactly one root.")

root_f = root.iloc[0][factor_cols].to_numpy(float)
F = df[factor_cols].to_numpy(float)

df["latent_distance_from_root"] = np.linalg.norm(F - root_f[None, :], axis=1)

reducer = umap.UMAP(
    n_neighbors=UMAP_N_NEIGHBORS,
    min_dist=UMAP_MIN_DIST,
    random_state=UMAP_RANDOM_STATE,
)

tip_mask = df["is_tip"] == 1
reducer.fit(df.loc[tip_mask, factor_cols].to_numpy(float))
embedding = reducer.transform(df[factor_cols].to_numpy(float))

df["umap1"] = embedding[:, 0]
df["umap2"] = embedding[:, 1]

root_row = df.loc[df["parent_id"] < 0].iloc[0]

df.to_csv(f"{args.out_prefix}.latent_tree_nodes.umap.tsv", sep="\t", index=False)

parent_df = df[["node_id", "pc1", "pc2", "umap1", "umap2"] + factor_cols].copy()
parent_df = parent_df.rename(
    columns={
        "node_id": "parent_id",
        "pc1": "parent_pc1",
        "pc2": "parent_pc2",
        "umap1": "parent_umap1",
        "umap2": "parent_umap2",
        **{c: f"parent_{c}" for c in factor_cols},
    }
)

edges = df[df["parent_id"] >= 0].merge(parent_df, on="parent_id", how="left")

child_f = edges[factor_cols].to_numpy(float)
parent_f = edges[[f"parent_{c}" for c in factor_cols]].to_numpy(float)
edges["latent_delta_norm"] = np.linalg.norm(child_f - parent_f, axis=1)

q = edges["latent_delta_norm"].quantile(args.min_delta_quantile)
edges = edges[edges["latent_delta_norm"] >= q]


def format_label(label):
    return label.replace("_", " ").capitalize()


def format_colorbar(cbar):
    formatter = mticker.StrMethodFormatter("{x:.2f}")
    cbar.ax.yaxis.set_major_formatter(formatter)
    cbar.update_ticks()
    cbar.outline.set_visible(False)


def get_plot_limits(xcol, ycol):
    xmin = float(df[xcol].min())
    xmax = float(df[xcol].max())
    ymin = float(df[ycol].min())
    ymax = float(df[ycol].max())
    xpad = 0.08 * (xmax - xmin)
    ypad = 0.08 * (ymax - ymin)
    return xmin - xpad, xmax + xpad, ymin - ypad, ymax + ypad


def format_axis(ax, xcol, ycol, prefix):
    xmin, xmax, ymin, ymax = get_plot_limits(xcol, ycol)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(xcol.upper())
    ax.set_ylabel(ycol.upper())
    if prefix == "pca":
        ax.set_aspect("equal", adjustable="box")
    sns.despine(ax=ax)


def draw_root(ax, xcol, ycol):
    ax.scatter(
        root_row[xcol], root_row[ycol],
        s=ROOT_SIZE, color=ROOT_COLOR,
        edgecolors="white", linewidths=0.5, zorder=10,
    )


def add_root_legend(ax):
    handles = [
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor=ROOT_COLOR, markeredgecolor="white",
               markeredgewidth=0.5, markersize=6, label="Root")
    ]
    ax.legend(
        handles=handles, loc="lower right",
        bbox_to_anchor=(0.98, 0.02), bbox_transform=ax.figure.transFigure,
        frameon=False, handletextpad=0.4, borderpad=0.2,
    )


def add_flow_legend(ax):
    handles = [
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor=ROOT_COLOR, markeredgecolor="white",
               markeredgewidth=0.5, markersize=6, label="Root"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="0.82", markeredgecolor="none",
               markersize=6, label="Node"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="0.35", markeredgecolor="none",
               markersize=6, label="Tip"),
    ]
    ax.legend(
        handles=handles, loc="center left",
        bbox_to_anchor=(1.03, 0.5), bbox_transform=ax.transAxes,
        frameon=False, handletextpad=0.05, borderpad=0.05,
    )


def draw_shortened_arrow(ax, x0, y0, x1, y1):
    dx = x1 - x0
    dy = y1 - y0
    if dx == 0 and dy == 0:
        return
    start_x = x0 + ARROW_GAP_FRAC * dx
    start_y = y0 + ARROW_GAP_FRAC * dy
    end_x = x1 - ARROW_GAP_FRAC * dx
    end_y = y1 - ARROW_GAP_FRAC * dy
    ax.annotate(
        "", xy=(end_x, end_y), xytext=(start_x, start_y),
        arrowprops=dict(
            arrowstyle="-|>", color=ARROW_COLOR, lw=ARROW_LW,
            alpha=ARROW_ALPHA, shrinkA=0, shrinkB=0, mutation_scale=8,
        ),
        zorder=3,
    )


def plot_colored_projection(xcol, ycol, prefix, color_col, suffix):
    tips = df[df["is_tip"] == 1]
    internal = df[df["is_tip"] == 0]
    vmax = float(df[color_col].max())

    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_axes(AX_RECT)

    ax.scatter(
        internal[xcol], internal[ycol],
        c=internal[color_col], cmap=CMAP, vmin=0, vmax=vmax,
        s=14, alpha=0.35, linewidths=0, zorder=1,
    )
    sc = ax.scatter(
        tips[xcol], tips[ycol],
        c=tips[color_col], cmap=CMAP, vmin=0, vmax=vmax,
        s=27, alpha=0.98, linewidths=0.25, edgecolors="white",
        rasterized=True, zorder=2,
    )

    draw_root(ax, xcol, ycol)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.035)
    cbar.set_label(format_label(color_col))
    format_colorbar(cbar)
    format_axis(ax, xcol, ycol, prefix)
    add_root_legend(ax)

    fig.savefig(f"{args.out_prefix}.{prefix}.{suffix}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_flow_projection(xcol, ycol, prefix):
    tips = df[df["is_tip"] == 1]
    internal = df[df["is_tip"] == 0]

    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_axes(AX_RECT)

    ax.scatter(
        internal[xcol], internal[ycol],
        color="0.82", s=14, alpha=0.55, linewidths=0, zorder=1,
    )
    ax.scatter(
        tips[xcol], tips[ycol],
        color="0.35", s=25, alpha=0.80, linewidths=0,
        rasterized=True, zorder=2,
    )

    for _, e in edges.iterrows():
        if prefix == "pca":
            x0, y0, x1, y1 = e["parent_pc1"], e["parent_pc2"], e["pc1"], e["pc2"]
        else:
            x0, y0, x1, y1 = e["parent_umap1"], e["parent_umap2"], e["umap1"], e["umap2"]
        draw_shortened_arrow(ax, x0, y0, x1, y1)

    draw_root(ax, xcol, ycol)
    format_axis(ax, xcol, ycol, prefix)
    add_flow_legend(ax)

    fig.savefig(f"{args.out_prefix}.{prefix}.flow_arrows.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_interpolated_vector_field(xcol, ycol, prefix):
    if prefix == "pca":
        x0 = edges["parent_pc1"].to_numpy(float)
        y0 = edges["parent_pc2"].to_numpy(float)
        x1 = edges["pc1"].to_numpy(float)
        y1 = edges["pc2"].to_numpy(float)
    else:
        x0 = edges["parent_umap1"].to_numpy(float)
        y0 = edges["parent_umap2"].to_numpy(float)
        x1 = edges["umap1"].to_numpy(float)
        y1 = edges["umap2"].to_numpy(float)

    xm = 0.5 * (x0 + x1)
    ym = 0.5 * (y0 + y1)
    u = x1 - x0
    v = y1 - y0

    xmin, xmax, ymin, ymax = get_plot_limits(xcol, ycol)
    gx = np.linspace(xmin, xmax, VECTOR_GRID_N)
    gy = np.linspace(ymin, ymax, VECTOR_GRID_N)
    GX, GY = np.meshgrid(gx, gy)

    U = np.zeros_like(GX)
    V = np.zeros_like(GY)
    bandwidth = VECTOR_BANDWIDTH_FRAC * max(xmax - xmin, ymax - ymin)

    for i in range(GX.shape[0]):
        for j in range(GX.shape[1]):
            dx = xm - GX[i, j]
            dy = ym - GY[i, j]
            dist2 = dx * dx + dy * dy
            w = np.exp(-0.5 * dist2 / (bandwidth * bandwidth))
            if w.sum() > 0:
                U[i, j] = np.sum(w * u) / np.sum(w)
                V[i, j] = np.sum(w * v) / np.sum(w)

    mag = np.sqrt(U * U + V * V)

    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_axes(AX_RECT)

    stream = ax.streamplot(
        gx, gy, U, V,
        density=1.5, color=mag, cmap="viridis",
        linewidth=1.0, arrowsize=1.0, maxlength=4.0,
        integration_direction="forward", zorder=2,
    )

    cbar = fig.colorbar(stream.lines, ax=ax, fraction=0.046, pad=0.035)
    cbar.set_label("Flow magnitude")
    format_colorbar(cbar)

    draw_root(ax, xcol, ycol)
    format_axis(ax, xcol, ycol, prefix)
    add_root_legend(ax)

    fig.savefig(f"{args.out_prefix}.{prefix}.vector_field.pdf", bbox_inches="tight")
    plt.close(fig)


def _make_white_center_cmap(base_cmap_name, n=256, white_frac=0.18):
    base = mpl.colormaps[base_cmap_name]
    colors = base(np.linspace(0, 1, n))
    mid = n // 2
    half_band = int(n * white_frac / 2)
    white = np.array([1.0, 1.0, 1.0, 1.0])
    for i in range(n):
        dist = abs(i - mid) / max(half_band, 1)
        if dist < 1.0:
            w = np.cos(dist * np.pi / 2) ** 2
            colors[i] = (1 - w) * colors[i] + w * white
    return mpl.colors.LinearSegmentedColormap.from_list(
        f"{base_cmap_name}_whitecenter", colors, N=n
    )


LANDSCAPE_CMAP_OBJ = _make_white_center_cmap(LANDSCAPE_CMAP)


def _ensure_parent_time_col(time_col):
    col = f"parent_{time_col}"
    if col not in edges.columns:
        lookup = df.set_index("node_id")[time_col].to_dict()
        edges[col] = edges["parent_id"].map(lookup)


def _parent_coord_col(coord_col):
    if coord_col.startswith("pc"):
        return f"parent_{coord_col}"
    if coord_col.startswith("umap"):
        return f"parent_{coord_col}"
    raise ValueError(f"Unsupported coordinate column: {coord_col}")


def _build_landscape_grid(x_col, time_col, z_col):
    x_pts = df[x_col].to_numpy(float)
    y_pts = df[time_col].to_numpy(float)
    z_pts = df[z_col].to_numpy(float)

    xpad = 0.08 * (x_pts.max() - x_pts.min())
    ypad = 0.08 * (y_pts.max() - y_pts.min())

    gx = np.linspace(x_pts.min() - xpad, x_pts.max() + xpad, LANDSCAPE_GRID_N)
    gy = np.linspace(y_pts.min() - ypad, y_pts.max() + ypad, LANDSCAPE_GRID_N)
    GX, GY = np.meshgrid(gx, gy)

    GZ = griddata(np.column_stack([x_pts, y_pts]), z_pts, (GX, GY), method="cubic")

    nan_mask = np.isnan(GZ)
    if nan_mask.any():
        GZ_lin = griddata(np.column_stack([x_pts, y_pts]), z_pts, (GX, GY), method="linear")
        GZ[nan_mask] = GZ_lin[nan_mask]

    nan_mask2 = np.isnan(GZ)
    if nan_mask2.any():
        GZ_near = griddata(np.column_stack([x_pts, y_pts]), z_pts, (GX, GY), method="nearest")
        GZ[nan_mask2] = GZ_near[nan_mask2]

    return gx, gy, GX, GY, GZ


def _draw_tree_trajectories_landscape(ax, x_col, time_col):
    parent_x_col = _parent_coord_col(x_col)

    for _, e in edges.iterrows():
        ax.plot(
            [e[parent_x_col], e[x_col]],
            [e[f"parent_{time_col}"], e[time_col]],
            color=LANDSCAPE_TRAJ_COLOR,
            lw=LANDSCAPE_TRAJ_LW,
            alpha=LANDSCAPE_TRAJ_ALPHA,
            solid_capstyle="round",
            zorder=3,
        )

    ax.scatter(
        root_row[x_col], root_row[time_col],
        s=ROOT_SIZE, color=ROOT_COLOR,
        edgecolors="white", linewidths=0.6, zorder=10,
    )


def _style_landscape_ax(ax, x_col, time_col):
    ax.invert_yaxis()
    ax.set_xlabel(x_col.upper(), fontsize=9)
    ax.set_ylabel(format_label(time_col), fontsize=9)
    sns.despine(ax=ax)


def plot_landscape_combined():
    time_cols = ["tree_depth", "latent_distance_from_root"]
    embedding_configs = [
        ("pca", "pc1", "pc2"),
        ("umap", "umap1", "umap2"),
    ]

    for time_col in time_cols:
        _ensure_parent_time_col(time_col)

    col_configs = [
        (time_cols[0], False),
        (time_cols[0], True),
        (time_cols[1], False),
        (time_cols[1], True),
    ]

    grids = {}
    norms = {}
    for prefix, x_col, z_col in embedding_configs:
        for time_col in time_cols:
            grids[(prefix, time_col)] = _build_landscape_grid(x_col, time_col, z_col)

        # Use a shared colour scale within each embedding row. This keeps the two
        # PCA panels comparable to each other and the two UMAP panels comparable
        # to each other without forcing PC2 and UMAP2 onto the same numeric scale.
        vabs = float(np.nanpercentile(
            np.abs(np.concatenate([
                grids[(prefix, time_col)][4].ravel()
                for time_col in time_cols
            ])),
            98,
        ))
        if not np.isfinite(vabs) or vabs <= 0:
            vabs = 1.0
        norms[prefix] = mpl.colors.Normalize(vmin=-vabs, vmax=vabs)

    fig = plt.figure(figsize=(20, 10))
    gs = fig.add_gridspec(
        2, 4,
        hspace=0.30, wspace=0.30,
        left=0.05, right=0.94, top=0.95, bottom=0.07,
    )

    for row_idx, (prefix, x_col, z_col) in enumerate(embedding_configs):
        norm = norms[prefix]

        for col_idx, (time_col, show_surface) in enumerate(col_configs):
            gx, gy, GX, GY, GZ = grids[(prefix, time_col)]
            ax = fig.add_subplot(gs[row_idx, col_idx])

            if show_surface:
                im = ax.pcolormesh(
                    GX, GY, GZ,
                    cmap=LANDSCAPE_CMAP_OBJ,
                    norm=norm,
                    shading="gouraud",
                    rasterized=True,
                    zorder=1,
                )
                ax.contour(
                    GX, GY, GZ,
                    levels=10,
                    colors="0.3",
                    linewidths=0.3,
                    alpha=0.35,
                    zorder=2,
                )

            _draw_tree_trajectories_landscape(ax, x_col, time_col)
            _style_landscape_ax(ax, x_col, time_col)

            if show_surface:
                cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink=0.88)
                cbar.set_label(z_col.upper(), fontsize=8)
                cbar.ax.tick_params(labelsize=7)
                cbar.outline.set_visible(False)

    fig.savefig(f"{args.out_prefix}.landscape_combined.pdf", bbox_inches="tight")
    plt.close(fig)



def _factor_density_grid(factor_col, time_col):
    """Build a smoothed 2-D density landscape for all nodes in factor-depth space."""
    x = df[factor_col].to_numpy(float)
    y = df[time_col].to_numpy(float)

    xmin = float(np.nanmin(x))
    xmax = float(np.nanmax(x))
    ymin = float(np.nanmin(y))
    ymax = float(np.nanmax(y))

    xpad = 0.10 * (xmax - xmin) if xmax > xmin else 0.5
    ypad = 0.10 * (ymax - ymin) if ymax > ymin else 0.5

    xedges = np.linspace(xmin - xpad, xmax + xpad, FACTOR_DENSITY_BINS + 1)
    yedges = np.linspace(ymin - ypad, ymax + ypad, FACTOR_DENSITY_BINS + 1)

    # Use every node: internal nodes + tips. This is the density background.
    H, _, _ = np.histogram2d(x, y, bins=[xedges, yedges])
    H = gaussian_filter(H, sigma=FACTOR_DENSITY_SMOOTH_SIGMA)

    # Convert from raw bin counts to a relative density in [0, 1]. This gives an
    # intuitive color scale shared across panels: 0 = sparse, 1 = densest region.
    hmax = float(np.nanmax(H))
    if np.isfinite(hmax) and hmax > 0:
        H = H / hmax

    return xedges, yedges, H.T


def _draw_factor_tree_trajectories(ax, factor_col, time_col):
    parent_factor_col = f"parent_{factor_col}"

    for _, e in edges.iterrows():
        ax.plot(
            [e[parent_factor_col], e[factor_col]],
            [e[f"parent_{time_col}"], e[time_col]],
            color=LANDSCAPE_TRAJ_COLOR,
            lw=LANDSCAPE_TRAJ_LW,
            alpha=LANDSCAPE_TRAJ_ALPHA,
            solid_capstyle="round",
            zorder=3,
        )

    # Light point overlay so the density remains the visual focus.
    internal = df[df["is_tip"] == 0]
    tips = df[df["is_tip"] == 1]

    ax.scatter(
        internal[factor_col], internal[time_col],
        color="0.25", s=FACTOR_NODE_SIZE, alpha=0.20,
        linewidths=0, rasterized=True, zorder=4,
    )
    ax.scatter(
        tips[factor_col], tips[time_col],
        color="0.05", s=FACTOR_TIP_SIZE, alpha=0.35,
        linewidths=0, rasterized=True, zorder=5,
    )
    ax.scatter(
        root_row[factor_col], root_row[time_col],
        s=ROOT_SIZE, color=ROOT_COLOR,
        edgecolors="white", linewidths=0.6, zorder=10,
    )


def plot_factor_landscapes(prefix):
    """Plot one factor-vs-tree-depth landscape panel for each latent factor.

    Each panel uses factor value on the x-axis, tree depth on the y-axis, and a
    smoothed density background computed from all nodes in that same plot space.
    The prefix only controls whether the output is named as the PCA or UMAP
    factor-landscape file; the plotted factor-depth landscapes are identical.
    """
    time_col = "tree_depth"
    _ensure_parent_time_col(time_col)

    n_factors = len(factor_cols)
    n_cols = min(FACTOR_LANDSCAPE_COLS, n_factors)
    n_rows = math.ceil(n_factors / FACTOR_LANDSCAPE_COLS)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(3.25 * n_cols + 0.8, 3.2 * n_rows),
        squeeze=False,
    )

    fig.subplots_adjust(
        left=0.06,
        right=0.92,
        bottom=0.07,
        top=0.94,
        wspace=0.34,
        hspace=0.48,
    )
    axes = axes.ravel()

    density_grids = {
        factor_col: _factor_density_grid(factor_col, time_col)
        for factor_col in factor_cols
    }

    norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
    last_im = None

    for ax_idx, factor_col in enumerate(factor_cols):
        ax = axes[ax_idx]
        xedges, yedges, H = density_grids[factor_col]

        last_im = ax.pcolormesh(
            xedges, yedges, H,
            cmap=FACTOR_DENSITY_CMAP,
            norm=norm,
            shading="auto",
            alpha=FACTOR_DENSITY_ALPHA,
            rasterized=True,
            zorder=1,
        )

        # Faint contours make the density landscape easier to read without
        # hiding the tree trajectories.
        xcenters = 0.5 * (xedges[:-1] + xedges[1:])
        ycenters = 0.5 * (yedges[:-1] + yedges[1:])
        if np.nanmax(H) > 0:
            ax.contour(
                xcenters, ycenters, H,
                levels=np.linspace(0.2, 0.9, 5),
                colors="0.25",
                linewidths=0.3,
                alpha=0.28,
                zorder=2,
            )

        _draw_factor_tree_trajectories(ax, factor_col, time_col)
        ax.set_xlim(float(xedges[0]), float(xedges[-1]))
        ax.set_ylim(float(yedges[-1]), float(yedges[0]))
        ax.set_title(format_label(factor_col), fontsize=10, pad=8)
        ax.set_xlabel("Factor value", fontsize=9, labelpad=5)
        ax.set_ylabel(format_label(time_col), fontsize=9, labelpad=6)
        ax.tick_params(labelsize=8)
        sns.despine(ax=ax)

    for ax in axes[n_factors:]:
        ax.axis("off")

    if last_im is not None:
        cax = fig.add_axes([0.935, 0.18, 0.015, 0.64])

        cbar = fig.colorbar(last_im, cax=cax)
        cbar.set_label("Relative node density", fontsize=9)
        cbar.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
        cbar.ax.tick_params(labelsize=8)
        cbar.outline.set_visible(False)

    fig.subplots_adjust(
        left=0.055, right=0.93, bottom=0.09, top=0.92,
        wspace=0.34, hspace=0.50,
    )
    fig.savefig(f"{args.out_prefix}.{prefix}.factor_landscapes.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_circular_factor_trees():
    node_ids = df["node_id"].astype(int).tolist()
    root_id = int(root_row["node_id"])

    children = {node_id: [] for node_id in node_ids}
    for _, row in df.iterrows():
        node_id = int(row["node_id"])
        parent_id = int(row["parent_id"])
        if parent_id >= 0:
            children[parent_id].append(node_id)

    for node_id in children:
        children[node_id] = sorted(children[node_id])

    def dfs_tip_order(node_id, ordered_tips):
        if len(children[node_id]) == 0:
            ordered_tips.append(node_id)
            return
        for child_id in children[node_id]:
            dfs_tip_order(child_id, ordered_tips)

    ordered_tips = []
    dfs_tip_order(root_id, ordered_tips)

    if len(ordered_tips) == 0:
        raise ValueError("No tips found for circular tree plotting.")

    theta = {}
    for i, tip_id in enumerate(ordered_tips):
        theta[tip_id] = 2.0 * np.pi * (i + 0.5) / len(ordered_tips)

    def assign_internal_angles(node_id):
        if node_id in theta:
            return theta[node_id]
        child_angles = [assign_internal_angles(child_id) for child_id in children[node_id]]
        theta[node_id] = float(np.mean(child_angles))
        return theta[node_id]

    assign_internal_angles(root_id)

    depth = dict(zip(df["node_id"].astype(int), df["tree_depth"].astype(float)))
    max_depth = max(depth.values())

    if max_depth <= 0:
        raise ValueError("Maximum tree_depth must be > 0 for circular tree plotting.")

    radius = {node_id: depth[node_id] / max_depth for node_id in node_ids}
    node_lookup = df.set_index("node_id")

    n_factors = len(factor_cols)
    n_cols = min(CIRCULAR_TREE_COLS, n_factors)
    n_rows = math.ceil(n_factors / CIRCULAR_TREE_COLS)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(3.1 * n_cols, 3.45 * n_rows),
        subplot_kw={"projection": "polar"},
    )
    axes = np.atleast_1d(axes).ravel()

    for ax_idx, factor_col in enumerate(factor_cols):
        ax = axes[ax_idx]

        values = dict(zip(df["node_id"].astype(int), df[factor_col].astype(float)))
        factor_abs = max(float(np.nanmax(np.abs(df[factor_col].to_numpy(float)))), 1e-12)
        norm = mpl.colors.TwoSlopeNorm(vmin=-factor_abs, vcenter=0.0, vmax=factor_abs)

        for child_id in node_ids:
            parent_id = int(node_lookup.loc[child_id, "parent_id"])
            if parent_id < 0:
                continue

            parent_theta = theta[parent_id]
            child_theta = theta[child_id]
            parent_r = radius[parent_id]
            child_r = radius[child_id]

            edge_value = 0.5 * (values[parent_id] + values[child_id])
            color = TREE_CMAP(norm(edge_value))

            if abs(child_theta - parent_theta) > 1e-12 and parent_r > 0:
                arc_theta = np.linspace(parent_theta, child_theta, 32)
                arc_r = np.full_like(arc_theta, parent_r)
                ax.plot(arc_theta, arc_r, color=color, linewidth=1.05, alpha=0.95,
                        solid_capstyle="round")

            radial_r = np.linspace(parent_r, child_r, 32)
            radial_theta = np.full_like(radial_r, child_theta)
            ax.plot(radial_theta, radial_r, color=color, linewidth=1.05, alpha=0.95,
                    solid_capstyle="round")

        ax.set_theta_direction(-1)
        ax.set_theta_offset(np.pi / 2.0)
        ax.set_ylim(0, 1.03)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        ax.spines["polar"].set_visible(False)
        ax.set_title(format_label(factor_col), pad=8)

        sm = mpl.cm.ScalarMappable(norm=norm, cmap=TREE_CMAP)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.03)
        cbar.set_label("", fontsize=8)
        cbar.ax.tick_params(labelsize=7)
        format_colorbar(cbar)

    for ax in axes[n_factors:]:
        ax.axis("off")

    fig.savefig(f"{args.out_prefix}.factor_trees.circular.pdf", bbox_inches="tight")
    plt.close(fig)


for prefix, xcol, ycol in [
    ("pca", "pc1", "pc2"),
    ("umap", "umap1", "umap2"),
]:
    plot_colored_projection(xcol, ycol, prefix, "tree_depth", "tree_depth")
    plot_colored_projection(xcol, ycol, prefix, "latent_distance_from_root", "latent_distance_from_root")
    plot_flow_projection(xcol, ycol, prefix)
    plot_interpolated_vector_field(xcol, ycol, prefix)

plot_circular_factor_trees()
plot_landscape_combined()
plot_factor_landscapes("pca")
plot_factor_landscapes("umap")
