#!/usr/bin/env python

import argparse

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
import umap


UMAP_N_NEIGHBORS = 15
UMAP_MIN_DIST = 0.3
UMAP_RANDOM_STATE = 1

CMAP = sns.color_palette("YlOrBr", as_cmap=True)

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

df["is_tip"] = df["is_tip"].astype(int)

root = df.loc[df["parent_id"] < 0]
if len(root) != 1:
    raise ValueError("Expected exactly one root.")

root_f = root.iloc[0][factor_cols].to_numpy(float)
F = df[factor_cols].to_numpy(float)

df["latent_distance_from_root"] = np.linalg.norm(
    F - root_f[None, :],
    axis=1,
)

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

df.to_csv(
    f"{args.out_prefix}.latent_tree_nodes.umap.tsv",
    sep="\t",
    index=False,
)

parent_df = df[
    ["node_id", "pc1", "pc2", "umap1", "umap2"] + factor_cols
].copy()

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

edges = df[df["parent_id"] >= 0].merge(
    parent_df,
    on="parent_id",
    how="left",
)

child_f = edges[factor_cols].to_numpy(float)
parent_f = edges[[f"parent_{c}" for c in factor_cols]].to_numpy(float)

edges["latent_delta_norm"] = np.linalg.norm(
    child_f - parent_f,
    axis=1,
)

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
        root_row[xcol],
        root_row[ycol],
        s=ROOT_SIZE,
        color=ROOT_COLOR,
        edgecolors="white",
        linewidths=0.5,
        zorder=10,
    )


def add_root_legend(ax):
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=ROOT_COLOR,
            markeredgecolor="white",
            markeredgewidth=0.5,
            markersize=6,
            label="Root",
        )
    ]

    ax.legend(
        handles=handles,
        loc="lower right",
        bbox_to_anchor=(0.98, 0.02),
        bbox_transform=ax.figure.transFigure,
        frameon=False,
        handletextpad=0.4,
        borderpad=0.2,
    )


def add_flow_legend(ax):
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=ROOT_COLOR,
            markeredgecolor="white",
            markeredgewidth=0.5,
            markersize=6,
            label="Root",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="0.82",
            markeredgecolor="none",
            markersize=6,
            label="Node",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="0.35",
            markeredgecolor="none",
            markersize=6,
            label="Tip",
        ),
    ]

    ax.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.03, 0.5),
        bbox_transform=ax.transAxes,
        frameon=False,
        handletextpad=0.05,
        borderpad=0.05,
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
        "",
        xy=(end_x, end_y),
        xytext=(start_x, start_y),
        arrowprops=dict(
            arrowstyle="-|>",
            color=ARROW_COLOR,
            lw=ARROW_LW,
            alpha=ARROW_ALPHA,
            shrinkA=0,
            shrinkB=0,
            mutation_scale=8,
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
        internal[xcol],
        internal[ycol],
        c=internal[color_col],
        cmap=CMAP,
        vmin=0,
        vmax=vmax,
        s=14,
        alpha=0.35,
        linewidths=0,
        zorder=1,
    )

    sc = ax.scatter(
        tips[xcol],
        tips[ycol],
        c=tips[color_col],
        cmap=CMAP,
        vmin=0,
        vmax=vmax,
        s=27,
        alpha=0.98,
        linewidths=0.25,
        edgecolors="white",
        rasterized=True,
        zorder=2,
    )

    draw_root(ax, xcol, ycol)

    cbar = fig.colorbar(
        sc,
        ax=ax,
        fraction=0.046,
        pad=0.035,
    )

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
        internal[xcol],
        internal[ycol],
        color="0.82",
        s=14,
        alpha=0.55,
        linewidths=0,
        zorder=1,
    )

    ax.scatter(
        tips[xcol],
        tips[ycol],
        color="0.35",
        s=25,
        alpha=0.80,
        linewidths=0,
        rasterized=True,
        zorder=2,
    )

    for _, e in edges.iterrows():
        if prefix == "pca":
            x0 = e["parent_pc1"]
            y0 = e["parent_pc2"]
            x1 = e["pc1"]
            y1 = e["pc2"]
        else:
            x0 = e["parent_umap1"]
            y0 = e["parent_umap2"]
            x1 = e["umap1"]
            y1 = e["umap2"]

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

    bandwidth = VECTOR_BANDWIDTH_FRAC * max(
        xmax - xmin,
        ymax - ymin,
    )

    for i in range(GX.shape[0]):
        for j in range(GX.shape[1]):
            dx = xm - GX[i, j]
            dy = ym - GY[i, j]

            dist2 = dx * dx + dy * dy

            w = np.exp(
                -0.5 * dist2 / (bandwidth * bandwidth)
            )

            if w.sum() > 0:
                U[i, j] = np.sum(w * u) / np.sum(w)
                V[i, j] = np.sum(w * v) / np.sum(w)

    mag = np.sqrt(U * U + V * V)

    fig = plt.figure(figsize=FIGSIZE)
    ax = fig.add_axes(AX_RECT)

    stream = ax.streamplot(
        gx,
        gy,
        U,
        V,
        density=1.5,
        color=mag,
        cmap="viridis",
        linewidth=1.0,
        arrowsize=1.0,
        maxlength=4.0,
        integration_direction="forward",
        zorder=2,
    )

    cbar = fig.colorbar(
        stream.lines,
        ax=ax,
        fraction=0.046,
        pad=0.035,
    )

    cbar.set_label("Flow magnitude")
    format_colorbar(cbar)

    draw_root(ax, xcol, ycol)

    format_axis(ax, xcol, ycol, prefix)
    add_root_legend(ax)

    fig.savefig(f"{args.out_prefix}.{prefix}.vector_field.pdf", bbox_inches="tight")
    plt.close(fig)


for prefix, xcol, ycol in [
    ("pca", "pc1", "pc2"),
    ("umap", "umap1", "umap2"),
]:
    plot_colored_projection(
        xcol,
        ycol,
        prefix,
        "tree_depth",
        "tree_depth",
    )

    plot_colored_projection(
        xcol,
        ycol,
        prefix,
        "latent_distance_from_root",
        "latent_distance_from_root",
    )

    plot_flow_projection(
        xcol,
        ycol,
        prefix,
    )

    plot_interpolated_vector_field(
        xcol,
        ycol,
        prefix,
    )