#!/usr/bin/env python

import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns

CMAP = sns.color_palette("YlOrBr", as_cmap=True)
ROOT_COLOR = "dodgerblue"
ROOT_SIZE = 60
FIGSIZE = (5.6, 4.7)
AX_RECT = [0.12, 0.14, 0.68, 0.78]

parser = argparse.ArgumentParser()
parser.add_argument("latent_tree_nodes_tsv")
parser.add_argument("out_prefix")
args = parser.parse_args()

sns.set_theme(context="paper", style="white", font_scale=1.15)
mpl.rcParams.update({"savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42})

df = pd.read_csv(args.latent_tree_nodes_tsv, sep="\t")
factor_cols = [c for c in df.columns if c.startswith("factor_")]
if not factor_cols:
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
root_row = df.loc[df["parent_id"] < 0].iloc[0]


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
    ax.set_aspect("auto")
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


plot_colored_projection("pc1", "pc2", "pca", "tree_depth", "tree_depth")
plot_colored_projection(
    "pc1", "pc2", "pca", "latent_distance_from_root", "latent_distance_from_root"
)

