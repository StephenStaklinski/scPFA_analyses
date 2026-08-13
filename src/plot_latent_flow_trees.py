#!/usr/bin/env python

import argparse
import math
import re
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

TREE_CMAP = mpl.colormaps["RdBu_r"]
CIRCULAR_TREE_COLS = 5

parser = argparse.ArgumentParser()
parser.add_argument("latent_tree_nodes_tsv")
parser.add_argument("out_prefix")
parser.add_argument("--fit-summary")
parser.add_argument("--active-variance-threshold", type=float, default=1e-5)
parser.add_argument("--factor", type=int, action="append",
                    help="Factor/module number to plot (repeatable); bypasses --fit-summary")
parser.add_argument("--all-factors", action="store_true",
                    help="Plot every factor/module column; bypasses --fit-summary")
parser.add_argument("--title-prefix", default="Factor")
parser.add_argument("--output",
                    help="Output path (default: OUT_PREFIX.factor_trees.pdf)")
parser.add_argument("--dpi", type=int, default=300)
args = parser.parse_args()

mpl.rcParams.update({"savefig.dpi": args.dpi, "pdf.fonttype": 42, "ps.fonttype": 42})

df = pd.read_csv(args.latent_tree_nodes_tsv, sep="\t")
factor_cols = [c for c in df.columns if c.startswith("factor_")]
if not factor_cols:
    raise ValueError("No factor_* columns found.")

active_factor_cols = []
if args.all_factors:
    active_factor_cols = factor_cols
elif args.factor:
    active_factor_cols = [f"factor_{number}" for number in args.factor]
else:
    if not args.fit_summary:
        parser.error("--fit-summary is required unless --factor or --all-factors is supplied")
    summary = pd.read_csv(args.fit_summary, sep="\t", index_col=0)
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
        factor_number = factor_col.removeprefix("factor_")
        ax.set_title(f"{args.title_prefix} {factor_number}", pad=8)

        sm = mpl.cm.ScalarMappable(norm=norm, cmap=TREE_CMAP)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.03)
        cbar.set_label("", fontsize=8)
        cbar.ax.tick_params(labelsize=7)
        format_colorbar(cbar)

    for ax in axes[n_factors:]:
        ax.axis("off")

    output = args.output or f"{args.out_prefix}.factor_trees.pdf"
    fig.savefig(output, bbox_inches="tight", dpi=args.dpi)
    plt.close(fig)


plot_circular_factor_trees()
