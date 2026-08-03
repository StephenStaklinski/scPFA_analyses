#!/usr/bin/env python3
"""Plot correlations between Hotspot module memberships across samples."""

import argparse
import re

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import pandas as pd
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample",
        action="append",
        nargs=3,
        metavar=("NAME", "MODULES_TSV", "BACKGROUND_TSV"),
        required=True,
        help="Sample name, Hotspot module assignments, and autocorrelation table",
    )
    parser.add_argument(
        "--pair-only",
        action="store_true",
        help="For exactly two samples, draw only sample 1 versus sample 2",
    )
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_sample(name, modules_path, background_path):
    modules_table = pd.read_table(modules_path, index_col=0)
    if "Module" not in modules_table.columns:
        raise ValueError(f"{modules_path} does not contain a Module column")
    modules = modules_table["Module"]
    background = pd.read_table(background_path, index_col=0).index
    module_ids = sorted(
        module_id
        for module_id in modules.dropna().unique()
        if module_id != -1
    )
    if not module_ids:
        raise ValueError(f"{name} has no assigned Hotspot modules")
    return {
        "name": name,
        "modules": modules,
        "background": background,
        "module_ids": module_ids,
    }


def membership_matrix(sample, genes):
    modules = sample["modules"]
    columns = {}
    for module_id in sample["module_ids"]:
        members = set(modules.index[modules == module_id])
        columns[f"M{int(module_id)}"] = [gene in members for gene in genes]
    return pd.DataFrame(columns, index=genes, dtype=float)


def pair_correlation(row_sample, column_sample):
    shared_genes = row_sample["background"].intersection(
        column_sample["background"], sort=False
    )
    if len(shared_genes) == 0:
        raise ValueError(
            f"{row_sample['name']} and {column_sample['name']} "
            "have no shared background genes"
        )
    row_membership = membership_matrix(row_sample, shared_genes)
    column_membership = membership_matrix(column_sample, shared_genes)
    row_membership.columns = [
        f"row::{column}" for column in row_membership.columns
    ]
    column_membership.columns = [
        f"column::{column}" for column in column_membership.columns
    ]
    combined = pd.concat([row_membership, column_membership], axis=1).corr()
    result = combined.loc[row_membership.columns, column_membership.columns]
    result.index = [label.replace("row::", "") for label in result.index]
    result.columns = [
        label.replace("column::", "") for label in result.columns
    ]
    row_order = result.max(axis=1).sort_values(ascending=False).index
    column_order = result.max(axis=0).sort_values(ascending=False).index
    result = result.loc[row_order, column_order]
    return result


def draw_heatmap(ax, correlation, show_x, show_y):
    sns.heatmap(
        correlation,
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        square=True,
        linewidths=0.2,
        linecolor="white",
        cbar=False,
        ax=ax,
    )
    ax.tick_params(axis="x", labelrotation=0, labelsize=10)
    ax.tick_params(axis="y", labelrotation=0, labelsize=10)
    if not show_x:
        ax.set_xticklabels([])
        ax.set_xlabel("")
    if not show_y:
        ax.set_yticklabels([])
        ax.set_ylabel("")


def draw_pair(samples, output):
    correlation = pair_correlation(samples[0], samples[1])
    width = max(7.5, 0.55 * correlation.shape[1] + 3)
    height = max(6.5, 0.45 * correlation.shape[0] + 2)
    fig, ax = plt.subplots(figsize=(width, height), constrained_layout=True)
    draw_heatmap(ax, correlation, True, True)
    ax.set_xlabel(samples[1]["name"])
    ax.set_ylabel(samples[0]["name"])
    fig.colorbar(
        ScalarMappable(norm=Normalize(-1, 1), cmap="RdBu_r"),
        ax=ax,
        label="Pearson correlation",
        shrink=0.55,
        fraction=0.035,
        pad=0.025,
    )
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def draw_grid(samples, output):
    n = len(samples)
    fig, axes = plt.subplots(
        n,
        n,
        figsize=(5.2 * n, 4.8 * n),
        constrained_layout=True,
        squeeze=False,
    )
    for row, row_sample in enumerate(samples):
        for column, column_sample in enumerate(samples):
            ax = axes[row, column]
            if row == column:
                ax.set_axis_off()
                continue
            correlation = pair_correlation(row_sample, column_sample)
            draw_heatmap(
                ax,
                correlation,
                show_x=True,
                show_y=True,
            )
            ax.tick_params(axis="x", labelrotation=0, labelsize=12)
            ax.tick_params(axis="y", labelrotation=0, labelsize=12)
    display_names = []
    for sample in samples:
        match = re.fullmatch(r"quinn_(\d+)", sample["name"])
        display_names.append(f"CP{match.group(1)}" if match else sample["name"])

    for column, display_name in enumerate(display_names):
        axes[0, column].set_title(display_name, fontsize=18)
    for row, display_name in enumerate(display_names):
        axes[row, 0].text(
            -0.18,
            0.5,
            display_name,
            transform=axes[row, 0].transAxes,
            rotation=90,
            ha="center",
            va="center",
            fontsize=17,
        )

    colorbar = fig.colorbar(
        ScalarMappable(norm=Normalize(-1, 1), cmap="RdBu_r"),
        ax=axes,
        shrink=0.28,
        fraction=0.012,
        pad=0.01,
        aspect=30,
    )
    colorbar.set_label("Pearson correlation", fontsize=14)
    colorbar.ax.tick_params(labelsize=12)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    sns.set_theme(style="white", context="paper", font_scale=1.3)
    samples = [read_sample(*sample_args) for sample_args in args.sample]
    if len(samples) < 2:
        raise ValueError("At least two --sample inputs are required")
    if args.pair_only:
        if len(samples) != 2:
            raise ValueError("--pair-only requires exactly two samples")
        draw_pair(samples, args.output)
    else:
        draw_grid(samples, args.output)


if __name__ == "__main__":
    main()
