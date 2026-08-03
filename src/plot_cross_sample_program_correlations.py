#!/usr/bin/env python3
"""Compare Hotspot modules and factor loadings between two samples."""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hotspot-modules-1", required=True)
    parser.add_argument("--hotspot-modules-2", required=True)
    parser.add_argument("--hotspot-background-1", required=True)
    parser.add_argument("--hotspot-background-2", required=True)
    parser.add_argument("--factor-loadings-1", required=True)
    parser.add_argument("--factor-loadings-2", required=True)
    parser.add_argument("--sample-1", required=True)
    parser.add_argument("--sample-2", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def read_modules(path):
    modules = pd.read_table(path, index_col=0)
    if "Module" not in modules.columns:
        raise ValueError(f"{path} does not contain a Module column")
    return modules["Module"]


def hotspot_membership_matrix(
    modules,
    shared_genes,
    sample,
):
    columns = {}
    module_ids = sorted(m for m in modules.dropna().unique() if m != -1)
    for module_id in module_ids:
        label = f"M{int(module_id)}"
        member_genes = set(modules.index[modules == module_id])
        columns[label] = [int(gene in member_genes) for gene in shared_genes]

    return pd.DataFrame(columns, index=shared_genes)


def factor_loading_matrix(loadings, shared_genes, sample):
    columns = {}
    for factor in loadings.index:
        factor_number = str(factor).replace("factor_", "")
        columns[f"F{factor_number}"] = loadings.loc[factor, shared_genes]

    return pd.DataFrame(columns, index=shared_genes)


def cross_correlations(matrix_1, matrix_2):
    labels_1 = list(matrix_1.columns)
    labels_2 = list(matrix_2.columns)
    unique_1 = matrix_1.copy()
    unique_2 = matrix_2.copy()
    unique_1.columns = [f"sample1::{label}" for label in labels_1]
    unique_2.columns = [f"sample2::{label}" for label in labels_2]

    combined = pd.concat([unique_1, unique_2], axis=1).corr()
    correlation = combined.loc[unique_1.columns, unique_2.columns]
    correlation.index = labels_1
    correlation.columns = labels_2

    row_order = correlation.max(axis=1).sort_values(ascending=False).index
    column_order = correlation.max(axis=0).sort_values(ascending=False).index
    return correlation.loc[row_order, column_order]


def draw_cross_sample_heatmap(ax, correlation, title):
    sns.heatmap(
        correlation,
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        linewidths=0.25,
        linecolor="white",
        cbar_kws={"label": "Pearson correlation", "shrink": 0.72},
        ax=ax,
    )
    ax.set_title(title, fontsize=14)
    ax.set_xlabel("TLS2", fontsize=12)
    ax.set_ylabel("TLS1", fontsize=12)
    ax.tick_params(axis="x", labelrotation=90, labelsize=9)
    ax.tick_params(axis="y", labelrotation=0, labelsize=9)


def main():
    args = parse_args()
    sns.set_theme(style="white", context="paper", font_scale=1.15)

    hotspot_1 = read_modules(args.hotspot_modules_1)
    hotspot_2 = read_modules(args.hotspot_modules_2)
    background_1 = pd.read_table(args.hotspot_background_1, index_col=0)
    background_2 = pd.read_table(args.hotspot_background_2, index_col=0)
    shared_hotspot_genes = background_1.index.intersection(background_2.index)
    if len(shared_hotspot_genes) == 0:
        raise ValueError("The Hotspot runs have no shared background genes")

    hotspot_matrix_1 = hotspot_membership_matrix(
        hotspot_1,
        shared_hotspot_genes,
        args.sample_1,
    )
    hotspot_matrix_2 = hotspot_membership_matrix(
        hotspot_2,
        shared_hotspot_genes,
        args.sample_2,
    )

    loadings_1 = pd.read_table(args.factor_loadings_1, index_col=0)
    loadings_2 = pd.read_table(args.factor_loadings_2, index_col=0)
    shared_factor_genes = loadings_1.columns.intersection(loadings_2.columns)
    if len(shared_factor_genes) == 0:
        raise ValueError("The factor runs have no shared genes")

    factor_matrix_1 = factor_loading_matrix(
        loadings_1,
        shared_factor_genes,
        args.sample_1,
    )
    factor_matrix_2 = factor_loading_matrix(
        loadings_2,
        shared_factor_genes,
        args.sample_2,
    )

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(12, 11),
        constrained_layout=True,
    )
    draw_cross_sample_heatmap(
        axes[0],
        cross_correlations(hotspot_matrix_1, hotspot_matrix_2),
        "Hotspot",
    )
    draw_cross_sample_heatmap(
        axes[1],
        cross_correlations(factor_matrix_1, factor_matrix_2),
        "Factor analysis",
    )
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
