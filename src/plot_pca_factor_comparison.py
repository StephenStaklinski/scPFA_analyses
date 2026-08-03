#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.optimize import linear_sum_assignment


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compare the top ten PCA gene loadings with fitted factor gene "
            "loadings for each sample."
        )
    )
    parser.add_argument(
        "--sample",
        action="append",
        nargs=3,
        metavar=("NAME", "PCA", "FACTORS"),
        required=True,
        help=(
            "Sample display name, PCA eigenvector table, and fitted factor-loading "
            "table. May be supplied more than once."
        ),
    )
    parser.add_argument("--columns", type=int, default=2)
    parser.add_argument("--out-prefix", required=True)
    return parser.parse_args()


def numeric_suffix(value):
    text = str(value)
    suffix = text.rsplit("_", 1)[-1]
    return (0, int(suffix)) if suffix.isdigit() else (1, text)


def order_correlations(correlations):
    """Jointly order factors and PCs around their strongest assignment."""
    scores = correlations.abs().fillna(0.0)
    matched_rows, matched_cols = linear_sum_assignment(
        -scores.to_numpy(float)
    )
    matched_pairs = sorted(
        zip(matched_rows, matched_cols),
        key=lambda pair: scores.iat[pair[0], pair[1]],
        reverse=True,
    )

    row_order = [row for row, _column in matched_pairs]
    column_order = [column for _row, column in matched_pairs]
    matched_row_set = set(row_order)
    matched_column_set = set(column_order)

    unmatched_rows = [
        row for row in range(scores.shape[0]) if row not in matched_row_set
    ]
    unmatched_rows.sort(
        key=lambda row: scores.iloc[row].max(),
        reverse=True,
    )
    unmatched_columns = [
        column
        for column in range(scores.shape[1])
        if column not in matched_column_set
    ]
    unmatched_columns.sort(
        key=lambda column: scores.iloc[:, column].max(),
        reverse=True,
    )

    row_order.extend(unmatched_rows)
    column_order.extend(unmatched_columns)
    return correlations.iloc[row_order, column_order]


def pc_factor_correlations(pca_file, factor_file):
    pca = pd.read_table(pca_file, index_col=0)
    factors = pd.read_table(factor_file, index_col=0)

    pc_names = [f"PC{i}" for i in range(1, 11)]
    missing_pcs = [pc for pc in pc_names if pc not in pca.columns]
    if missing_pcs:
        raise ValueError(f"{pca_file} is missing: {', '.join(missing_pcs)}")

    genes = pca.index.intersection(factors.columns, sort=False)
    if len(genes) < 2:
        raise ValueError(
            f"Fewer than two genes overlap between {pca_file} and {factor_file}"
        )

    factor_ids = sorted(factors.index, key=numeric_suffix)
    if not factor_ids:
        raise ValueError(f"No fitted factors found in {factor_file}")

    pc_loadings = pca.loc[genes, pc_names].apply(pd.to_numeric, errors="raise")
    factor_loadings = factors.loc[factor_ids, genes].T
    factor_loadings.columns = [f"F{i + 1}" for i in range(len(factor_ids))]
    factor_loadings = factor_loadings.apply(pd.to_numeric, errors="raise")

    correlations = pd.DataFrame(
        {
            factor: pc_loadings.corrwith(factor_loadings[factor])
            for factor in factor_loadings.columns
        }
    ).T
    correlations.index.name = "Factor"
    return order_correlations(correlations), len(genes)


def safe_name(name):
    return "".join(character if character.isalnum() else "_" for character in name)


def main():
    args = parse_args()
    if args.columns < 1:
        raise ValueError("--columns must be at least 1")

    out_prefix = Path(args.out_prefix)
    results = []
    for name, pca_file, factor_file in args.sample:
        correlations, shared_gene_count = pc_factor_correlations(
            pca_file, factor_file
        )
        correlations.to_csv(
            f"{out_prefix}.{safe_name(name)}.correlations.tsv",
            sep="\t",
        )
        results.append((name, correlations, shared_gene_count))

    ncols = min(args.columns, len(results))
    nrows = math.ceil(len(results) / ncols)
    sns.set_theme(style="white", context="paper", font_scale=1.35)
    figure, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(6.2 * ncols, 4.8 * nrows),
        squeeze=False,
    )
    axes = axes.ravel()

    for axis, (name, correlations, _shared_gene_count) in zip(axes, results):
        sns.heatmap(
            correlations,
            ax=axis,
            cmap="RdBu_r",
            vmin=-1,
            vmax=1,
            center=0,
            cbar=False,
            linewidths=0.25,
            linecolor="white",
            xticklabels=True,
            yticklabels=True,
        )
        axis.set_title(name, pad=8)
        axis.set_xlabel("Principal component")
        axis.set_ylabel("Fitted factor")
        axis.tick_params(axis="x", labelrotation=0)
        axis.tick_params(axis="y", labelrotation=0)

    for axis in axes[len(results) :]:
        axis.set_visible(False)

    figure.subplots_adjust(
        left=0.08,
        right=0.86,
        bottom=0.10,
        top=0.91,
        wspace=0.32,
        hspace=0.36,
    )
    colorbar_axis = figure.add_axes([0.89, 0.28, 0.014, 0.44])
    colorbar = figure.colorbar(
        ScalarMappable(norm=Normalize(vmin=-1, vmax=1), cmap="RdBu_r"),
        cax=colorbar_axis,
    )
    colorbar.set_label("Pearson correlation")
    figure.savefig(f"{out_prefix}.pdf", bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
