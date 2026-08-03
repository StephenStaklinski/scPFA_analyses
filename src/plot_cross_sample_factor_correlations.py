#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.optimize import linear_sum_assignment


def parse_args():
    parser = argparse.ArgumentParser(
        description="Correlate active factor gene loadings between two samples."
    )
    parser.add_argument("--row-sample", nargs=3,
                        metavar=("LABEL", "LOADINGS", "SUMMARY"))
    parser.add_argument("--column-sample", nargs=3,
                        metavar=("LABEL", "LOADINGS", "SUMMARY"))
    parser.add_argument(
        "--comparison",
        nargs=6,
        action="append",
        metavar=(
            "ROW_LABEL", "ROW_LOADINGS", "ROW_SUMMARY",
            "COLUMN_LABEL", "COLUMN_LOADINGS", "COLUMN_SUMMARY",
        ),
        help="Pair of samples to plot; may be supplied more than once.",
    )
    parser.add_argument("--active-variance-threshold", type=float, default=1e-5)
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()
    if args.comparison:
        if args.row_sample or args.column_sample:
            parser.error("Use either --comparison or --row-sample/--column-sample")
    elif not (args.row_sample and args.column_sample):
        parser.error("--row-sample and --column-sample are required together")
    return args


def active_loadings(loadings_path, summary_path, threshold):
    loadings = pd.read_table(loadings_path, index_col=0)
    summary = pd.read_table(summary_path, index_col=0)
    active = []
    for parameter, value in summary["value"].items():
        match = re.fullmatch(r"sigma2_latent_LF([0-9]+)", parameter)
        if match and float(value) > threshold:
            active.append((int(match.group(1)), f"factor_{int(match.group(1))}"))
    missing = [name for _number, name in active if name not in loadings.index]
    if missing:
        raise ValueError(
            f"Active factors absent from {loadings_path}: {', '.join(missing)}"
        )
    if not active:
        raise ValueError(f"No active factors found in {summary_path}")
    result = loadings.loc[[name for _number, name in active]].copy()
    result.index = [f"Factor {number}" for number, _name in active]
    return result


def jointly_order(correlations):
    scores = correlations.abs().fillna(0.0)
    rows, columns = linear_sum_assignment(-scores.to_numpy(float))
    matched = sorted(
        zip(rows, columns),
        key=lambda pair: scores.iat[pair[0], pair[1]],
        reverse=True,
    )
    row_order = [row for row, _column in matched]
    column_order = [column for _row, column in matched]
    row_order.extend(
        row for row in range(scores.shape[0]) if row not in set(row_order)
    )
    column_order.extend(
        column
        for column in range(scores.shape[1])
        if column not in set(column_order)
    )
    return correlations.iloc[row_order, column_order]


def correlate(comparison, threshold):
    (
        row_label, row_loadings_path, row_summary_path,
        column_label, column_loadings_path, column_summary_path,
    ) = comparison
    row_loadings = active_loadings(
        row_loadings_path, row_summary_path, threshold
    )
    column_loadings = active_loadings(
        column_loadings_path, column_summary_path, threshold
    )

    genes = row_loadings.columns.intersection(column_loadings.columns, sort=False)
    if len(genes) < 2:
        raise ValueError("Fewer than two genes overlap between loading matrices")

    correlations = pd.DataFrame(
        [
            [
                row_loadings.loc[row_factor, genes].corr(
                    column_loadings.loc[column_factor, genes]
                )
                for column_factor in column_loadings.index
            ]
            for row_factor in row_loadings.index
        ],
        index=row_loadings.index,
        columns=column_loadings.index,
    )
    return row_label, column_label, jointly_order(correlations)


def main():
    args = parse_args()
    comparisons = args.comparison or [args.row_sample + args.column_sample]
    results = [
        correlate(comparison, args.active_variance_threshold)
        for comparison in comparisons
    ]

    out_prefix = Path(args.out_prefix)
    for index, (row_label, column_label, correlations) in enumerate(results):
        suffix = "" if len(results) == 1 else f".{row_label.lower()}_vs_{column_label.lower()}"
        correlations.to_csv(
            f"{out_prefix}{suffix}.tsv", sep="\t", index_label=row_label
        )

    sns.set_theme(style="white", context="paper", font_scale=1.25)
    panel_widths = [
        max(3.4, 0.72 * correlations.shape[1] + 1.2)
        for _row_label, _column_label, correlations in results
    ]
    height = max(
        4.2,
        max(0.65 * correlations.shape[0] + 1.8
            for _row_label, _column_label, correlations in results),
    )
    figure, axes = plt.subplots(
        1, len(results),
        figsize=(sum(panel_widths) + 2.0, height),
        squeeze=False,
    )
    for index, (axis, (row_label, column_label, correlations)) in enumerate(
        zip(axes[0], results)
    ):
        sns.heatmap(
            correlations,
            ax=axis,
            cmap="RdBu_r",
            vmin=-1,
            vmax=1,
            center=0,
            linewidths=0.35,
            linecolor="white",
            cbar=index == len(results) - 1,
            cbar_kws={"label": "Pearson correlation", "shrink": 0.8},
        )
        axis.set_xlabel(column_label)
        axis.set_ylabel(row_label)
        axis.tick_params(axis="x", labelrotation=0)
        axis.tick_params(axis="y", labelrotation=0)
    figure.tight_layout()
    figure.savefig(f"{out_prefix}.pdf", bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
