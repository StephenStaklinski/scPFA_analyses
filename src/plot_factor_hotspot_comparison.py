#!/usr/bin/env python3

import argparse
import math
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.optimize import linear_sum_assignment


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Correlate Hotspot module membership with fitted factor gene loadings."
        )
    )
    parser.add_argument(
        "--sample",
        action="append",
        nargs=4,
        metavar=("NAME", "MODULES", "BACKGROUND", "LOADINGS"),
        help=(
            "Sample display name, Hotspot module table, Hotspot autocorrelation "
            "table, and factor-loading table. May be supplied more than once."
        ),
    )
    parser.add_argument(
        "--active-sample",
        action="append",
        nargs=5,
        metavar=("NAME", "MODULES", "BACKGROUND", "LOADINGS", "SUMMARY"),
        help=(
            "As --sample, with a fit summary used to exclude factors whose "
            "Brownian variance has collapsed."
        ),
    )
    parser.add_argument(
        "--active-variance-threshold",
        type=float,
        default=1e-5,
        help="Minimum Brownian variance for an active factor (default: 1e-5).",
    )
    parser.add_argument("--columns", type=int, default=2)
    parser.add_argument("--out-prefix", required=True)
    return parser.parse_args()


def numeric_suffix(value):
    text = str(value)
    suffix = text.rsplit("_", 1)[-1]
    return (0, int(suffix)) if suffix.isdigit() else (1, text)


def order_correlations(correlations):
    """Jointly order modules and factors around their strongest assignment."""
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


def module_factor_correlations(
    modules_file,
    background_file,
    loadings_file,
    summary_file=None,
    active_variance_threshold=1e-5,
):
    modules = pd.read_table(modules_file, index_col=0)
    background = pd.read_table(background_file, index_col=0)
    loadings = pd.read_table(loadings_file, index_col=0)

    if "Module" not in modules.columns:
        raise ValueError(f"{modules_file} does not contain a Module column")

    genes = background.index.intersection(loadings.columns, sort=False)
    if len(genes) < 2:
        raise ValueError(
            f"Fewer than two genes overlap between {background_file} "
            f"and {loadings_file}"
        )

    assigned = modules["Module"].dropna()
    assigned = assigned[assigned.astype(str) != "-1"]
    module_ids = sorted(assigned.unique(), key=numeric_suffix)
    if not module_ids:
        raise ValueError(f"No assigned Hotspot modules found in {modules_file}")

    factor_ids = sorted(loadings.index, key=numeric_suffix)
    if summary_file is not None:
        summary = pd.read_table(summary_file, index_col=0)
        active_ids = []
        for parameter, value in summary["value"].items():
            match = re.fullmatch(r"sigma2_latent_LF([0-9]+)", parameter)
            if match and float(value) > active_variance_threshold:
                active_ids.append(f"factor_{int(match.group(1))}")
        missing = [factor for factor in active_ids if factor not in loadings.index]
        if missing:
            raise ValueError(
                f"Active factors from {summary_file} are absent from "
                f"{loadings_file}: {', '.join(missing)}"
            )
        if not active_ids:
            raise ValueError(f"No active factors found in {summary_file}")
        factor_ids = active_ids
    factor_names = [f"F{numeric_suffix(factor)[1]}" for factor in factor_ids]
    factor_loadings = loadings.loc[factor_ids, genes].T
    factor_loadings.columns = factor_names
    factor_loadings = factor_loadings.apply(pd.to_numeric, errors="raise")

    correlations = pd.DataFrame(
        index=module_ids, columns=factor_names, dtype=float
    )
    module_by_gene = modules["Module"].reindex(genes)
    for module_id in module_ids:
        membership = (module_by_gene == module_id).astype(float)
        if membership.nunique() < 2:
            correlations.loc[module_id] = np.nan
        else:
            correlations.loc[module_id] = factor_loadings.corrwith(membership)

    correlations.index = [f"M{module_id}" for module_id in correlations.index]
    return order_correlations(correlations)


def safe_name(name):
    return "".join(character if character.isalnum() else "_" for character in name)


def main():
    args = parse_args()
    if args.columns < 1:
        raise ValueError("--columns must be at least 1")

    out_prefix = Path(args.out_prefix)
    results = []
    sample_args = [
        (*sample, None) for sample in (args.sample or [])
    ] + [
        tuple(sample) for sample in (args.active_sample or [])
    ]
    if not sample_args:
        raise ValueError("At least one --sample or --active-sample is required")
    for name, modules_file, background_file, loadings_file, summary_file in sample_args:
        correlations = module_factor_correlations(
            modules_file,
            background_file,
            loadings_file,
            summary_file,
            args.active_variance_threshold,
        )
        correlations.to_csv(
            f"{out_prefix}.{safe_name(name)}.correlations.tsv",
            sep="\t",
            index_label="Hotspot_module",
        )
        results.append((name, correlations))

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

    for axis, (name, correlations) in zip(axes, results):
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
        axis.set_xlabel("Factor")
        axis.set_ylabel("Hotspot module")
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
