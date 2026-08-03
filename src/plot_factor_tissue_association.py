#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


TISSUE_ORDER = ["M1", "M2", "RW", "RE", "LL", "Liv"]
TISSUE_COLORS = {
    "M1": "#4c78a8",
    "M2": "#f58518",
    "RW": "#54a24b",
    "RE": "#e45756",
    "LL": "#b279a2",
    "Liv": "#ffbf79",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize active factor-score associations with tissue."
    )
    parser.add_argument(
        "--sample",
        nargs=3,
        action="append",
        required=True,
        metavar=("LABEL", "FACTOR_SCORES", "FIT_SUMMARY"),
    )
    parser.add_argument("--active-variance-threshold", type=float, default=1e-5)
    parser.add_argument("--output", required=True)
    parser.add_argument("--table-output", required=True)
    return parser.parse_args()


def active_factor_names(summary_path, threshold):
    summary = pd.read_table(summary_path, index_col=0)
    active = []
    for parameter, value in summary["value"].items():
        match = re.fullmatch(r"sigma2_latent_LF([0-9]+)", parameter)
        if match and float(value) > threshold:
            number = int(match.group(1))
            active.append((number, f"factor_{number}"))
    if not active:
        raise ValueError(f"No active factors found in {summary_path}")
    return active


def summarize_sample(label, scores_path, summary_path, threshold):
    scores = pd.read_table(scores_path, index_col=0)
    active = active_factor_names(summary_path, threshold)
    missing = [name for _number, name in active if name not in scores.columns]
    if missing:
        raise ValueError(
            f"Active factors absent from {scores_path}: {', '.join(missing)}"
        )

    scores = scores[[name for _number, name in active]].copy()
    scores["tissue"] = scores.index.astype(str).str.split(".", n=1).str[0]
    tissues = [tissue for tissue in TISSUE_ORDER if tissue in set(scores["tissue"])]
    tissues += sorted(set(scores["tissue"]) - set(tissues))

    standardized_means = []
    standardized_stddevs = []
    eta_squared = []
    records = []
    for number, factor in active:
        values = scores[factor].astype(float)
        standard_deviation = values.std(ddof=0)
        means = scores.groupby("tissue", sort=False)[factor].mean()
        if standard_deviation > 0:
            z_means = (means - values.mean()) / standard_deviation
            z_stddevs = (
                scores.groupby("tissue", sort=False)[factor].std(ddof=0)
                / standard_deviation
            )
        else:
            z_means = means * 0.0
            z_stddevs = means * 0.0

        total_ss = float(((values - values.mean()) ** 2).sum())
        group_counts = scores.groupby("tissue", sort=False)[factor].size()
        between_ss = float(
            sum(
                group_counts[tissue] * (means[tissue] - values.mean()) ** 2
                for tissue in means.index
            )
        )
        eta2 = between_ss / total_ss if total_ss > 0 else np.nan

        standardized_means.append([z_means[tissue] for tissue in tissues])
        standardized_stddevs.append([z_stddevs[tissue] for tissue in tissues])
        eta_squared.append(eta2)
        for tissue in tissues:
            records.append({
                "sample": label,
                "factor": f"Factor {number}",
                "tissue": tissue,
                "standardized_mean": z_means[tissue],
                "standardized_stddev": z_stddevs[tissue],
                "eta_squared": eta2,
                "n_cells": int(group_counts[tissue]),
            })

    heatmap = pd.DataFrame(
        standardized_means,
        index=[f"Factor {number}" for number, _factor in active],
        columns=tissues,
    )
    stddevs = pd.DataFrame(
        standardized_stddevs,
        index=heatmap.index,
        columns=tissues,
    )
    return heatmap, stddevs, np.asarray(eta_squared), records


def main():
    args = parse_args()
    summaries = [
        (
            label,
            *summarize_sample(
                label, scores_path, summary_path,
                args.active_variance_threshold,
            ),
        )
        for label, scores_path, summary_path in args.sample
    ]
    records = [
        record
        for _label, _means, _stddevs, _eta_squared, sample_records in summaries
        for record in sample_records
    ]
    pd.DataFrame(records).to_csv(args.table_output, sep="\t", index=False)

    sns.set_theme(style="white", context="paper", font_scale=1.15)
    mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})

    tissues = [
        tissue
        for tissue in TISSUE_ORDER
        if any(tissue in means.columns for _label, means, *_rest in summaries)
    ]
    factor_order = list(summaries[0][1].index)
    y_positions = np.arange(len(tissues))
    x_limit = max(
        float(
            np.nanmax(
                np.abs(means[tissues].to_numpy(float))
                + stddevs[tissues].to_numpy(float)
            )
        )
        for _label, means, stddevs, _eta_squared, _records in summaries
    )
    x_limit = max(1.0, np.ceil(2.0 * x_limit) / 2.0)

    figure, axes = plt.subplots(
        len(summaries),
        len(factor_order),
        figsize=(3.35 * len(factor_order), 2.65 * len(summaries)),
        sharex=True,
        sharey=True,
        squeeze=False,
    )

    for row, (label, means, stddevs, _eta_squared, _records) in enumerate(summaries):
        for column, factor in enumerate(factor_order):
            axis = axes[row, column]
            for tissue_index, tissue in enumerate(tissues):
                axis.errorbar(
                    means.loc[factor, tissue],
                    tissue_index,
                    xerr=stddevs.loc[factor, tissue],
                    fmt="o",
                    markersize=7.0,
                    color=TISSUE_COLORS.get(tissue, "0.35"),
                    ecolor=TISSUE_COLORS.get(tissue, "0.35"),
                    elinewidth=1.7,
                    capsize=3.2,
                    capthick=1.7,
                    markeredgecolor="white",
                    markeredgewidth=0.6,
                    zorder=2,
                )
            axis.set_xlim(-x_limit, x_limit)
            axis.set_ylim(len(tissues) - 0.5, -0.5)
            axis.set_yticks(y_positions)
            axis.set_yticklabels(tissues)
            axis.grid(axis="x", color="0.92", linewidth=0.7)
            axis.set_axisbelow(True)
            axis.spines[["top", "right"]].set_visible(False)

            if row == 0:
                axis.set_title(factor, fontweight="normal", pad=7)
            if column == 0:
                axis.set_ylabel(label, fontweight="normal", labelpad=12)

    tissue_handles = [
        mpl.lines.Line2D(
            [0], [0],
            marker="o",
            linestyle="-",
            color=TISSUE_COLORS[tissue],
            markerfacecolor=TISSUE_COLORS[tissue],
            markeredgecolor="white",
            linewidth=1.4,
            markersize=6,
            label=tissue,
        )
        for tissue in tissues
    ]
    figure.legend(
        handles=tissue_handles,
        loc="center left",
        bbox_to_anchor=(0.91, 0.5),
        ncol=1,
        frameon=False,
        title="Tissue",
    )
    figure.subplots_adjust(
        top=0.90,
        bottom=0.14,
        left=0.08,
        right=0.90,
        hspace=0.20,
        wspace=0.16,
    )
    figure.savefig(args.output, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
