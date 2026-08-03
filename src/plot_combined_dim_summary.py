#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import seaborn as sns

from plot_dim_log_likelihood import (
    BROWNIAN_PLOT_THRESHOLDS,
    BROWNIAN_THRESHOLD_COLORS,
    COMBINED_COLOR,
    COUNT_COLOR,
    NEW_GO_COLOR,
    OBS_COLOR,
    compute_new_go_terms,
    set_integer_k_axis,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine dimensionality-sweep diagnostics into a multi-row figure."
    )
    parser.add_argument(
        "--sample",
        action="append",
        nargs=3,
        metavar=("LABEL", "SUMMARY_TSV", "LOG_DIRECTORY"),
        required=True,
        help="Row label, dimension-summary TSV, and directory containing fit logs/GSEA.",
    )
    parser.add_argument("--gsea-padj-cutoff", type=float, default=0.05)
    parser.add_argument("output_pdf")
    return parser.parse_args()


def load_summary(summary_path, log_directory, padj_cutoff):
    df = pd.read_table(summary_path).sort_values("k").copy()
    log_directory = Path(log_directory)
    df["log_path"] = [str(log_directory / name) for name in df["log"]]
    df["n_new_go_terms"] = compute_new_go_terms(df, padj_cutoff)
    increments = df["delta_k"].fillna(df["k"])
    df["new_go_terms_per_factor"] = df["n_new_go_terms"] / increments
    return df


def plot_row(axes, df, label, padj_cutoff):
    delta = df.dropna(
        subset=["delta_observation_log_likelihood", "delta_combined_log_likelihood"]
    )
    axes[0].axhline(0, linewidth=0.8, color="0.45", linestyle=":")
    axes[0].plot(
        delta["k"],
        delta["delta_observation_log_likelihood"],
        marker="o",
        linewidth=1.4,
        color=OBS_COLOR,
        label="Observation log-likelihood",
    )
    axes[0].plot(
        delta["k"],
        delta["delta_combined_log_likelihood"],
        marker="o",
        linewidth=1.0,
        color=COMBINED_COLOR,
        label="Observation + Brownian\nlog-likelihoods",
    )
    axes[0].set_ylabel("Log-likelihood gain\nper dimension")
    axes[0].legend(frameon=False, fontsize=7)

    go = df.dropna(subset=["n_go_enriched_factors"])
    axes[1].plot(
        go["k"],
        go["n_go_enriched_factors"],
        marker="s",
        linewidth=1.4,
        color=COUNT_COLOR,
    )
    axes[1].set_ylabel(
        f"Factors with significant\nGO enrichment ($q < {padj_cutoff:g}$)"
    )
    axes[1].yaxis.set_major_locator(mticker.MaxNLocator(integer=True))

    new_go = df.dropna(subset=["new_go_terms_per_factor"])
    axes[2].plot(
        new_go["k"],
        new_go["new_go_terms_per_factor"],
        marker="s",
        linewidth=1.4,
        color=NEW_GO_COLOR,
    )
    axes[2].set_ylabel(
        f"New GO terms per\nadded factor ($q < {padj_cutoff:g}$)"
    )

    brownian_columns = [
        f"n_factors_bv_{int(round(threshold * 100))}pct"
        for threshold in BROWNIAN_PLOT_THRESHOLDS
    ]
    for threshold, column, color in zip(
        BROWNIAN_PLOT_THRESHOLDS,
        brownian_columns,
        BROWNIAN_THRESHOLD_COLORS,
    ):
        subset = df.dropna(subset=[column])
        axes[3].plot(
            subset["k"],
            subset[column],
            marker="o",
            linewidth=1.2,
            color=color,
            label=f"{int(round(threshold * 100))}%",
        )
    axes[3].set_ylabel("Factors explaining X%\nof Brownian variance")
    axes[3].yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    axes[3].legend(frameon=False, fontsize=7)

    for axis in axes:
        axis.set_xlabel("Latent dimension (k)")
        set_integer_k_axis(axis, df["k"])
        sns.despine(ax=axis)

def main():
    args = parse_args()
    if not args.sample:
        raise ValueError("At least one --sample argument is required.")

    sns.set_theme(style="ticks", context="paper")
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.labelsize": 10,
            "legend.fontsize": 8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    nrows = len(args.sample)
    figure, axes = plt.subplots(
        nrows, 4, figsize=(13.6, 3.25 * nrows), squeeze=False
    )
    for row, (label, summary_path, log_directory) in enumerate(args.sample):
        summary = load_summary(
            summary_path, log_directory, args.gsea_padj_cutoff
        )
        plot_row(axes[row], summary, label, args.gsea_padj_cutoff)

    figure.tight_layout(h_pad=3.2, w_pad=1.4, rect=(0, 0, 1, 0.98))
    for row, (label, _summary_path, _log_directory) in enumerate(args.sample):
        row_left = axes[row, 0].get_position().x0
        row_right = axes[row, -1].get_position().x1
        row_top = max(axis.get_position().y1 for axis in axes[row])
        figure.text(
            (row_left + row_right) / 2,
            row_top + 0.012,
            label,
            ha="center",
            va="bottom",
            fontsize=12,
            fontweight="normal",
        )
    output = Path(args.output_pdf)
    figure.savefig(output, bbox_inches="tight")
    plt.close(figure)
    print(f"Saved combined dimension-summary plot to: {output}")


if __name__ == "__main__":
    main()
