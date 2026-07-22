#!/usr/bin/env python3

import argparse
import os
import sys

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
import seaborn as sns


sns.set_theme(style="white", context="paper")

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Make a penalty-sweep heatmap focused on factor overlap."
    )
    parser.add_argument("summary_tsv")
    parser.add_argument("output_pdf")
    parser.add_argument("--clone", default=None)
    parser.add_argument(
        "--no-annot",
        action="store_true",
        help="Do not print metric values inside heatmap cells.",
    )
    return parser.parse_args()


def format_penalty(value):
    value = float(value)
    return "0" if value == 0 else f"{value:g}"


def metric_heatmap(
    ax,
    df,
    value_col,
    cbar_label,
    cmap,
    fmt=".3f",
    aggfunc="mean",
    annot=True,
):
    pivot = df.pivot_table(
        index="L_loading_overlap_strength",
        columns="L_l1_strength",
        values=value_col,
        aggfunc=aggfunc,
    ).sort_index(ascending=False)
    pivot = pivot.reindex(sorted(pivot.columns), axis=1)

    sns.heatmap(
        pivot,
        ax=ax,
        cmap=cmap,
        annot=annot,
        annot_kws={"fontsize": 11},
        fmt=fmt,
        linewidths=0.6,
        linecolor="white",
        cbar_kws={"label": cbar_label},
    )
    ax.set_xlabel("L1 penalty")
    ax.set_ylabel("L-overlap penalty")
    ax.set_xticklabels([format_penalty(x.get_text()) for x in ax.get_xticklabels()])
    ax.set_yticklabels([format_penalty(y.get_text()) for y in ax.get_yticklabels()], rotation=0)


def main():
    args = parse_args()
    df = pd.read_csv(args.summary_tsv, sep="\t")
    if args.clone:
        df = df[df["clone"] == args.clone].copy()
    if df.empty:
        sys.stderr.write("No rows to plot.\n")
        return 1

    required = {
        "condition",
        "remove_ribo_mito",
        "L_l1_strength",
        "L_loading_overlap_strength",
        "mean_offdiag_abs_L_pearson",
        "max_offdiag_abs_L_pearson",
        "top_gene_mean_pairwise_signed_jaccard",
        "top_gene_max_pairwise_signed_jaccard",
        "combined_log_likelihood",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        sys.stderr.write(f"Missing required columns: {missing}\n")
        return 1

    for col in [
        "L_l1_strength",
        "L_loading_overlap_strength",
        "remove_ribo_mito",
        "mean_offdiag_abs_L_pearson",
        "max_offdiag_abs_L_pearson",
        "top_gene_mean_pairwise_signed_jaccard",
        "top_gene_max_pairwise_signed_jaccard",
        "combined_log_likelihood",
    ]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=list(required - {"condition"})).copy()
    df = df[df["remove_ribo_mito"] == 0].copy()
    if df.empty:
        sys.stderr.write("No default no-ribo/mito-filter rows to plot.\n")
        return 1

    cmap = LinearSegmentedColormap.from_list(
        "low_overlap_to_high_overlap",
        [
            "#08306B",
            "#6BAED6",
            "#F7FBFF",
            "#FDD0A2",
            "#99000D",
        ],
    )

    fig, heat_axes = plt.subplots(2, 3, figsize=(18.0, 11.2), constrained_layout=True)

    metric_heatmap(
        heat_axes[0, 0],
        df,
        "mean_offdiag_abs_L_pearson",
        "Mean pairwise |Pearson r|\nbetween L rows",
        cmap,
        annot=not args.no_annot,
    )
    metric_heatmap(
        heat_axes[0, 1],
        df,
        "top_gene_mean_pairwise_signed_jaccard",
        "Mean pairwise signed Jaccard\ntop-gene overlap",
        cmap,
        annot=not args.no_annot,
    )
    metric_heatmap(
        heat_axes[0, 2],
        df,
        "combined_log_likelihood",
        "Best-state\ncombined log likelihood",
        cmap.reversed(),
        fmt=".0f",
        annot=not args.no_annot,
    )

    metric_heatmap(
        heat_axes[1, 0],
        df,
        "max_offdiag_abs_L_pearson",
        "Max pairwise |Pearson r|\nbetween L rows",
        cmap,
        annot=not args.no_annot,
    )
    metric_heatmap(
        heat_axes[1, 1],
        df,
        "top_gene_max_pairwise_signed_jaccard",
        "Max pairwise signed Jaccard\ntop-gene overlap",
        cmap,
        annot=not args.no_annot,
    )
    heat_axes[1, 2].set_visible(False)

    fig.savefig(args.output_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved figure to: {args.output_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
