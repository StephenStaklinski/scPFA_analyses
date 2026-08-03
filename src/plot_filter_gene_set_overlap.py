#!/usr/bin/env python3
"""Plot cross-sample overlap categories for genes passing a lineage filter."""

import argparse
import warnings

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_a_tsv")
    parser.add_argument("sample_b_tsv")
    parser.add_argument("output_prefix")
    parser.add_argument("--sample-a", default="TLS1")
    parser.add_argument("--sample-b", default="TLS2")
    return parser.parse_args()


def read_results(path):
    results = pd.read_table(path)
    required = {"gene", "keep"}
    missing = required - set(results.columns)
    if missing:
        raise ValueError(f"{path} is missing: {', '.join(sorted(missing))}")
    keep = results["keep"]
    if keep.dtype != bool:
        normalized = keep.astype(str).str.lower()
        if not normalized.isin({"true", "false"}).all():
            raise ValueError(f"{path} contains non-Boolean keep values")
        keep = normalized.eq("true")

    gene_results = pd.DataFrame(
        {
            "gene": results["gene"].astype(str),
            "keep": keep.to_numpy(),
        }
    )
    duplicate_rows = int(gene_results["gene"].duplicated().sum())
    if duplicate_rows:
        conflicting_symbols = int(
            (gene_results.groupby("gene", sort=False)["keep"].nunique() > 1).sum()
        )
        warnings.warn(
            f"{path}: collapsed {duplicate_rows} duplicate row(s) by gene symbol "
            "using keep=any; "
            f"{conflicting_symbols} symbol(s) had conflicting keep values"
        )

    return gene_results.groupby("gene", sort=False)["keep"].any()


def main():
    args = parse_args()
    sample_a = read_results(args.sample_a_tsv)
    sample_b = read_results(args.sample_b_tsv)

    genes_a = set(sample_a.index)
    genes_b = set(sample_b.index)
    shared = genes_a & genes_b
    only_a = genes_a - genes_b
    only_b = genes_b - genes_a

    shared_a = sample_a.loc[list(shared)]
    shared_b = sample_b.loc[list(shared)]
    both_failing = int((~shared_a & ~shared_b).sum())
    only_a_failing = int((~sample_a.loc[list(only_a)]).sum())
    only_b_failing = int((~sample_b.loc[list(only_b)]).sum())

    rows = [
        {
            "category": "heritable_a_not_present_b",
            "description": (
                f"Heritable in {args.sample_a}; not present in {args.sample_b}"
            ),
            "count": int(sample_a.loc[list(only_a)].sum()),
        },
        {
            "category": "heritable_a_not_heritable_b",
            "description": (
                f"Heritable in {args.sample_a}; not heritable in {args.sample_b}"
            ),
            "count": int((shared_a & ~shared_b).sum()),
        },
        {
            "category": "heritable_both",
            "description": f"Heritable in both {args.sample_a} and {args.sample_b}",
            "count": int((shared_a & shared_b).sum()),
        },
        {
            "category": "heritable_b_not_heritable_a",
            "description": (
                f"Heritable in {args.sample_b}; not heritable in {args.sample_a}"
            ),
            "count": int((~shared_a & shared_b).sum()),
        },
        {
            "category": "heritable_b_not_present_a",
            "description": (
                f"Heritable in {args.sample_b}; not present in {args.sample_a}"
            ),
            "count": int(sample_b.loc[list(only_b)].sum()),
        },
        {
            "category": "not_heritable_a_not_present_b",
            "description": (
                f"Not heritable in {args.sample_a}; not present in {args.sample_b}"
            ),
            "count": only_a_failing,
        },
        {
            "category": "not_heritable_both",
            "description": (
                f"Not heritable in either {args.sample_a} or {args.sample_b}"
            ),
            "count": both_failing,
        },
        {
            "category": "not_heritable_b_not_present_a",
            "description": (
                f"Not heritable in {args.sample_b}; not present in {args.sample_a}"
            ),
            "count": only_b_failing,
        },
    ]
    counts = pd.DataFrame(rows)
    counts.to_csv(f"{args.output_prefix}.tsv", sep="\t", index=False)

    sky_blue = "#56B4E9"
    teal = "#009E73"
    amber = "#E69F00"
    pink = "#CC79A7"
    gray = "#BDBDBD"

    tls1_only_total = rows[0]["count"] + rows[5]["count"]
    tls2_only_total = rows[4]["count"] + rows[7]["count"]
    shared_total = (
        rows[1]["count"]
        + rows[2]["count"]
        + rows[3]["count"]
        + rows[6]["count"]
    )
    union_total = tls1_only_total + shared_total + tls2_only_total
    sample_a_heritable = int(sample_a.sum())
    sample_b_heritable = int(sample_b.sum())
    sample_a_total = len(sample_a)
    sample_b_total = len(sample_b)

    bar_data = [
        (
            args.sample_a,
            sample_a_total,
            [
                ("Heritable", sample_a_heritable, teal),
                ("Not heritable", sample_a_total - sample_a_heritable, gray),
            ],
        ),
        (
            args.sample_b,
            sample_b_total,
            [
                ("Heritable", sample_b_heritable, teal),
                ("Not heritable", sample_b_total - sample_b_heritable, gray),
            ],
        ),
        (
            "Genes recovered across samples",
            union_total,
            [
                (f"{args.sample_a} only", tls1_only_total, sky_blue),
                ("Present in both", shared_total, pink),
                (f"{args.sample_b} only", tls2_only_total, amber),
            ],
        ),
        (
            "Heritability of genes recovered in both samples",
            shared_total,
            [
                (
                    f"Heritable in {args.sample_a} only",
                    rows[1]["count"],
                    sky_blue,
                ),
                ("Heritable in both", rows[2]["count"], pink),
                (
                    f"Heritable in {args.sample_b} only",
                    rows[3]["count"],
                    amber,
                ),
                ("Not heritable in either", rows[6]["count"], gray),
            ],
        ),
    ]

    figure, axes = plt.subplots(
        len(bar_data),
        1,
        figsize=(16, 8.5),
        sharex=False,
    )
    for axis, (bar_name, total, segments) in zip(axes, bar_data):
        left = 0.0
        for label, count, color in segments:
            fraction = count / total
            axis.barh(
                0,
                fraction,
                left=left,
                height=0.62,
                color=color,
                edgecolor="white",
                linewidth=1.2,
            )
            left += fraction

        axis.set_xlim(0, 1)
        axis.set_ylim(-0.5, 0.5)
        axis.set_yticks([])
        axis.set_title(bar_name, loc="center", fontsize=14, pad=-2)
        axis.spines[["top", "right", "left"]].set_visible(False)
        tick_fractions = (0, 0.25, 0.5, 0.75, 1)
        count_ticks = [round(total * fraction) for fraction in tick_fractions]
        axis.set_xticks(tick_fractions)
        axis.set_xticklabels([f"{value:,}" for value in count_ticks])
        axis.set_xlabel("Gene count", fontsize=11, labelpad=4)
        axis.tick_params(axis="x", labelsize=10, pad=3)
        axis.spines["bottom"].set_color("0.45")
        axis.legend(
            handles=[
                Patch(facecolor=color, edgecolor="none", label=label)
                for label, _, color in segments
            ],
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            frameon=False,
            fontsize=10,
            handlelength=1.3,
            handleheight=1.0,
            labelspacing=0.5,
            borderaxespad=0,
        )

    percent_axis = figure.add_axes(
        [0.20, 0.055, 0.59, 0.001],
    )
    percent_axis.set_xlim(0, 1)
    percent_axis.set_ylim(0, 1)
    percent_axis.set_yticks([])
    percent_axis.spines[["top", "right", "left"]].set_visible(False)
    percent_axis.spines["bottom"].set_color("0.25")
    percent_axis.spines["bottom"].set_linewidth(0.8)
    percent_axis.patch.set_alpha(0)
    percent_axis.set_xticks([0, 0.25, 0.5, 0.75, 1])
    percent_axis.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
    percent_axis.tick_params(axis="x", length=4, labelsize=12)
    figure.subplots_adjust(
        left=0.20,
        right=0.79,
        bottom=0.13,
        top=0.96,
        hspace=1.05,
    )

    figure.savefig(
        f"{args.output_prefix}.pdf",
        bbox_inches="tight",
    )
    plt.close(figure)


if __name__ == "__main__":
    main()
