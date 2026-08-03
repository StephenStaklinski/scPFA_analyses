#!/usr/bin/env python3
"""Plot heritability of ribosomal and mitochondrial genes."""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Patch


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("output_pdf")
    parser.add_argument(
        "--sample",
        action="append",
        nargs=3,
        metavar=("LABEL", "ALL_GENES_TSV", "NO_RIBO_MITO_TSV"),
        required=True,
    )
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

    return pd.DataFrame(
        {"gene": results["gene"].astype(str), "keep": keep.to_numpy()}
    ).groupby("gene", sort=False)["keep"].any()


def main():
    args = parse_args()
    teal = "#009E73"
    gray = "#BDBDBD"
    figure, axes = plt.subplots(
        2,
        len(args.sample),
        figsize=(3.5 * len(args.sample) + 1.4, 6.2),
        squeeze=False,
    )

    for sample_index, (sample, all_path, retained_path) in enumerate(args.sample):
        all_genes = read_results(all_path)
        retained_genes = set(read_results(retained_path).index)
        removed = all_genes.loc[~all_genes.index.isin(retained_genes)]
        upper_names = removed.index.str.upper()
        categories = [
            ("Ribosomal", removed.loc[~upper_names.str.startswith("MT")]),
            ("Mitochondrial", removed.loc[upper_names.str.startswith("MT")]),
        ]

        for category_index, (title, results) in enumerate(categories):
            axis = axes[category_index, sample_index]
            heritable = int(results.sum())
            total = len(results)
            not_heritable = total - heritable
            if total == 0:
                axis.text(0.5, 0.5, "No genes", ha="center", va="center")
                axis.axis("off")
            else:
                axis.pie(
                    [heritable, not_heritable],
                    colors=[teal, gray],
                    startangle=90,
                    counterclock=False,
                    autopct=lambda percent: f"{percent:.1f}%",
                    wedgeprops={"edgecolor": "white", "linewidth": 1.2},
                    textprops={"fontsize": 11},
                )
            axis.text(
                0.5,
                0.90,
                f"n = {total:,}",
                transform=axis.transAxes,
                ha="center",
                va="bottom",
                fontsize=11,
            )
            if sample_index == 0:
                axis.set_ylabel(title, fontsize=13, labelpad=18)

        axes[0, sample_index].set_title(sample, fontsize=15, pad=14)

    figure.legend(
        handles=[
            Patch(facecolor=teal, edgecolor="none", label="Heritable"),
            Patch(facecolor=gray, edgecolor="none", label="Not heritable"),
        ],
        loc="center left",
        bbox_to_anchor=(0.845, 0.5),
        ncol=1,
        frameon=False,
        fontsize=11,
    )
    figure.subplots_adjust(
        left=0.10,
        right=0.84,
        bottom=0.08,
        top=0.90,
        hspace=-0.02,
        wspace=0.22,
    )
    figure.savefig(args.output_pdf, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
