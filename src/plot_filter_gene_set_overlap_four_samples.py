#!/usr/bin/env python3
"""Plot gene recovery and heritability overlap across four samples."""

import argparse

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from plot_filter_gene_set_overlap import read_results


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_tsvs", nargs=4)
    parser.add_argument("output_prefix")
    parser.add_argument(
        "--sample-names",
        nargs=4,
        default=["CP1", "CP2", "CP3", "CP4"],
    )
    return parser.parse_args()


def main():
    args = parse_args()
    samples = [
        read_results(path) for path in args.sample_tsvs
    ]
    names = args.sample_names

    union = set().union(*(set(sample.index) for sample in samples))
    shared = set(samples[0].index).intersection(
        *(set(sample.index) for sample in samples[1:])
    )

    recovery_counts = {number: 0 for number in range(1, 5)}
    for gene in union:
        number_present = sum(gene in sample.index for sample in samples)
        recovery_counts[number_present] += 1

    heritability_counts = {number: 0 for number in range(5)}
    for gene in shared:
        number_heritable = sum(bool(sample.loc[gene]) for sample in samples)
        heritability_counts[number_heritable] += 1

    sky_blue = "#56B4E9"
    teal = "#009E73"
    amber = "#E69F00"
    pink = "#CC79A7"
    blue = "#0072B2"
    vermillion = "#D55E00"
    gray = "#BDBDBD"

    bar_data = []
    for name, sample in zip(names, samples):
        heritable = int(sample.sum())
        total = len(sample)
        bar_data.append(
            (
                name,
                total,
                [
                    ("Heritable", heritable, teal),
                    ("Not heritable", total - heritable, gray),
                ],
            )
        )

    bar_data.extend(
        [
            (
                "Genes recovered across samples",
                len(union),
                [
                    ("Recovered in 1 sample", recovery_counts[1], sky_blue),
                    ("Recovered in 2 samples", recovery_counts[2], amber),
                    ("Recovered in 3 samples", recovery_counts[3], blue),
                    ("Recovered in all 4 samples", recovery_counts[4], pink),
                ],
            ),
            (
                "Heritability of genes recovered in all samples",
                len(shared),
                [
                    ("Heritable in all 4 samples", heritability_counts[4], teal),
                    ("Heritable in 3 samples", heritability_counts[3], blue),
                    ("Heritable in 2 samples", heritability_counts[2], pink),
                    ("Heritable in 1 sample", heritability_counts[1], vermillion),
                    ("Not heritable in any sample", heritability_counts[0], gray),
                ],
            ),
        ]
    )

    figure, axes = plt.subplots(
        len(bar_data),
        1,
        figsize=(16, 12),
        sharex=False,
    )
    for axis, (bar_name, total, segments) in zip(axes, bar_data):
        left = 0.0
        for _, count, color in segments:
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
        count_ticks = [
            round(total * fraction) for fraction in tick_fractions
        ]
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

    percent_axis = figure.add_axes([0.20, 0.045, 0.59, 0.001])
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
        bottom=0.11,
        top=0.97,
        hspace=1.05,
    )
    figure.savefig(f"{args.output_prefix}.pdf", bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
