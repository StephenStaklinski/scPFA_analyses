#!/usr/bin/env python3

import argparse
import re

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot fitted Brownian variance for each latent factor."
    )
    parser.add_argument("fit_summary_tsv")
    parser.add_argument("output_pdf")
    return parser.parse_args()


def main():
    args = parse_args()
    summary = pd.read_table(args.fit_summary_tsv, index_col=0)

    variances = []
    for parameter, value in summary["value"].items():
        match = re.fullmatch(r"sigma2_latent_LF([0-9]+)", parameter)
        if match:
            factor = int(match.group(1))
            variances.append((factor, float(value)))
    variances.sort()
    if not variances:
        raise ValueError(
            f"No sigma2_latent_LF parameters found in {args.fit_summary_tsv}"
        )

    sns.set_theme(style="ticks", context="paper")
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    labels = [f"Factor {factor}" for factor, _value in variances]
    values = [value for _factor, value in variances]
    figure, axis = plt.subplots(figsize=(4.5, 3.5))
    axis.bar(labels, values, color=sns.color_palette("deep")[0])
    axis.set_xlabel("Factors")
    axis.set_ylabel(r"Brownian variance ($\sigma^2$)")
    axis.tick_params(axis="x", rotation=90)
    sns.despine(ax=axis)
    figure.tight_layout()
    figure.savefig(args.output_pdf, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
