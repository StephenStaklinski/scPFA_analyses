#!/usr/bin/env python3

"""Compare active phylogenetic factors across development and cancer by GO GSEA."""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from scipy.optimize import linear_sum_assignment


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample",
        action="append",
        nargs=4,
        required=True,
        metavar=("LABEL", "CONTEXT", "GSEA", "SUMMARY"),
    )
    parser.add_argument("--adjusted-p-threshold", type=float, default=0.05)
    parser.add_argument("--active-variance-threshold", type=float, default=1e-5)
    parser.add_argument("--out-prefix", required=True)
    return parser.parse_args()


def active_factor_numbers(summary_path, threshold):
    summary = pd.read_table(summary_path)
    active = []
    for row in summary.itertuples(index=False):
        match = re.fullmatch(r"sigma2_latent_LF([0-9]+)", str(row.parameter))
        if match and float(row.value) > threshold:
            active.append(int(match.group(1)))
    if not active:
        raise ValueError(f"No active factors found in {summary_path}")
    return active


def read_sample(specification, p_threshold, variance_threshold):
    label, context, gsea_path, summary_path = specification
    active = active_factor_numbers(summary_path, variance_threshold)
    gsea = pd.read_table(gsea_path)
    required = {"ID", "p.adjust", "factor"}
    missing = required.difference(gsea.columns)
    if missing:
        raise ValueError(f"Missing columns in {gsea_path}: {sorted(missing)}")

    factor_numbers = gsea["factor"].str.extract(r"factor_([0-9]+)")[0].astype(int)
    gsea = gsea.loc[factor_numbers.isin(active)].copy()
    gsea["factor_number"] = factor_numbers.loc[gsea.index]
    gsea = gsea.loc[gsea["p.adjust"] <= p_threshold].copy()
    gsea["sample"] = label
    gsea["context"] = context
    gsea["factor_label"] = label + " F" + gsea["factor_number"].astype(str)

    factor_labels = [f"{label} F{number}" for number in active]
    return {
        "label": label,
        "context": context,
        "active": active,
        "factor_labels": factor_labels,
        "gsea": gsea,
    }


def factor_go_sets(samples):
    result = {}
    for sample in samples:
        for factor_label in sample["factor_labels"]:
            result[factor_label] = set(
                sample["gsea"].loc[
                    sample["gsea"]["factor_label"] == factor_label, "ID"
                ].dropna()
            )
    return result


def jaccard(left, right):
    union = left | right
    if not union:
        return 0.0
    return len(left & right) / len(union)


def cross_context_similarity(samples, go_sets):
    contexts = []
    for sample in samples:
        if sample["context"] not in contexts:
            contexts.append(sample["context"])
    if len(contexts) != 2:
        raise ValueError("Exactly two contexts are required")
    row_factors = [
        label for sample in samples if sample["context"] == contexts[0]
        for label in sample["factor_labels"]
    ]
    column_factors = [
        label for sample in samples if sample["context"] == contexts[1]
        for label in sample["factor_labels"]
    ]
    similarity = pd.DataFrame(
        [
            [jaccard(go_sets[row], go_sets[column])
             for column in column_factors]
            for row in row_factors
        ],
        index=row_factors,
        columns=column_factors,
    )
    return contexts, jointly_order(similarity)


def jointly_order(similarity):
    """Order factors around the strongest one-to-one assignments."""
    scores = similarity.fillna(0.0)
    matched_rows, matched_columns = linear_sum_assignment(-scores.to_numpy(float))
    matched = sorted(
        zip(matched_rows, matched_columns),
        key=lambda pair: scores.iat[pair[0], pair[1]],
        reverse=True,
    )
    row_order = [row for row, _column in matched]
    column_order = [column for _row, column in matched]
    unmatched_rows = [row for row in range(scores.shape[0]) if row not in row_order]
    unmatched_columns = [
        column for column in range(scores.shape[1]) if column not in column_order
    ]
    unmatched_rows.sort(key=lambda row: scores.iloc[row].max(), reverse=True)
    unmatched_columns.sort(
        key=lambda column: scores.iloc[:, column].max(), reverse=True
    )
    return similarity.iloc[row_order + unmatched_rows, column_order + unmatched_columns]


def plot_similarity(similarity, output_path):
    figure, axis = plt.subplots(figsize=(6.9, 6.1))
    # Use the positive (white-to-red) half of the same RdBu_r palette used by
    # the paper's signed correlation matrices. Jaccard similarity is unsigned.
    rd_bu = plt.get_cmap("RdBu_r")
    positive_rd_bu = LinearSegmentedColormap.from_list(
        "positive_RdBu_r", [rd_bu(value) for value in [0.5, 0.625, 0.75, 0.875, 1.0]]
    )
    sns.heatmap(
        similarity,
        ax=axis,
        cmap=positive_rd_bu,
        vmin=0,
        vmax=0.20,
        linewidths=0.45,
        linecolor="white",
        cbar_kws={
            "label": "Jaccard similarity of significant GO term sets",
            "shrink": 0.78,
            "ticks": [0.00, 0.05, 0.10, 0.15, 0.20],
        },
    )
    axis.set_xlabel("Lung cancer factors")
    axis.set_ylabel("Developmental factors")
    axis.tick_params(axis="x", rotation=45)
    axis.tick_params(axis="y", rotation=0)
    figure.tight_layout()
    figure.savefig(output_path, bbox_inches="tight")
    plt.close(figure)


def main():
    args = parse_args()
    samples = [
        read_sample(spec, args.adjusted_p_threshold, args.active_variance_threshold)
        for spec in args.sample
    ]
    sns.set_theme(style="white", context="paper", font_scale=1.08)

    go_sets = factor_go_sets(samples)
    _contexts, similarity = cross_context_similarity(samples, go_sets)

    prefix = Path(args.out_prefix)
    similarity.to_csv(
        f"{prefix}.cross_context_factor_similarity.tsv", sep="\t", index_label="factor"
    )
    plot_similarity(similarity, f"{prefix}.cross_context_factor_similarity.pdf")


if __name__ == "__main__":
    main()
