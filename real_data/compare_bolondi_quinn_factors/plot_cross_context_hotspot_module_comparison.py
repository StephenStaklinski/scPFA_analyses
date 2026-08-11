#!/usr/bin/env python3

"""Compare Hotspot modules across development and cancer using saved GO results."""

import argparse
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
        nargs=3,
        required=True,
        metavar=("LABEL", "CONTEXT", "GO_ENRICHMENT"),
    )
    parser.add_argument("--adjusted-p-threshold", type=float, default=0.05)
    parser.add_argument("--out-prefix", required=True)
    return parser.parse_args()


def read_sample(specification, p_threshold):
    label, context, enrichment_path = specification
    enrichment = pd.read_table(enrichment_path)
    required = {"module", "ID", "p.adjust"}
    missing = required.difference(enrichment.columns)
    if missing:
        raise ValueError(f"Missing columns in {enrichment_path}: {sorted(missing)}")
    enrichment = enrichment.loc[enrichment["p.adjust"] <= p_threshold].copy()
    enrichment["module"] = enrichment["module"].astype(int)
    enrichment["module_label"] = (
        label + " M" + enrichment["module"].astype(str)
    )
    module_labels = [
        f"{label} M{number}" for number in sorted(enrichment["module"].unique())
    ]
    return {
        "context": context,
        "module_labels": module_labels,
        "enrichment": enrichment,
    }


def go_term_sets(samples):
    result = {}
    for sample in samples:
        for module_label in sample["module_labels"]:
            result[module_label] = set(
                sample["enrichment"].loc[
                    sample["enrichment"]["module_label"] == module_label, "ID"
                ].dropna()
            )
    return result


def jaccard(left, right):
    union = left | right
    return len(left & right) / len(union) if union else 0.0


def jointly_order(similarity):
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


def cross_context_similarity(samples, term_sets):
    contexts = []
    for sample in samples:
        if sample["context"] not in contexts:
            contexts.append(sample["context"])
    if len(contexts) != 2:
        raise ValueError("Exactly two contexts are required")
    rows = [
        label for sample in samples if sample["context"] == contexts[0]
        for label in sample["module_labels"]
    ]
    columns = [
        label for sample in samples if sample["context"] == contexts[1]
        for label in sample["module_labels"]
    ]
    similarity = pd.DataFrame(
        [[jaccard(term_sets[row], term_sets[column]) for column in columns]
         for row in rows],
        index=rows,
        columns=columns,
    )
    return jointly_order(similarity)


def plot_similarity(similarity, output_path):
    figure, axis = plt.subplots(figsize=(8.2, 7.0))
    rd_bu = plt.get_cmap("RdBu_r")
    positive_rd_bu = LinearSegmentedColormap.from_list(
        "positive_RdBu_r", [rd_bu(value) for value in [0.5, 0.625, 0.75, 0.875, 1.0]]
    )
    sns.heatmap(
        similarity,
        ax=axis,
        cmap=positive_rd_bu,
        vmin=0,
        vmax=0.15,
        linewidths=0.45,
        linecolor="white",
        cbar_kws={
            "label": "Jaccard similarity of significant GO term sets",
            "shrink": 0.78,
            "ticks": [0.00, 0.05, 0.10, 0.15],
        },
    )
    axis.set_xlabel("Lung cancer Hotspot modules")
    axis.set_ylabel("Developmental Hotspot modules")
    axis.tick_params(axis="x", rotation=45)
    axis.tick_params(axis="y", rotation=0)
    figure.tight_layout()
    figure.savefig(output_path, bbox_inches="tight")
    plt.close(figure)


def main():
    args = parse_args()
    samples = [read_sample(spec, args.adjusted_p_threshold) for spec in args.sample]
    sns.set_theme(style="white", context="paper", font_scale=1.08)
    similarity = cross_context_similarity(samples, go_term_sets(samples))
    prefix = Path(args.out_prefix)
    similarity.to_csv(
        f"{prefix}.cross_context_hotspot_module_similarity.tsv",
        sep="\t",
        index_label="module",
    )
    plot_similarity(
        similarity, f"{prefix}.cross_context_hotspot_module_similarity.pdf"
    )


if __name__ == "__main__":
    main()
