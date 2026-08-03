#!/usr/bin/env python3

import argparse
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, confusion_matrix
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


TISSUE_ORDER = ["M1", "M2", "RW", "RE", "LL", "Liv"]
SAMPLE_COLORS = ["#4c78a8", "#f58518"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Cross-validated prediction of tissue from active factor scores."
    )
    parser.add_argument(
        "--sample",
        nargs=3,
        action="append",
        required=True,
        metavar=("LABEL", "FACTOR_SCORES", "FIT_SUMMARY"),
    )
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--active-variance-threshold", type=float, default=1e-5)
    parser.add_argument("--output", required=True)
    parser.add_argument("--metrics-output", required=True)
    parser.add_argument("--predictions-output", required=True)
    return parser.parse_args()


def active_factors(summary_path, threshold):
    summary = pd.read_table(summary_path, index_col=0)
    factors = []
    for parameter, value in summary["value"].items():
        match = re.fullmatch(r"sigma2_latent_LF([0-9]+)", parameter)
        if match and float(value) > threshold:
            factors.append(f"factor_{int(match.group(1))}")
    if not factors:
        raise ValueError(f"No active factors found in {summary_path}")
    return factors


def classify_sample(label, scores_path, summary_path, args):
    scores = pd.read_table(scores_path, index_col=0)
    factors = active_factors(summary_path, args.active_variance_threshold)
    missing = [factor for factor in factors if factor not in scores.columns]
    if missing:
        raise ValueError(
            f"Active factors absent from {scores_path}: {', '.join(missing)}"
        )

    tissues = scores.index.astype(str).str.split(".", n=1).str[0]
    observed_tissues = [
        tissue for tissue in TISSUE_ORDER if tissue in set(tissues)
    ]
    X = scores[factors].to_numpy(float)
    y = tissues.to_numpy(str)

    model = make_pipeline(
        StandardScaler(),
        LogisticRegression(
            solver="lbfgs",
            class_weight="balanced",
            max_iter=2000,
            random_state=args.seed,
        ),
    )
    folds = StratifiedKFold(
        n_splits=args.folds,
        shuffle=True,
        random_state=args.seed,
    )
    predicted = cross_val_predict(model, X, y, cv=folds, method="predict")

    matrix = confusion_matrix(y, predicted, labels=observed_tissues)
    recalls = np.divide(
        np.diag(matrix),
        matrix.sum(axis=1),
        out=np.zeros(len(observed_tissues), dtype=float),
        where=matrix.sum(axis=1) > 0,
    )
    balanced_accuracy = balanced_accuracy_score(y, predicted)

    metrics = [{
        "sample": label,
        "metric": "Overall",
        "accuracy": balanced_accuracy,
        "n_cells": len(y),
        "n_active_factors": len(factors),
    }]
    metrics.extend({
        "sample": label,
        "metric": tissue,
        "accuracy": recall,
        "n_cells": int((y == tissue).sum()),
        "n_active_factors": len(factors),
    } for tissue, recall in zip(observed_tissues, recalls))

    predictions = pd.DataFrame({
        "sample": label,
        "cell": scores.index.astype(str),
        "observed_tissue": y,
        "predicted_tissue": predicted,
        "correct": y == predicted,
    })
    return metrics, predictions, len(observed_tissues)


def main():
    args = parse_args()
    all_metrics = []
    all_predictions = []
    class_counts = []
    for label, scores_path, summary_path in args.sample:
        metrics, predictions, n_classes = classify_sample(
            label, scores_path, summary_path, args
        )
        all_metrics.extend(metrics)
        all_predictions.append(predictions)
        class_counts.append(n_classes)

    metrics = pd.DataFrame(all_metrics)
    predictions = pd.concat(all_predictions, ignore_index=True)
    metrics.to_csv(args.metrics_output, sep="\t", index=False)
    predictions.to_csv(args.predictions_output, sep="\t", index=False)

    metric_order = ["Overall"] + [
        tissue for tissue in TISSUE_ORDER
        if tissue in set(metrics["metric"])
    ]
    y_positions = np.arange(len(metric_order))
    offsets = np.linspace(-0.10, 0.10, len(args.sample))

    sns.set_theme(style="white", context="paper", font_scale=1.15)
    mpl.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42})
    figure, axis = plt.subplots(figsize=(6.5, 4.8))

    for index, (label, _scores_path, _summary_path) in enumerate(args.sample):
        sample_metrics = metrics[metrics["sample"] == label].set_index("metric")
        values = np.array([
            100.0 * sample_metrics.loc[metric, "accuracy"]
            for metric in metric_order
        ])
        axis.scatter(
            values,
            y_positions + offsets[index],
            s=62,
            color=SAMPLE_COLORS[index],
            edgecolor="white",
            linewidth=0.7,
            label=label,
            zorder=3,
        )

    chance = 100.0 / class_counts[0]
    axis.axvline(
        chance,
        color="0.55",
        linestyle="--",
        linewidth=1.0,
        label=f"Chance ({chance:.0f}%)",
        zorder=1,
    )
    axis.set_xlim(0, 100)
    axis.set_ylim(len(metric_order) - 0.5, -0.5)
    axis.set_yticks(y_positions)
    axis.set_yticklabels(metric_order)
    axis.set_xlabel("Cross-validated tissue prediction accuracy")
    axis.set_ylabel("")
    axis.xaxis.set_major_formatter(mpl.ticker.PercentFormatter(xmax=100))
    axis.grid(axis="x", color="0.9", linewidth=0.8)
    axis.set_axisbelow(True)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False, loc="lower right")
    figure.tight_layout()
    figure.savefig(args.output, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
