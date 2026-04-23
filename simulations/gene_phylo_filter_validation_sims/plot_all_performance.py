#!/usr/bin/env python3

import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set consistent styling for figures
plt.rcParams.update(
    {
        "font.size": 10,
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

input_tsv = sys.argv[1]
output_pdf = sys.argv[2]

method_order = ["path", "moran", "lrt_full", "lrt_lambda"]
method_labels = {
    "path": "PATH",
    "moran": "PATH re-implementation",
    "lrt_full": "LRT full covariance",
    "lrt_lambda": "LRT Pagel's lambda covariance",
}

# Colorblind-safe palette
method_colors = {
    "path": "#F0E442",
    "moran": "#E69F00",
    "lrt_full": "#56B4E9",
    "lrt_lambda": "#009E73",
}

df = pd.read_csv(input_tsv, sep="\t")

required_cols = {"ntaxa", "method", "sim_num", "TP", "FN", "TN", "FP"}
missing_cols = required_cols - set(df.columns)
if missing_cols:
    sys.stderr.write(f"Missing required columns: {sorted(missing_cols)}\n")
    sys.exit(1)

tp = df["TP"].to_numpy(dtype=float)
fn = df["FN"].to_numpy(dtype=float)
fp = df["FP"].to_numpy(dtype=float)

precision = np.full(len(df), np.nan, dtype=float)
recall = np.full(len(df), np.nan, dtype=float)

precision_denom = tp + fp
recall_denom = tp + fn

precision[precision_denom != 0] = tp[precision_denom != 0] / precision_denom[precision_denom != 0]
recall[recall_denom != 0] = tp[recall_denom != 0] / recall_denom[recall_denom != 0]

df["precision"] = precision
df["recall"] = recall
df["method"] = pd.Categorical(df["method"], categories=method_order, ordered=True)

plot_df = df.melt(
    id_vars=["ntaxa", "method", "sim_num"],
    value_vars=["precision", "recall"],
    var_name="metric",
    value_name="value",
)

plot_df["metric"] = pd.Categorical(
    plot_df["metric"],
    categories=["precision", "recall"],
    ordered=True,
)

plot_df["method_label"] = plot_df["method"].map(method_labels)

fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.5))
plt.subplots_adjust(right=0.80)

for ax, metric, ylabel in zip(axes, ["precision", "recall"], ["Precision", "Recall"]):
    sub = plot_df.loc[plot_df["metric"] == metric].copy()

    sns.barplot(
        data=sub,
        x="ntaxa",
        y="value",
        hue="method",
        order=sorted(df["ntaxa"].unique()),
        hue_order=method_order,
        palette=method_colors,
        errorbar="sd",
        capsize=0.1,
        err_kws={"linewidth": 1},
        ax=ax,
    )

    ax.set_xlabel("Number of taxa")
    for label in ax.get_xticklabels():
        label.set_rotation(-30)
        label.set_horizontalalignment("center")
    ax.set_ylabel(ylabel)
    ax.set_title("")
    ax.set_ylim(0, 1.05)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linewidth=0.5, alpha=0.3)
    ax.set_axisbelow(True)
    ax.legend_.remove()

handles, labels = axes[0].get_legend_handles_labels()
labels = [method_labels[m] for m in method_order]
fig.legend(
    handles,
    labels,
    loc="center left",
    frameon=False,
    bbox_to_anchor=(0.82, 0.5),
)

fig.savefig(output_pdf, bbox_inches="tight")
plt.close(fig)

print(f"Saved figure to: {output_pdf}")