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

performance_tsv = sys.argv[1]
time_tsv = sys.argv[2]
output_pdf = sys.argv[3]

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

df = pd.read_csv(performance_tsv, sep="\t")

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

time_df = pd.read_csv(time_tsv, sep="\t")

required_time_cols = {"ntaxa", "method", "sim_num", "time_sec"}
missing_time_cols = required_time_cols - set(time_df.columns)
if missing_time_cols:
    sys.stderr.write(f"Missing required time columns: {sorted(missing_time_cols)}\n")
    sys.exit(1)

time_df["method"] = pd.Categorical(
    time_df["method"], categories=method_order, ordered=True
)
time_df["method_label"] = time_df["method"].map(method_labels)

ntaxa_order = sorted(df["ntaxa"].unique())
method_label_order = [method_labels[m] for m in method_order]
palette_by_label = {method_labels[m]: method_colors[m] for m in method_order}

fig, axes = plt.subplots(1, 3, figsize=(12.25, 3.5))
plt.subplots_adjust(right=0.86, wspace=0.35)

for ax, metric, ylabel in zip(axes, ["precision", "recall"], ["Precision", "Recall"]):
    sub = plot_df.loc[plot_df["metric"] == metric].copy()

    sns.barplot(
        data=sub,
        x="ntaxa",
        y="value",
        hue="method",
        order=ntaxa_order,
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

runtime_ax = axes[2]
sns.barplot(
    data=time_df,
    x="ntaxa",
    y="time_sec",
    hue="method_label",
    order=ntaxa_order,
    hue_order=method_label_order,
    palette=palette_by_label,
    errorbar="sd",
    capsize=0.1,
    err_kws={"linewidth": 1},
    ax=runtime_ax,
)

runtime_ax.set_xlabel("Number of taxa")
for label in runtime_ax.get_xticklabels():
    label.set_rotation(-30)
    label.set_horizontalalignment("center")
runtime_ax.set_ylabel("Runtime (s)")
runtime_ax.set_title("")
runtime_ax.set_yscale("log")
runtime_ax.spines["top"].set_visible(False)
runtime_ax.spines["right"].set_visible(False)
runtime_ax.grid(axis="y", linewidth=0.5, alpha=0.3)
runtime_ax.set_axisbelow(True)
runtime_ax.legend_.remove()

handles, labels = axes[0].get_legend_handles_labels()
labels = [method_labels[m] for m in method_order]
fig.legend(
    handles,
    labels,
    loc="center left",
    frameon=False,
    bbox_to_anchor=(0.87, 0.5),
)

fig.savefig(output_pdf, bbox_inches="tight")
plt.close(fig)

print(f"Saved figure to: {output_pdf}")
