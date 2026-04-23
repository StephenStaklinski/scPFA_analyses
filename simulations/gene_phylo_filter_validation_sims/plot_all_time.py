#!/usr/bin/env python3

import sys

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

required_cols = {"ntaxa", "method", "sim_num", "time_sec"}
missing_cols = required_cols - set(df.columns)
if missing_cols:
    sys.stderr.write(f"Missing required columns: {sorted(missing_cols)}\n")
    sys.exit(1)

df["method"] = pd.Categorical(df["method"], categories=method_order, ordered=True)
df["method_label"] = df["method"].map(method_labels)

ntaxa_order = sorted(df["ntaxa"].unique())
method_label_order = [method_labels[m] for m in method_order]
palette_by_label = {method_labels[m]: method_colors[m] for m in method_order}

fig, ax = plt.subplots(1, 1, figsize=(4.75, 3.5))
plt.subplots_adjust(right=0.80)

sns.barplot(
    data=df,
    x="ntaxa",
    y="time_sec",
    hue="method_label",
    order=ntaxa_order,
    hue_order=method_label_order,
    palette=palette_by_label,
    errorbar="sd",
    capsize=0.1,
    err_kws={"linewidth": 1},
    ax=ax,
)

ax.set_xlabel("Number of taxa")
for label in ax.get_xticklabels():
        label.set_rotation(-30)
        label.set_horizontalalignment("center")
ax.set_ylabel("Runtime (s)")
ax.set_title("")
ax.set_yscale("log")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="y", linewidth=0.5, alpha=0.3)
ax.set_axisbelow(True)

handles, labels = ax.get_legend_handles_labels()
ax.legend_.remove()
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