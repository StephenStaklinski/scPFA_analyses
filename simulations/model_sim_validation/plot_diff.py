#!/usr/bin/env python3

import sys

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.labelsize": 12,
        "axes.titlesize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

input_tsv = sys.argv[1]
output_pdf = sys.argv[2]

df = pd.read_csv(input_tsv, sep="\t")

required_cols = {"NTAXA", "NGENES", "sim"}
missing_cols = required_cols - set(df.columns)
if missing_cols:
    sys.stderr.write(f"Missing required columns: {sorted(missing_cols)}\n")
    sys.exit(1)

# Make sure grouping columns are numeric
df["NTAXA"] = pd.to_numeric(df["NTAXA"])
df["NGENES"] = pd.to_numeric(df["NGENES"])

# Metrics = everything except grouping / replicate columns
exclude_cols = {"NTAXA", "NGENES", "sim"}
metric_cols = [c for c in df.columns if c not in exclude_cols]

if len(metric_cols) == 0:
    sys.stderr.write("No metric columns found to plot.\n")
    sys.exit(1)

# Convert metrics to numeric
for col in metric_cols:
    df[col] = pd.to_numeric(df[col], errors="coerce")

ntaxa_order = sorted(df["NTAXA"].dropna().unique())
ngenes_order = sorted(df["NGENES"].dropna().unique())

# Same style of palette as before
ngenes_colors = {
    500: "#E69F00",
    1000: "#56B4E9",
    5000: "#009E73",
}
fallback_colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", "#D55E00"]

palette_by_ngenes = {}
for i, ngenes in enumerate(ngenes_order):
    palette_by_ngenes[ngenes] = ngenes_colors.get(ngenes, fallback_colors[i % len(fallback_colors)])

metric_labels = {
    "z_rmse": "Z RMSE",
    "x_rmse": "X RMSE",
    "cell_cov_correlation": "Cell covariance correlation",
    "gene_cov_correlation": "Gene covariance correlation",
    "F_col_mean_abs_corr": "F column mean abs corr",
    "L_row_mean_abs_corr": "L row mean abs corr",
}

fig, axes = plt.subplots(2, 3, figsize=(4.5 * 3, 3.6 * 2))

for i, metric in enumerate(metric_cols):
    ax = axes[i // 3, i % 3]

    plot_df = df[["NTAXA", "NGENES", metric]].dropna()

    sns.barplot(
        data=plot_df,
        x="NTAXA",
        y=metric,
        hue="NGENES",
        order=ntaxa_order,
        hue_order=ngenes_order,
        palette=palette_by_ngenes,
        errorbar="sd",
        capsize=0.1,
        err_kws={"linewidth": 1},
        ax=ax,
    )
    
    ymin, ymax = ax.get_ylim()
    if ymax < 1:
        ax.set_ylim(0, 1)

    ax.set_title("")
    ax.set_xlabel("Number of taxa")

    ax.set_ylabel(metric_labels.get(metric, metric.replace("_", " ")))

    for label in ax.get_xticklabels():
        label.set_rotation(0)
        label.set_horizontalalignment("center")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linewidth=0.5, alpha=0.3)
    ax.set_axisbelow(True)

    # Keep only one shared legend
    if ax.legend_ is not None:
        ax.legend_.remove()

handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(
    handles,
    [str(x) for x in ngenes_order],
    title="Number of genes",
    loc="center left",
    frameon=False,
    bbox_to_anchor=(0.9, 0.5),
)

fig.savefig(output_pdf, bbox_inches="tight")
plt.close(fig)

print(f"Saved figure to: {output_pdf}")