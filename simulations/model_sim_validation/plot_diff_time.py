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

diff_tsv = sys.argv[1]
time_tsv = sys.argv[2]
output_pdf = sys.argv[3]

metric_order = [
    "F_col_mean_corr",
    "L_row_mean_corr",
    "cell_cov_correlation",
    "gene_cov_correlation",
]
metric_labels = {
    "F_col_mean_corr": "F column mean corr",
    "L_row_mean_corr": "L row mean corr",
    "cell_cov_correlation": "Cell covariance correlation",
    "gene_cov_correlation": "Gene covariance correlation",
}

df = pd.read_csv(diff_tsv, sep="\t")
required_cols = {"NTAXA", "NGENES", "sim", *metric_order}
missing_cols = required_cols - set(df.columns)
if missing_cols:
    sys.stderr.write(f"Missing required columns: {sorted(missing_cols)}\n")
    sys.exit(1)

time_df = pd.read_csv(time_tsv, sep="\t")
required_time_cols = {"NTAXA", "NGENES", "sim_num", "sec"}
missing_time_cols = required_time_cols - set(time_df.columns)
if missing_time_cols:
    sys.stderr.write(f"Missing required time columns: {sorted(missing_time_cols)}\n")
    sys.exit(1)

for frame in (df, time_df):
    frame["NTAXA"] = pd.to_numeric(frame["NTAXA"])
    frame["NGENES"] = pd.to_numeric(frame["NGENES"])
for metric in metric_order:
    df[metric] = pd.to_numeric(df[metric], errors="coerce")
time_df["sec"] = pd.to_numeric(time_df["sec"], errors="coerce")

ntaxa_order = sorted(df["NTAXA"].dropna().unique())
ngenes_order = sorted(df["NGENES"].dropna().unique())

ngenes_colors = {
    500: "#E69F00",
    1000: "#56B4E9",
    5000: "#009E73",
}
fallback_colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", "#D55E00"]
palette_by_ngenes = {
    ngenes: ngenes_colors.get(ngenes, fallback_colors[i % len(fallback_colors)])
    for i, ngenes in enumerate(ngenes_order)
}

fig, axes = plt.subplots(1, 5, figsize=(22.5, 3.6))
plt.subplots_adjust(right=0.94, wspace=0.35)

for ax, metric in zip(axes[:4], metric_order):
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
    if ymin < 0 and ymax <= 1:
        ax.set_ylim(-1, 1)
    elif ymax < 1:
        ax.set_ylim(0, 1)

    ax.set_ylabel(metric_labels[metric])

runtime_ax = axes[4]
sns.barplot(
    data=time_df,
    x="NTAXA",
    y="sec",
    hue="NGENES",
    order=ntaxa_order,
    hue_order=ngenes_order,
    palette=palette_by_ngenes,
    errorbar="sd",
    capsize=0.1,
    err_kws={"linewidth": 1},
    ax=runtime_ax,
)
runtime_ax.set_ylabel("Runtime (s)")

for ax in axes:
    ax.set_title("")
    ax.set_xlabel("Number of taxa")
    for label in ax.get_xticklabels():
        label.set_rotation(0)
        label.set_horizontalalignment("center")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linewidth=0.5, alpha=0.3)
    ax.set_axisbelow(True)
    if ax.legend_ is not None:
        ax.legend_.remove()

handles, _ = axes[0].get_legend_handles_labels()
fig.legend(
    handles,
    [str(x) for x in ngenes_order],
    title="Number of genes",
    loc="center left",
    frameon=False,
    bbox_to_anchor=(0.945, 0.5),
)

fig.savefig(output_pdf, bbox_inches="tight")
plt.close(fig)

print(f"Saved figure to: {output_pdf}")
