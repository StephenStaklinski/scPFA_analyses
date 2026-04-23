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

required_cols = {
    "sim_num",
    "n_cells_simulated",
    "n_genes_simulated",
    "brownian_negll_fit",
    "brownian_negll_simulated",
    "observation_negll_fit",
    "observation_negll_simulated",
}
missing_cols = required_cols - set(df.columns)
if missing_cols:
    sys.stderr.write(f"Missing required columns: {sorted(missing_cols)}\n")
    sys.exit(1)

numeric_cols = [
    "sim_num",
    "n_cells_simulated",
    "n_genes_simulated",
    "brownian_negll_fit",
    "brownian_negll_simulated",
    "observation_negll_fit",
    "observation_negll_simulated",
]

for col in numeric_cols:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df = df.dropna(subset=numeric_cols).copy()

# Convert negative log-likelihood (negll) to log-likelihood (ll)
df["brownian_ll_fit"] = -df["brownian_negll_fit"]
df["brownian_ll_simulated"] = -df["brownian_negll_simulated"]

df["observation_ll_fit"] = -df["observation_negll_fit"]
df["observation_ll_simulated"] = -df["observation_negll_simulated"]

# Relative percent change in log-likelihood:
# 100 * (fit - simulated) / |simulated|
df["brownian_pct_change"] = (
    100.0
    * (df["brownian_ll_fit"] - df["brownian_ll_simulated"])
    / df["brownian_ll_simulated"].abs()
)

df["observation_pct_change"] = (
    100.0
    * (df["observation_ll_fit"] - df["observation_ll_simulated"])
    / df["observation_ll_simulated"].abs()
)

cell_order = sorted(df["n_cells_simulated"].unique())
gene_order = sorted(df["n_genes_simulated"].unique())

gene_colors = {
    500: "#E69F00",
    1000: "#56B4E9",
    5000: "#009E73",
}
fallback_colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", "#D55E00"]

palette_by_gene = {}
for i, ngenes in enumerate(gene_order):
    palette_by_gene[ngenes] = gene_colors.get(ngenes, fallback_colors[i % len(fallback_colors)])

fig, axes = plt.subplots(1, 2, figsize=(4.5 * 2, 3.6 * 1))

plot_specs = [
    ("brownian_pct_change", "Brownian log-likelihood\nrelative difference (%)"),
    ("observation_pct_change", "Observation log-likelihood\nrelative difference (%)")
]

for i, (metric, ylabel) in enumerate(plot_specs):
    ax = axes[i]

    sns.barplot(
        data=df,
        x="n_cells_simulated",
        y=metric,
        hue="n_genes_simulated",
        order=cell_order,
        hue_order=gene_order,
        palette=palette_by_gene,
        errorbar="sd",
        capsize=0.1,
        err_kws={"linewidth": 1},
        ax=ax,
    )

    ymin, ymax = ax.get_ylim()
    if ymin > -1 and ymax < 1:
        ax.set_ylim(-1, 1)
    ax.axhline(0, color="black", linewidth=0.8)

    ax.set_title("")
    ax.set_xlabel("Number of taxa")
    ax.set_ylabel(ylabel)

    for label in ax.get_xticklabels():
        label.set_rotation(0)
        label.set_horizontalalignment("center")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linewidth=0.5, alpha=0.3)
    ax.set_axisbelow(True)

    if ax.legend_ is not None:
        ax.legend_.remove()

handles, labels = axes[0].get_legend_handles_labels()
fig.legend(
    handles,
    [str(x) for x in gene_order],
    title="Number of genes",
    loc="center left",
    frameon=False,
    bbox_to_anchor=(0.9, 0.5),
)

plt.subplots_adjust(wspace=0.4)
fig.savefig(output_pdf, bbox_inches="tight")
plt.close(fig)

print(f"Saved figure to: {output_pdf}")