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

# Colorblind-safe palette
ngenes_colors = {
    500: "#E69F00",
    1000: "#56B4E9",
    5000: "#009E73",
}

df = pd.read_csv(input_tsv, sep="\t")

required_cols = {"NTAXA", "NGENES", "sim_num", "sec"}
missing_cols = required_cols - set(df.columns)
if missing_cols:
    sys.stderr.write(f"Missing required columns: {sorted(missing_cols)}\n")
    sys.exit(1)

df["NTAXA"] = pd.to_numeric(df["NTAXA"])
df["NGENES"] = pd.to_numeric(df["NGENES"])
df["sec"] = pd.to_numeric(df["sec"])

ntaxa_order = sorted(df["NTAXA"].unique())
ngenes_order = sorted(df["NGENES"].unique())

palette_by_ngenes = {}
fallback_colors = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", "#D55E00"]
for i, ngenes in enumerate(ngenes_order):
    palette_by_ngenes[ngenes] = ngenes_colors.get(ngenes, fallback_colors[i % len(fallback_colors)])

fig, ax = plt.subplots(1, 1, figsize=(4.75, 3.5))
plt.subplots_adjust(right=0.80)

sns.barplot(
    data=df,
    x="NTAXA",
    y="sec",
    hue="NGENES",
    order=ntaxa_order,
    hue_order=ngenes_order,
    palette=palette_by_ngenes,
    errorbar="sd",
    capsize=0.1,
    err_kws={"linewidth": 1},
    ax=ax,
)

ax.set_xlabel("Number of taxa")
for label in ax.get_xticklabels():
    label.set_rotation(0)
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
    [str(x) for x in ngenes_order],
    title="Number of genes",
    loc="center left",
    frameon=False,
    bbox_to_anchor=(0.82, 0.5),
)

fig.savefig(output_pdf, bbox_inches="tight")
plt.close(fig)

print(f"Saved figure to: {output_pdf}")