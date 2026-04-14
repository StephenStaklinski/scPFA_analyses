#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.stats import spearmanr

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

perf_path = sys.argv[1]
stats_path = sys.argv[2]
outfile = sys.argv[3]

perf = pd.read_csv(perf_path, sep="\t")
stats = pd.read_csv(stats_path, sep="\t")

# Keep only moran/path results
perf = perf[perf["method"].isin(["moran", "path"])].copy()

# Compute metrics
perf["recall"] = perf["TP"] / (perf["TP"] + perf["FN"])
perf["precision"] = perf["TP"] / (perf["TP"] + perf["FP"])

# Collapse identical methods
perf = perf.groupby(["ntaxa", "sim_num"], as_index=False)[["recall", "precision"]].mean()

# Merge
merged = perf.merge(
    stats[["ntaxa", "sim_num", "pairwise_distance_min"]],
    on=["ntaxa", "sim_num"]
)

# Layout
ntaxa_values = sorted(merged["ntaxa"].unique())
ncols = len(ntaxa_values)

fig, axes = plt.subplots(
    nrows=2,
    ncols=ncols,
    figsize=(3 * ncols, 6),
    sharex=True,
    sharey=True
)

if ncols == 1:
    axes = [[axes[0]], [axes[1]]]

colors = plt.cm.viridis(
    [i / (len(ntaxa_values) - 1) if len(ntaxa_values) > 1 else 0.5
     for i in range(len(ntaxa_values))]
)

for col_idx, (ntaxa, color) in enumerate(zip(ntaxa_values, colors)):
    sub = merged[merged["ntaxa"] == ntaxa].sort_values("pairwise_distance_min")

    # ----- Recall -----
    ax = axes[0][col_idx]
    ax.scatter(sub["pairwise_distance_min"], sub["recall"],
               color=color, s=60, zorder=3, alpha=0.85)

    if len(sub) > 1:
        rho, _ = spearmanr(sub["pairwise_distance_min"], sub["recall"])
        label = f"ρ = {rho:.2f}"
    else:
        label = "ρ = NA"

    ax.text(0.98, 0.02, label, transform=ax.transAxes,
            ha="right", va="bottom", fontsize=9)

    ax.axhline(1.0, color="grey", linewidth=0.8, linestyle="--", zorder=1)
    ax.set_ylim(-0.05, 1.1)
    ax.set_xscale("log")
    ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.6)
    ax.set_ylabel("Recall")
    ax.set_title(f"{ntaxa} taxa")

    # ----- Precision -----
    ax = axes[1][col_idx]
    ax.scatter(sub["pairwise_distance_min"], sub["precision"],
               color=color, s=60, zorder=3, alpha=0.85)

    if len(sub) > 1:
        rho, _ = spearmanr(sub["pairwise_distance_min"], sub["precision"])
        label = f"ρ = {rho:.2f}"
    else:
        label = "ρ = NA"

    ax.text(0.98, 0.02, label, transform=ax.transAxes,
            ha="right", va="bottom", fontsize=9)

    ax.axhline(1.0, color="grey", linewidth=0.8, linestyle="--", zorder=1)
    ax.set_ylim(-0.05, 1.1)
    ax.set_xscale("log")
    ax.set_ylabel("Precision")
    ax.grid(True, which="both", linestyle=":", linewidth=0.5, alpha=0.6)

# Global labels
fig.supxlabel("Minimum pairwise tip distance")

plt.tight_layout()
plt.savefig(outfile)
