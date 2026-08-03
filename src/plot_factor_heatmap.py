#!/usr/bin/env python

import argparse
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument("L_tsv", help="Factor loading matrix TSV: rows=factors, columns=genes")
parser.add_argument("outfile", help="Output figure path")
parser.add_argument("--top_n", type=int, default=5,
                    help="Top up/down genes per factor")
parser.add_argument(
    "--fit-summary",
    help="Fit summary TSV used to retain only factors with active Brownian variance",
)
parser.add_argument(
    "--active-variance-threshold",
    type=float,
    default=1e-5,
    help="Minimum Brownian variance for an active factor (default: 1e-5)",
)

args = parser.parse_args()

# Read loading matrix
L = pd.read_csv(args.L_tsv, sep="\t", index_col=0)

if args.fit_summary:
    summary = pd.read_csv(args.fit_summary, sep="\t", index_col=0)
    active_factors = []
    for parameter, value in summary["value"].items():
        match = re.fullmatch(r"sigma2_latent_LF([0-9]+)", parameter)
        if match and float(value) > args.active_variance_threshold:
            active_factors.append(f"factor_{int(match.group(1))}")
    missing = [factor for factor in active_factors if factor not in L.index]
    if missing:
        raise ValueError(
            f"Active factors from {args.fit_summary} are absent from {args.L_tsv}: "
            + ", ".join(missing)
        )
    if not active_factors:
        raise ValueError(
            "No factors have Brownian variance above "
            f"{args.active_variance_threshold:g}"
        )
    L = L.loc[active_factors]

factors = list(L.index)

# Order genes by factor, using top up then top down, skipping genes already assigned earlier
gene_order = []
seen_genes = set()

for factor in factors:

    vals = L.loc[factor]

    top_up = vals.nlargest(args.top_n).index.tolist()
    top_down = vals.nsmallest(args.top_n).index.tolist()

    for gene in top_up + top_down:
        if gene not in seen_genes:
            gene_order.append(gene)
            seen_genes.add(gene)

# Subset matrix
plot_mat = L[gene_order].copy()

# Z-score each gene across factors
plot_mat_z = plot_mat.copy()

for gene in plot_mat_z.columns:

    std = plot_mat_z[gene].std()

    if std > 0:
        plot_mat_z[gene] = (
            plot_mat_z[gene] - plot_mat_z[gene].mean()
        ) / std

# Plot settings
sns.set_context("talk")

plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "xtick.labelsize": 10,
    "ytick.labelsize": 13,
    "figure.titlesize": 18
})

# Plot heatmap
plt.figure(figsize=(18, 6))

ax = sns.heatmap(
    plot_mat_z,
    cmap="RdBu_r",
    center=0,
    xticklabels=True,
    yticklabels=True,
    cbar_kws={"label": "Gene loading z-score"}
)

ax.set_xlabel("Genes")
ax.set_ylabel("Factors")
ax.set_title(f"Top {args.top_n} Up/Down Genes Per Factor")

plt.setp(
    ax.get_xticklabels(),
    rotation=90,
    ha="center",
    va="top"
)

plt.setp(
    ax.get_yticklabels(),
    rotation=0
)

plt.tight_layout()

plt.savefig(
    args.outfile,
    dpi=300,
    bbox_inches="tight"
)

plt.close()
