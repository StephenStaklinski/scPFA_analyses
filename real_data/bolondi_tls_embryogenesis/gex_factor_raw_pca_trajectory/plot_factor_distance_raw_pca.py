#!/usr/bin/env python

import argparse
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import seaborn as sns


parser = argparse.ArgumentParser(
    description=(
        "Plot fitted active-factor distance from the root on existing raw-data "
        "PCA coordinates."
    )
)
parser.add_argument("raw_pca_tsv")
parser.add_argument("latent_flow_tsv")
parser.add_argument("fit_summary_tsv")
parser.add_argument("output_pdf")
parser.add_argument("--active-variance-threshold", type=float, default=1e-5)
args = parser.parse_args()


sns.set_theme(
    context="paper",
    style="white",
    font_scale=1.15,
    rc={
        "axes.edgecolor": "0.15",
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
    },
)
mpl.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


raw_pca = pd.read_csv(args.raw_pca_tsv, sep="\t")
latent = pd.read_csv(args.latent_flow_tsv, sep="\t")
summary = pd.read_csv(args.fit_summary_tsv, sep="\t", index_col=0)

required_pca = {"cell", "PC1", "PC2"}
missing_pca = required_pca - set(raw_pca.columns)
if missing_pca:
    raise ValueError("Raw PCA table is missing: " + ", ".join(sorted(missing_pca)))

available_factors = [column for column in latent.columns if column.startswith("factor_")]
active_factors = []
for parameter, value in summary["value"].items():
    match = re.fullmatch(r"sigma2_latent_LF([0-9]+)", parameter)
    if match and float(value) > args.active_variance_threshold:
        active_factors.append((int(match.group(1)), f"factor_{int(match.group(1))}"))
active_factors = [factor for _, factor in sorted(active_factors)]

missing_factors = [factor for factor in active_factors if factor not in available_factors]
if missing_factors:
    raise ValueError(
        "Active factors absent from latent-flow table: " + ", ".join(missing_factors)
    )
if not active_factors:
    raise ValueError("No active factors found in fit summary.")

root = latent.loc[latent["parent_id"] < 0]
if len(root) != 1:
    raise ValueError("Expected exactly one reconstructed root.")
root_values = root.iloc[0][active_factors].to_numpy(float)

tips = latent.loc[latent["is_tip"] == 1, ["node_name"] + active_factors].copy()
tip_values = tips[active_factors].to_numpy(float)
tips["latent_distance_from_root"] = np.linalg.norm(
    tip_values - root_values[None, :], axis=1
)
tips = tips.rename(columns={"node_name": "cell"})

plot_df = raw_pca[["cell", "PC1", "PC2"]].merge(
    tips[["cell", "latent_distance_from_root"]],
    on="cell",
    how="inner",
    validate="one_to_one",
)
if plot_df.empty:
    raise ValueError("No cell barcodes match between raw PCA and fitted tree tips.")
if len(plot_df) != len(tips):
    raise ValueError(
        f"Matched {len(plot_df)} of {len(tips)} fitted tips to raw PCA coordinates."
    )

fig = plt.figure(figsize=(5.6, 4.7))
ax = fig.add_axes([0.12, 0.14, 0.68, 0.78])
scatter = ax.scatter(
    plot_df["PC1"],
    plot_df["PC2"],
    c=plot_df["latent_distance_from_root"],
    cmap=sns.color_palette("YlOrBr", as_cmap=True),
    vmin=0.0,
    vmax=float(plot_df["latent_distance_from_root"].max()),
    s=27,
    alpha=0.98,
    linewidths=0.25,
    edgecolors="white",
    rasterized=True,
)

ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
sns.despine(ax=ax)

cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.035)
cbar.set_label("Euclidean distance from root")
cbar.ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("{x:.2f}"))
cbar.update_ticks()
cbar.outline.set_visible(False)

fig.savefig(args.output_pdf)
plt.close(fig)
