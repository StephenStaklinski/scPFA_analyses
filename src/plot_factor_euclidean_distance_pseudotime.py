#!/usr/bin/env python

import argparse
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns


parser = argparse.ArgumentParser()
parser.add_argument("pca_tsv")
parser.add_argument("umap_tsv")
parser.add_argument("distance_tsv")
parser.add_argument("out_prefix")
args = parser.parse_args()


sns.set_theme(
    context="paper",
    style="white",
    font_scale=1.2,
    rc={
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    },
)

mpl.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


dist = pd.read_csv(args.distance_tsv, sep="\t")
pca = pd.read_csv(args.pca_tsv, sep="\t")
umap = pd.read_csv(args.umap_tsv, sep="\t")

if "cell" not in dist.columns:
    raise ValueError("distance_tsv must contain a 'cell' column.")
if "total_euclidean_distance" not in dist.columns:
    raise ValueError("distance_tsv must contain a 'total_euclidean_distance' column.")

for name, df, cols in [
    ("PCA", pca, ["cell", "PC1", "PC2"]),
    ("UMAP", umap, ["cell", "UMAP1", "UMAP2"]),
]:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{name} file is missing columns: {missing}")

dist = dist[["cell", "total_euclidean_distance"]].copy()

pca_df = pca[["cell", "PC1", "PC2"]].merge(dist, on="cell", how="inner")
umap_df = umap[["cell", "UMAP1", "UMAP2"]].merge(dist, on="cell", how="inner")

if len(pca_df) == 0:
    raise ValueError("No matching cells found between distance file and PCA file.")
if len(umap_df) == 0:
    raise ValueError("No matching cells found between distance file and UMAP file.")

print(f"PCA: plotted {len(pca_df)} matching cells.")
print(f"UMAP: plotted {len(umap_df)} matching cells.")

vmin = dist["total_euclidean_distance"].min()
vmax = dist["total_euclidean_distance"].max()

fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.2), constrained_layout=True)

sc1 = axes[0].scatter(
    pca_df["PC1"],
    pca_df["PC2"],
    c=pca_df["total_euclidean_distance"],
    cmap="Spectral_r",
    vmin=vmin,
    vmax=vmax,
    s=90,
    alpha=0.95,
    linewidths=0,
    rasterized=True,
)

axes[0].set_title("PCA")
axes[0].set_xlabel("PC1")
axes[0].set_ylabel("PC2")
sns.despine(ax=axes[0])

sc2 = axes[1].scatter(
    umap_df["UMAP1"],
    umap_df["UMAP2"],
    c=umap_df["total_euclidean_distance"],
    cmap="Spectral_r",
    vmin=vmin,
    vmax=vmax,
    s=90,
    alpha=0.95,
    linewidths=0,
    rasterized=True,
)

axes[1].set_title("UMAP")
axes[1].set_xlabel("UMAP1")
axes[1].set_ylabel("UMAP2")
sns.despine(ax=axes[1])

cbar = fig.colorbar(sc2, ax=axes, fraction=0.046, pad=0.03)
cbar.set_label("Euclidean factor pseudotime")
cbar.outline.set_visible(False)

fig.savefig(f"{args.out_prefix}.euclidean_factor_pseudotime.pdf")
plt.close(fig)