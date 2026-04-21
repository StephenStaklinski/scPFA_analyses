#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import umap

sns.set_theme(style="ticks", context="paper")

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

sim_f_tsv = sys.argv[1]
sim_l_tsv = sys.argv[2]
sim_x_tsv = sys.argv[3]
fit_f_tsv = sys.argv[4]
fit_l_tsv = sys.argv[5]
fit_x_tsv = sys.argv[6]
output_file = sys.argv[7]

# Read matrices
sim_F_df = pd.read_csv(sim_f_tsv, sep="\t", index_col=0)
sim_L_df = pd.read_csv(sim_l_tsv, sep="\t", index_col=0)
sim_X_df = pd.read_csv(sim_x_tsv, sep="\t", index_col=0)
fit_F_df = pd.read_csv(fit_f_tsv, sep="\t", index_col=0)
fit_L_df = pd.read_csv(fit_l_tsv, sep="\t", index_col=0)
fit_X_df = pd.read_csv(fit_x_tsv, sep="\t", index_col=0)

# Calculate Z = F * L
sim_Z_df = sim_F_df @ sim_L_df
fit_Z_df = fit_F_df @ fit_L_df

sim_F = sim_F_df.to_numpy()
fit_F = fit_F_df.to_numpy()

sim_cell_colors = np.arange(sim_F.shape[0])
fit_cell_colors = np.arange(fit_F.shape[0])

# PCA
sim_pca = PCA(n_components=2).fit_transform(sim_F)
fit_pca = PCA(n_components=2).fit_transform(fit_F)

# UMAP
sim_umap = umap.UMAP(n_components=2, random_state=42).fit_transform(sim_F)
fit_umap = umap.UMAP(n_components=2, random_state=42).fit_transform(fit_F)

fig, axes = plt.subplots(2, 6, figsize=(24, 8))

# -------- top row: sim --------
sns.heatmap(
    sim_F_df,
    ax=axes[0, 0],
    cmap="coolwarm",
    linewidths=0,
    cbar=True,
    xticklabels=False,
    yticklabels=False
)
axes[0, 0].set_title("Sim F")
axes[0, 0].set_xlabel("Factors")
axes[0, 0].set_ylabel("Cells")

sns.heatmap(
    sim_L_df,
    ax=axes[0, 1],
    cmap="coolwarm",
    linewidths=0,
    cbar=True,
    xticklabels=False,
    yticklabels=False
)
axes[0, 1].set_title("Sim L")
axes[0, 1].set_xlabel("Genes")
axes[0, 1].set_ylabel("Factors")

sns.heatmap(
    sim_Z_df,
    ax=axes[0, 2],
    cmap="coolwarm",
    linewidths=0,
    cbar=True,
    xticklabels=False,
    yticklabels=False
)
axes[0, 2].set_title("Sim Z")
axes[0, 2].set_xlabel("Cells")
axes[0, 2].set_ylabel("Genes")

sns.heatmap(
    sim_X_df,
    ax=axes[0, 3],
    cmap="coolwarm",
    linewidths=0,
    cbar=True,
    xticklabels=False,
    yticklabels=False
)
axes[0, 3].set_title("Sim X")
axes[0, 3].set_xlabel("Cells")
axes[0, 3].set_ylabel("Genes")

axes[0, 4].scatter(sim_pca[:, 0], sim_pca[:, 1], s=8, c=sim_cell_colors, cmap="viridis")
axes[0, 4].set_title("Sim F PCA")
axes[0, 4].set_xlabel("PC1")
axes[0, 4].set_ylabel("PC2")

axes[0, 5].scatter(sim_umap[:, 0], sim_umap[:, 1], s=8, c=sim_cell_colors, cmap="viridis")
axes[0, 5].set_title("Sim F UMAP")
axes[0, 5].set_xlabel("UMAP1")
axes[0, 5].set_ylabel("UMAP2")

# -------- bottom row: fit --------
sns.heatmap(
    fit_F_df,
    ax=axes[1, 0],
    cmap="coolwarm",
    linewidths=0,
    cbar=True,
    xticklabels=False,
    yticklabels=False
)
axes[1, 0].set_title("Fit F")
axes[1, 0].set_xlabel("Factors")
axes[1, 0].set_ylabel("Cells")

sns.heatmap(
    fit_L_df,
    ax=axes[1, 1],
    cmap="coolwarm",
    linewidths=0,
    cbar=True,
    xticklabels=False,
    yticklabels=False
)
axes[1, 1].set_title("Fit L")
axes[1, 1].set_xlabel("Genes")
axes[1, 1].set_ylabel("Factors")

sns.heatmap(
    fit_Z_df,
    ax=axes[1, 2],
    cmap="coolwarm",
    linewidths=0,
    cbar=True,
    xticklabels=False,
    yticklabels=False
)
axes[1, 2].set_title("Fit Z")
axes[1, 2].set_xlabel("Cells")
axes[1, 2].set_ylabel("Genes")

sns.heatmap(
    fit_X_df,
    ax=axes[1, 3],
    cmap="coolwarm",
    linewidths=0,
    cbar=True,
    xticklabels=False,
    yticklabels=False
)
axes[1, 3].set_title("Fit X")
axes[1, 3].set_xlabel("Cells")
axes[1, 3].set_ylabel("Genes")

axes[1, 4].scatter(fit_pca[:, 0], fit_pca[:, 1], s=8, c=fit_cell_colors, cmap="viridis")
axes[1, 4].set_title("Fit F PCA")
axes[1, 4].set_xlabel("PC1")
axes[1, 4].set_ylabel("PC2")

axes[1, 5].scatter(fit_umap[:, 0], fit_umap[:, 1], s=8, c=fit_cell_colors, cmap="viridis")
axes[1, 5].set_title("Fit F UMAP")
axes[1, 5].set_xlabel("UMAP1")
axes[1, 5].set_ylabel("UMAP2")

plt.tight_layout()
plt.savefig(output_file)
plt.close()