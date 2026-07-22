#!/usr/bin/env python3

import os

# Prevent script from consuming all CPU cores
N_THREADS = "1"
os.environ["OMP_NUM_THREADS"] = N_THREADS
os.environ["OPENBLAS_NUM_THREADS"] = N_THREADS
os.environ["MKL_NUM_THREADS"] = N_THREADS
os.environ["VECLIB_MAXIMUM_THREADS"] = N_THREADS
os.environ["NUMEXPR_NUM_THREADS"] = N_THREADS

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import umap

# Set plotting defaults
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

fit_f_tsv = sys.argv[1]
fit_l_tsv = sys.argv[2]
fit_x_tsv = sys.argv[3]
output_file = sys.argv[4]

# Read matrices
fit_F_df = pd.read_csv(fit_f_tsv, sep="\t", index_col=0)
fit_L_df = pd.read_csv(fit_l_tsv, sep="\t", index_col=0)
fit_X_df = pd.read_csv(fit_x_tsv, sep="\t", index_col=0)

# Calculate Z = F * L
fit_Z_df = fit_F_df @ fit_L_df

fit_F = fit_F_df.to_numpy()
fit_Z = fit_Z_df.to_numpy()
fit_X = fit_X_df.to_numpy()

fit_cell_colors = np.arange(fit_F.shape[0])


def two_dim_pca(matrix):
    n_components = min(2, matrix.shape[0], matrix.shape[1])
    if n_components <= 0:
        return np.zeros((matrix.shape[0], 2))
    coords = PCA(n_components=n_components).fit_transform(matrix)
    if coords.shape[1] == 1:
        coords = np.column_stack([coords[:, 0], np.zeros(coords.shape[0])])
    return coords


def two_dim_umap(matrix):
    if matrix.shape[1] < 2:
        centered = matrix[:, 0] - np.mean(matrix[:, 0])
        return np.column_stack([centered, np.zeros(matrix.shape[0])])
    return umap.UMAP(n_components=2, random_state=42).fit_transform(matrix)


# PCA
fit_F_pca = two_dim_pca(fit_F)
fit_Z_pca = two_dim_pca(fit_Z)
fit_X_pca = two_dim_pca(fit_X)

# UMAP
fit_F_umap = two_dim_umap(fit_F)
fit_Z_umap = two_dim_umap(fit_Z)
fit_X_umap = two_dim_umap(fit_X)

fig, axes = plt.subplots(1, 10, figsize=(40, 4))

# -------- fit row --------
sns.heatmap(
    fit_F_df,
    ax=axes[0],
    cmap="coolwarm",
    linewidths=0,
    cbar=True,
    xticklabels=False,
    yticklabels=False
)
axes[0].set_title("Fit F")
axes[0].set_xlabel("Factors")
axes[0].set_ylabel("Cells")

sns.heatmap(
    fit_L_df,
    ax=axes[1],
    cmap="coolwarm",
    linewidths=0,
    cbar=True,
    xticklabels=False,
    yticklabels=False
)
axes[1].set_title("Fit L")
axes[1].set_xlabel("Genes")
axes[1].set_ylabel("Factors")

sns.heatmap(
    fit_Z_df,
    ax=axes[2],
    cmap="coolwarm",
    linewidths=0,
    cbar=True,
    xticklabels=False,
    yticklabels=False
)
axes[2].set_title("Fit Z")
axes[2].set_xlabel("Cells")
axes[2].set_ylabel("Genes")

sns.heatmap(
    fit_X_df,
    ax=axes[3],
    cmap="coolwarm",
    linewidths=0,
    cbar=True,
    xticklabels=False,
    yticklabels=False
)
axes[3].set_title("Fit X")
axes[3].set_xlabel("Cells")
axes[3].set_ylabel("Genes")

axes[4].scatter(
    fit_F_pca[:, 0],
    fit_F_pca[:, 1],
    s=8,
    c=fit_cell_colors,
    cmap="viridis"
)
axes[4].set_title("Fit F PCA")
axes[4].set_xlabel("PC1")
axes[4].set_ylabel("PC2")

axes[5].scatter(
    fit_F_umap[:, 0],
    fit_F_umap[:, 1],
    s=8,
    c=fit_cell_colors,
    cmap="viridis"
)
axes[5].set_title("Fit F UMAP")
axes[5].set_xlabel("UMAP1")
axes[5].set_ylabel("UMAP2")

axes[6].scatter(
    fit_Z_pca[:, 0],
    fit_Z_pca[:, 1],
    s=8,
    c=fit_cell_colors,
    cmap="viridis"
)
axes[6].set_title("Fit Z PCA")
axes[6].set_xlabel("PC1")
axes[6].set_ylabel("PC2")

axes[7].scatter(
    fit_Z_umap[:, 0],
    fit_Z_umap[:, 1],
    s=8,
    c=fit_cell_colors,
    cmap="viridis"
)
axes[7].set_title("Fit Z UMAP")
axes[7].set_xlabel("UMAP1")
axes[7].set_ylabel("UMAP2")

axes[8].scatter(
    fit_X_pca[:, 0],
    fit_X_pca[:, 1],
    s=8,
    c=fit_cell_colors,
    cmap="viridis"
)
axes[8].set_title("Fit X PCA")
axes[8].set_xlabel("PC1")
axes[8].set_ylabel("PC2")

axes[9].scatter(
    fit_X_umap[:, 0],
    fit_X_umap[:, 1],
    s=8,
    c=fit_cell_colors,
    cmap="viridis"
)
axes[9].set_title("Fit X UMAP")
axes[9].set_xlabel("UMAP1")
axes[9].set_ylabel("UMAP2")

plt.tight_layout()
plt.savefig(output_file)
plt.close()
