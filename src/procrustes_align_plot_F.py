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
fit_f_tsv = sys.argv[2]
output_file = sys.argv[3]

sim_F_df = pd.read_csv(sim_f_tsv, sep="\t", index_col=0)
fit_F_df = pd.read_csv(fit_f_tsv, sep="\t", index_col=0)

common_rows = sim_F_df.index.intersection(fit_F_df.index)
if len(common_rows) == 0:
    raise ValueError("No overlapping row names between sim F and fit F.")

sim_F_df = sim_F_df.loc[common_rows]
fit_F_df = fit_F_df.loc[common_rows]

if sim_F_df.shape[1] != fit_F_df.shape[1]:
    raise ValueError(
        f"Sim and fit F must have same number of columns after row matching, got {sim_F_df.shape[1]} and {fit_F_df.shape[1]}"
    )

sim_F = sim_F_df.to_numpy(dtype=float)
fit_F = fit_F_df.to_numpy(dtype=float)

sim_F_centered = sim_F - sim_F.mean(axis=0, keepdims=True)
fit_F_centered = fit_F - fit_F.mean(axis=0, keepdims=True)

# Orthogonal Procrustes
M = fit_F_centered.T @ sim_F_centered
U, s, Vt = np.linalg.svd(M, full_matrices=False)
R = U @ Vt
fit_F_aligned = fit_F_centered @ R

fit_F_aligned_df = pd.DataFrame(
    fit_F_aligned,
    index=common_rows,
    columns=sim_F_df.columns,
)

sim_colors = np.arange(sim_F.shape[0])
fit_colors = np.arange(fit_F.shape[0])

# PCA
sim_pca = PCA(n_components=2).fit_transform(sim_F)
fit_pca = PCA(n_components=2).fit_transform(fit_F)
fit_aligned_pca = PCA(n_components=2).fit_transform(fit_F_aligned)

# UMAP
sim_umap = umap.UMAP(n_components=2, random_state=42).fit_transform(sim_F)
fit_umap = umap.UMAP(n_components=2, random_state=42).fit_transform(fit_F)
fit_aligned_umap = umap.UMAP(n_components=2, random_state=42).fit_transform(fit_F_aligned)

fig, axes = plt.subplots(3, 3, figsize=(12, 12))

# ---------------- top row: sim ----------------
sns.heatmap(
    sim_F_df,
    ax=axes[0, 0],
    cmap="coolwarm",
    cbar=True,
    xticklabels=False,
    yticklabels=False,
)
axes[0, 0].set_title("Sim F")
axes[0, 0].set_xlabel("Factors")
axes[0, 0].set_ylabel("Cells")

axes[0, 1].scatter(sim_pca[:, 0], sim_pca[:, 1], s=8, c=sim_colors, cmap="viridis")
axes[0, 1].set_title("Sim F PCA")
axes[0, 1].set_xlabel("PC1")
axes[0, 1].set_ylabel("PC2")

axes[0, 2].scatter(sim_umap[:, 0], sim_umap[:, 1], s=8, c=sim_colors, cmap="viridis")
axes[0, 2].set_title("Sim F UMAP")
axes[0, 2].set_xlabel("UMAP1")
axes[0, 2].set_ylabel("UMAP2")

# ---------------- middle row: original fit ----------------
sns.heatmap(
    fit_F_df,
    ax=axes[1, 0],
    cmap="coolwarm",
    cbar=True,
    xticklabels=False,
    yticklabels=False,
)
axes[1, 0].set_title("Fit F")
axes[1, 0].set_xlabel("Factors")
axes[1, 0].set_ylabel("Cells")

axes[1, 1].scatter(fit_pca[:, 0], fit_pca[:, 1], s=8, c=fit_colors, cmap="viridis")
axes[1, 1].set_title("Fit F PCA")
axes[1, 1].set_xlabel("PC1")
axes[1, 1].set_ylabel("PC2")

axes[1, 2].scatter(fit_umap[:, 0], fit_umap[:, 1], s=8, c=fit_colors, cmap="viridis")
axes[1, 2].set_title("Fit F UMAP")
axes[1, 2].set_xlabel("UMAP1")
axes[1, 2].set_ylabel("UMAP2")

# ---------------- bottom row: aligned fit ----------------
sns.heatmap(
    fit_F_aligned_df,
    ax=axes[2, 0],
    cmap="coolwarm",
    cbar=True,
    xticklabels=False,
    yticklabels=False,
)
axes[2, 0].set_title("Fit F aligned")
axes[2, 0].set_xlabel("Factors")
axes[2, 0].set_ylabel("Cells")

axes[2, 1].scatter(
    fit_aligned_pca[:, 0],
    fit_aligned_pca[:, 1],
    s=8,
    c=fit_colors,
    cmap="viridis",
)
axes[2, 1].set_title("Fit F aligned PCA")
axes[2, 1].set_xlabel("PC1")
axes[2, 1].set_ylabel("PC2")

axes[2, 2].scatter(
    fit_aligned_umap[:, 0],
    fit_aligned_umap[:, 1],
    s=8,
    c=fit_colors,
    cmap="viridis",
)
axes[2, 2].set_title("Fit F aligned UMAP")
axes[2, 2].set_xlabel("UMAP1")
axes[2, 2].set_ylabel("UMAP2")

obj_before = np.linalg.norm(fit_F_centered - sim_F_centered, ord="fro") ** 2
obj_after = np.linalg.norm(fit_F_aligned - sim_F_centered, ord="fro") ** 2
mse_before = np.mean((fit_F_centered - sim_F_centered) ** 2)
mse_after = np.mean((fit_F_aligned - sim_F_centered) ** 2)

fig.suptitle(
    f"Orthogonal Procrustes on F\n"
    f"Objective before = {obj_before:.4g}, after = {obj_after:.4g}; "
    f"MSE before = {mse_before:.4g}, after = {mse_after:.4g}",
    y=0.995,
)

plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(output_file, bbox_inches="tight")
plt.close()