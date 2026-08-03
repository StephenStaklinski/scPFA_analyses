#!/usr/bin/env python3

import os
import re

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
if len(sys.argv) == 5:
    output_file = sys.argv[4]
    if not fit_f_tsv.endswith(".F.tsv"):
        raise ValueError(
            "Cannot infer the fit summary path because the F matrix path does "
            f"not end in '.F.tsv': {fit_f_tsv}"
        )
    fit_summary_tsv = fit_f_tsv.removesuffix(".F.tsv") + ".summary.tsv"
elif len(sys.argv) == 6:
    fit_summary_tsv = sys.argv[4]
    output_file = sys.argv[5]
else:
    raise SystemExit(
        "Usage: plot_fit_panel.py F.tsv L.tsv X.tsv [summary.tsv] output"
    )

# Read matrices
fit_F_df = pd.read_csv(fit_f_tsv, sep="\t", index_col=0)
fit_L_df = pd.read_csv(fit_l_tsv, sep="\t", index_col=0)
fit_X_df = pd.read_csv(fit_x_tsv, sep="\t", index_col=0)
fit_summary_df = pd.read_csv(fit_summary_tsv, sep="\t", index_col=0)

brownian_variance = []
for parameter, value in fit_summary_df["value"].items():
    match = re.fullmatch(r"sigma2_latent_LF([0-9]+)", parameter)
    if match:
        brownian_variance.append((int(match.group(1)), float(value)))
brownian_variance.sort()

if len(brownian_variance) != fit_F_df.shape[1]:
    raise ValueError(
        f"Expected {fit_F_df.shape[1]} sigma2_latent values in "
        f"{fit_summary_tsv}, found {len(brownian_variance)}"
    )

# Calculate Z = F * L
fit_Z_df = fit_F_df @ fit_L_df

fit_F = fit_F_df.to_numpy()
fit_Z = fit_Z_df.to_numpy()
fit_X = fit_X_df.to_numpy()

fit_cell_colors = np.arange(fit_F.shape[0])

# Standardize each factor across cells for the F heatmap.
fit_F_std = fit_F_df.std(axis=0, ddof=0).replace(0, 1)
fit_F_z_df = (fit_F_df - fit_F_df.mean(axis=0)) / fit_F_std


def zero_centered_heatmap(matrix, ax):
    max_abs = np.nanmax(np.abs(matrix.to_numpy()))
    if not np.isfinite(max_abs) or max_abs == 0:
        max_abs = 1.0
    return sns.heatmap(
        matrix,
        ax=ax,
        cmap="coolwarm",
        vmin=-max_abs,
        vmax=max_abs,
        center=0,
        linewidths=0,
        cbar=True,
        xticklabels=False,
        yticklabels=False,
    )


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

fig, axes = plt.subplots(1, 11, figsize=(44, 4))

# -------- fit row --------
factor_labels = list(fit_F_df.columns)
axes[0].bar(
    factor_labels,
    [value for _, value in brownian_variance],
    color=sns.color_palette("deep")[0],
)
axes[0].set_title("Brownian Variance")
axes[0].set_xlabel("Factors")
axes[0].set_ylabel(r"Brownian variance ($\sigma^2$)")
axes[0].tick_params(axis="x", rotation=90)

zero_centered_heatmap(fit_F_z_df, axes[1])
axes[1].set_title("Fit F (z-score)")
axes[1].set_xlabel("Factors")
axes[1].set_ylabel("Cells")

zero_centered_heatmap(fit_L_df, axes[2])
axes[2].set_title("Fit L")
axes[2].set_xlabel("Genes")
axes[2].set_ylabel("Factors")

zero_centered_heatmap(fit_Z_df, axes[3])
axes[3].set_title("Fit Z")
axes[3].set_xlabel("Cells")
axes[3].set_ylabel("Genes")

zero_centered_heatmap(fit_X_df, axes[4])
axes[4].set_title("Fit X")
axes[4].set_xlabel("Cells")
axes[4].set_ylabel("Genes")

axes[5].scatter(
    fit_F_pca[:, 0],
    fit_F_pca[:, 1],
    s=8,
    c=fit_cell_colors,
    cmap="viridis"
)
axes[5].set_title("Fit F PCA")
axes[5].set_xlabel("PC1")
axes[5].set_ylabel("PC2")

axes[6].scatter(
    fit_F_umap[:, 0],
    fit_F_umap[:, 1],
    s=8,
    c=fit_cell_colors,
    cmap="viridis"
)
axes[6].set_title("Fit F UMAP")
axes[6].set_xlabel("UMAP1")
axes[6].set_ylabel("UMAP2")

axes[7].scatter(
    fit_Z_pca[:, 0],
    fit_Z_pca[:, 1],
    s=8,
    c=fit_cell_colors,
    cmap="viridis"
)
axes[7].set_title("Fit Z PCA")
axes[7].set_xlabel("PC1")
axes[7].set_ylabel("PC2")

axes[8].scatter(
    fit_Z_umap[:, 0],
    fit_Z_umap[:, 1],
    s=8,
    c=fit_cell_colors,
    cmap="viridis"
)
axes[8].set_title("Fit Z UMAP")
axes[8].set_xlabel("UMAP1")
axes[8].set_ylabel("UMAP2")

axes[9].scatter(
    fit_X_pca[:, 0],
    fit_X_pca[:, 1],
    s=8,
    c=fit_cell_colors,
    cmap="viridis"
)
axes[9].set_title("Fit X PCA")
axes[9].set_xlabel("PC1")
axes[9].set_ylabel("PC2")

axes[10].scatter(
    fit_X_umap[:, 0],
    fit_X_umap[:, 1],
    s=8,
    c=fit_cell_colors,
    cmap="viridis"
)
axes[10].set_title("Fit X UMAP")
axes[10].set_xlabel("UMAP1")
axes[10].set_ylabel("UMAP2")

plt.tight_layout()
plt.savefig(output_file)
plt.close()
