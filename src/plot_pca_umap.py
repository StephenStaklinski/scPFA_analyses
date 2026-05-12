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


def preprocess_counts(X, target_sum=1e4, n_variable_genes=2000):
    # Row/library-size normalize each cell
    library_sizes = X.sum(axis=1, keepdims=True)

    # Avoid divide-by-zero for empty rows
    library_sizes[library_sizes == 0] = 1.0

    X_norm = X / library_sizes * target_sum

    # Standard log transform for count data
    X_log = np.log1p(X_norm)

    # Keep top variable genes
    n_keep = min(n_variable_genes, X_log.shape[1])
    gene_vars = np.var(X_log, axis=0, ddof=1)
    top_gene_idx = np.argsort(gene_vars)[-n_keep:]

    X_hvg = X_log[:, top_gene_idx]

    return X_hvg


def scale_genes(X):
    # Standard single-cell preprocessing usually scales genes before PCA
    gene_means = X.mean(axis=0, keepdims=True)
    gene_stds = X.std(axis=0, keepdims=True)

    gene_stds[gene_stds == 0] = 1.0

    return (X - gene_means) / gene_stds


x_tsv = sys.argv[1]
output_file = sys.argv[2]

preprocess = True
scale_before_pca = True
n_pcs_for_umap = 30

# Read X matrix
X_df = pd.read_csv(x_tsv, sep="\t", index_col=0)
X = X_df.to_numpy(dtype=float)

if preprocess:
    X = preprocess_counts(X)

if scale_before_pca:
    X = scale_genes(X)

cell_colors = np.arange(X.shape[0])

# PCA for plotting
X_pca = PCA(n_components=2).fit_transform(X)

# PCA denoising before UMAP
n_pcs = min(n_pcs_for_umap, X.shape[0] - 1, X.shape[1])
X_pca_for_umap = PCA(n_components=n_pcs).fit_transform(X)

# UMAP on PCs, not raw genes
X_umap = umap.UMAP(
    n_components=2,
    n_neighbors=15,
    min_dist=0.1,
    random_state=42
).fit_transform(X_pca_for_umap)

fig, axes = plt.subplots(1, 2, figsize=(8, 4))

axes[0].scatter(
    X_pca[:, 0],
    X_pca[:, 1],
    s=8,
    c=cell_colors,
    cmap="viridis"
)
axes[0].set_title("X PCA")
axes[0].set_xlabel("PC1")
axes[0].set_ylabel("PC2")

axes[1].scatter(
    X_umap[:, 0],
    X_umap[:, 1],
    s=8,
    c=cell_colors,
    cmap="viridis"
)
axes[1].set_title("X UMAP on PCs")
axes[1].set_xlabel("UMAP1")
axes[1].set_ylabel("UMAP2")

plt.tight_layout()
plt.savefig(output_file)
plt.close()