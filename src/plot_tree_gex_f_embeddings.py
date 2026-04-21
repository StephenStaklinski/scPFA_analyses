#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import umap
from Bio import Phylo

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

nexus_file = sys.argv[1]
gex_tsv = sys.argv[2]
f_tsv = sys.argv[3]
output_file = sys.argv[4]


def average_brownian_covariance_from_file(file_path, cell_names):
    if file_path.endswith(".nex") or file_path.endswith(".nexus"):
        trees = list(Phylo.parse(file_path, "nexus"))
    else:
        trees = list(Phylo.parse(file_path, "newick"))

    if len(trees) == 0:
        raise ValueError(f"No trees found in file: {file_path}")

    n = len(cell_names)
    avg_cov = np.zeros((n, n), dtype=float)
    used_trees = 0

    for tree in trees:
        tips = list(tree.get_terminals())
        tip_names = {str(tip.name) for tip in tips}
        missing = [name for name in cell_names if name not in tip_names]
        if missing:
            raise ValueError(
                f"Tree is missing {len(missing)} required cell(s), "
                f"example missing: {missing[:5]}"
            )

        # Map tip names to tip objects once
        name_to_tip = {str(tip.name): tip for tip in tips}
        ordered_tips = [name_to_tip[name] for name in cell_names]

        # Precompute root-to-node distances in one traversal
        root_dists = {}

        def assign_root_dists(clade, parent_dist):
            branch_len = clade.branch_length if clade.branch_length is not None else 0.0
            cur_dist = parent_dist + branch_len
            root_dists[clade] = cur_dist
            for child in clade.clades:
                assign_root_dists(child, cur_dist)

        root_dists[tree.root] = 0.0
        for child in tree.root.clades:
            assign_root_dists(child, 0.0)

        cov = np.zeros((n, n), dtype=float)

        # Diagonal from cached root distances
        for i, tip in enumerate(ordered_tips):
            cov[i, i] = root_dists[tip]

        # Off-diagonal from MRCA root distance lookup
        for i in range(n):
            tip_i = ordered_tips[i]
            for j in range(i + 1, n):
                tip_j = ordered_tips[j]
                mrca = tree.common_ancestor(tip_i, tip_j)
                shared = root_dists[mrca]
                cov[i, j] = shared
                cov[j, i] = shared

        avg_cov += cov
        used_trees += 1

    avg_cov /= used_trees
    return avg_cov


# Read inputs
gex_df = pd.read_csv(gex_tsv, sep="\t", index_col=0)
F_df = pd.read_csv(f_tsv, sep="\t", index_col=0)

# Use shared cells only, preserving gene expression order
gex_df.index = gex_df.index.astype(str)
F_df.index = F_df.index.astype(str)
shared_cells = [c for c in gex_df.index if c in F_df.index]
if len(shared_cells) == 0:
    raise ValueError("No overlapping cell names between gene expression and F matrices.")

gex_df = gex_df.loc[shared_cells]
F_df = F_df.loc[shared_cells]

# Build average Brownian covariance across trees for these cells
tree_cov = average_brownian_covariance_from_file(nexus_file, shared_cells)
tree_cov_df = pd.DataFrame(tree_cov, index=shared_cells, columns=shared_cells)

# Convert to arrays
tree_X = tree_cov_df.to_numpy()
gex_X = gex_df.to_numpy()
F_X = F_df.to_numpy()

cell_colors = np.arange(len(shared_cells))

# Embeddings    
tree_pca = PCA(n_components=2).fit_transform(tree_X)
tree_umap = umap.UMAP(n_components=2, random_state=42).fit_transform(tree_X)

gex_pca = PCA(n_components=2).fit_transform(gex_X)
gex_umap = umap.UMAP(n_components=2, random_state=42).fit_transform(gex_X)

F_pca = PCA(n_components=2).fit_transform(F_X)
F_umap = umap.UMAP(n_components=2, random_state=42).fit_transform(F_X)


fig, axes = plt.subplots(3, 2, figsize=(8, 12))

# Row 1: tree covariance
axes[0, 0].scatter(tree_pca[:, 0], tree_pca[:, 1], s=8, c=cell_colors, cmap="viridis")
axes[0, 0].set_title("Tree covariance PCA")
axes[0, 0].set_xlabel("PC1")
axes[0, 0].set_ylabel("PC2")

axes[0, 1].scatter(tree_umap[:, 0], tree_umap[:, 1], s=8, c=cell_colors, cmap="viridis")
axes[0, 1].set_title("Tree covariance UMAP")
axes[0, 1].set_xlabel("UMAP1")
axes[0, 1].set_ylabel("UMAP2")

# Row 2: gene expression
axes[1, 0].scatter(gex_pca[:, 0], gex_pca[:, 1], s=8, c=cell_colors, cmap="viridis")
axes[1, 0].set_title("Gene expression PCA")
axes[1, 0].set_xlabel("PC1")
axes[1, 0].set_ylabel("PC2")

axes[1, 1].scatter(gex_umap[:, 0], gex_umap[:, 1], s=8, c=cell_colors, cmap="viridis")
axes[1, 1].set_title("Gene expression UMAP")
axes[1, 1].set_xlabel("UMAP1")
axes[1, 1].set_ylabel("UMAP2")

# Row 3: latent F
axes[2, 0].scatter(F_pca[:, 0], F_pca[:, 1], s=8, c=cell_colors, cmap="viridis")
axes[2, 0].set_title("Latent F PCA")
axes[2, 0].set_xlabel("PC1")
axes[2, 0].set_ylabel("PC2")

axes[2, 1].scatter(F_umap[:, 0], F_umap[:, 1], s=8, c=cell_colors, cmap="viridis")
axes[2, 1].set_title("Latent F UMAP")
axes[2, 1].set_xlabel("UMAP1")
axes[2, 1].set_ylabel("UMAP2")

plt.tight_layout()
plt.savefig(output_file)
plt.close()