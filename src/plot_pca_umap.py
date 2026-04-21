#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import umap

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

f_tsv = sys.argv[1]
outprefix = sys.argv[2]

# Read F matrix: rows = cells, columns = latent factors
df = pd.read_csv(f_tsv, sep="\t", index_col=0)
F = df.to_numpy()

n_cells = F.shape[0]
cell_colors = np.arange(n_cells)

# PCA
pca = PCA(n_components=2)
F_pca = pca.fit_transform(F)

plt.figure(figsize=(4, 3.5))
plt.scatter(F_pca[:, 0], F_pca[:, 1], s=8, c=cell_colors, cmap='viridis')
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.tight_layout()
plt.savefig(f"{outprefix}.pca.pdf")
plt.close()

# UMAP
umap_model = umap.UMAP(n_components=2, random_state=42)
F_umap = umap_model.fit_transform(F)

plt.figure(figsize=(4, 3.5))
plt.scatter(F_umap[:, 0], F_umap[:, 1], s=8, c=cell_colors, cmap='viridis')
plt.xlabel("UMAP1")
plt.ylabel("UMAP2")
plt.tight_layout()
plt.savefig(f"{outprefix}.umap.pdf")
plt.close()