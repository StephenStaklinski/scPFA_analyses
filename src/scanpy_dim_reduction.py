#!/usr/bin/env python

import argparse
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

sc.settings.n_jobs = 1


parser = argparse.ArgumentParser()
parser.add_argument("expr_tsv", help="Input expression TSV: rows=cells, columns=genes")
parser.add_argument("out_prefix", help="Output prefix, e.g. results/sample")
parser.add_argument("--min_genes", type=int, default=0)
parser.add_argument("--min_cells", type=int, default=0)
parser.add_argument("--n_hvg", type=int, default=2000)
parser.add_argument("--n_pcs", type=int, default=50)
parser.add_argument("--neighbors", type=int, default=15)

args = parser.parse_args()

# Read expression matrix (cells x genes)
X = pd.read_csv(args.expr_tsv, sep="\t", index_col=0)

# Build AnnData: rows = cells, columns = genes
adata = sc.AnnData(X)

# Basic QC filtering
sc.pp.filter_cells(adata, min_genes=args.min_genes)
sc.pp.filter_genes(adata, min_cells=args.min_cells)

# Preserve raw counts
adata.layers["counts"] = adata.X.copy()

# Standard preprocessing
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=args.n_hvg,
    flavor="seurat"
)

adata = adata[:, adata.var["highly_variable"]].copy()

# Scale and PCA
sc.pp.scale(adata, max_value=10)

max_pcs = min(adata.n_obs, adata.n_vars) - 1
n_pcs = min(args.n_pcs, max_pcs)

if n_pcs < 1:
    raise ValueError(
        f"Not enough cells/genes after filtering for PCA: "
        f"{adata.n_obs} cells x {adata.n_vars} genes"
    )

print(f"Running PCA with n_pcs={n_pcs} "
      f"(requested {args.n_pcs}; matrix is {adata.n_obs} x {adata.n_vars})")

sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")

# Neighbors and UMAP
n_neighbors = min(args.neighbors, adata.n_obs - 1)
if n_neighbors < 2:
    raise ValueError(
        f"Not enough cells for neighbors/UMAP: {adata.n_obs} cells"
    )
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
sc.tl.umap(adata)

# Export coordinates
pca_df = pd.DataFrame(
    adata.obsm["X_pca"],
    index=adata.obs_names,
    columns=[f"PC{i+1}" for i in range(adata.obsm["X_pca"].shape[1])]
)
pca_df.index.name = "cell"
pca_df.to_csv(f"{args.out_prefix}.pca.tsv", sep="\t")

umap_df = pd.DataFrame(
    adata.obsm["X_umap"],
    index=adata.obs_names,
    columns=["UMAP1", "UMAP2"]
)
umap_df.index.name = "cell"
umap_df.to_csv(f"{args.out_prefix}.umap.tsv", sep="\t")

# Save PCA plot
sc.pl.pca(
    adata,
    show=False,
    save=None
)
plt.savefig(f"{args.out_prefix}.pca.pdf", bbox_inches="tight")
plt.close()

# Save UMAP plot
sc.pl.umap(
    adata,
    show=False,
    save=None
)
plt.savefig(f"{args.out_prefix}.umap.pdf", bbox_inches="tight")
plt.close()

# # Save processed AnnData
# adata.write_h5ad(f"{args.out_prefix}.processed.h5ad")
