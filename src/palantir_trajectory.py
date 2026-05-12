#!/usr/bin/env python
import argparse
import pandas as pd
import scanpy as sc
import palantir
import matplotlib.pyplot as plt

sc.settings.n_jobs = 1


parser = argparse.ArgumentParser()
parser.add_argument("expr_tsv", help="Input TSV: rows=cells, columns=genes")
parser.add_argument("out_prefix", help="Output prefix")

parser.add_argument(
    "--start_cell",
    default=None,
    help="Starting/root cell ID for Palantir (default: first cell in matrix)"
)

parser.add_argument("--min_genes", type=int, default=0)
parser.add_argument("--min_cells", type=int, default=0)
parser.add_argument("--n_hvg", type=int, default=2000)
parser.add_argument("--n_pcs", type=int, default=50)
parser.add_argument("--neighbors", type=int, default=15)
parser.add_argument("--n_components", type=int, default=10)
parser.add_argument("--num_waypoints", type=int, default=500)

args = parser.parse_args()


# Read expression matrix
X = pd.read_csv(args.expr_tsv, sep="\t", index_col=0)

# Build AnnData: rows = cells, columns = genes
adata = sc.AnnData(X)
adata.obs_names = adata.obs_names.astype(str)
adata.var_names = adata.var_names.astype(str)


# Optional filtering
if args.min_genes > 0:
    sc.pp.filter_cells(adata, min_genes=args.min_genes)

if args.min_cells > 0:
    sc.pp.filter_genes(adata, min_cells=args.min_cells)

# Preserve raw counts
adata.layers["counts"] = adata.X.copy()

# Standard preprocessing
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


# Highly variable genes
n_hvg = min(args.n_hvg, adata.n_vars)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=n_hvg,
    flavor="seurat"
)
adata = adata[:, adata.var["highly_variable"]].copy()


# Scale and PCA
sc.pp.scale(adata, max_value=10)
max_pcs = min(adata.n_obs, adata.n_vars) - 1
n_pcs = min(args.n_pcs, max_pcs)
if n_pcs < 1:
    raise ValueError(
        f"Not enough cells/genes after preprocessing for PCA: "
        f"{adata.n_obs} cells x {adata.n_vars} genes"
    )
sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")


# Neighbors and UMAP
n_neighbors = min(args.neighbors, adata.n_obs - 1)
if n_neighbors < 2:
    raise ValueError(
        f"Not enough cells for neighbors/UMAP: {adata.n_obs} cells"
    )

sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
sc.tl.umap(adata)


# Determine root/start cell
if args.start_cell is None:
    start_cell = adata.obs_names[0]
    print(f"No --start_cell provided; using first cell: {start_cell}")
else:
    start_cell = args.start_cell

if start_cell not in adata.obs_names:
    raise ValueError(
        f"Start cell '{start_cell}' not found after filtering/preprocessing. "
        f"Example valid cell: {adata.obs_names[0]}"
    )


# Palantir diffusion maps
n_components = min(args.n_components, n_pcs, adata.n_obs - 1)
if n_components < 2:
    raise ValueError(
        f"Not enough dimensions/cells for Palantir diffusion maps: "
        f"n_components={n_components}, n_pcs={n_pcs}, n_cells={adata.n_obs}"
    )

palantir.utils.run_diffusion_maps(adata, n_components=n_components)
palantir.utils.determine_multiscale_space(adata)


# Palantir trajectory
num_waypoints = min(args.num_waypoints, adata.n_obs)
palantir_knn = min(30, adata.n_obs - 1)
if palantir_knn < 2:
    raise ValueError(
        f"Not enough cells for Palantir: {adata.n_obs} cells"
    )

pr_res = palantir.core.run_palantir(
    adata,
    start_cell,
    num_waypoints=num_waypoints,
    knn=palantir_knn
)

# Store Palantir results
adata.obs["palantir_pseudotime"] = pr_res.pseudotime
adata.obs["palantir_entropy"] = pr_res.entropy

branch_probs = pr_res.branch_probs
branch_probs.to_csv(f"{args.out_prefix}.palantir.branch_probs.tsv", sep="\t")

palantir_summary = pd.DataFrame({
    "cell": adata.obs_names,
    "palantir_pseudotime": adata.obs["palantir_pseudotime"].values,
    "palantir_entropy": adata.obs["palantir_entropy"].values,
})
palantir_summary.to_csv(
    f"{args.out_prefix}.palantir.pseudotime.tsv",
    sep="\t",
    index=False
)


# # Export coordinates
# pca_df = pd.DataFrame(
#     adata.obsm["X_pca"],
#     index=adata.obs_names,
#     columns=[f"PC{i+1}" for i in range(adata.obsm["X_pca"].shape[1])]
# )
# pca_df.index.name = "cell"
# pca_df.to_csv(f"{args.out_prefix}.pca.tsv", sep="\t")

# umap_df = pd.DataFrame(
#     adata.obsm["X_umap"],
#     index=adata.obs_names,
#     columns=["UMAP1", "UMAP2"]
# )
# umap_df.index.name = "cell"
# umap_df.to_csv(f"{args.out_prefix}.umap.tsv", sep="\t")


# Plots
sc.pl.umap(
    adata,
    color="palantir_pseudotime",
    show=False
)
plt.savefig(f"{args.out_prefix}.palantir_pseudotime.umap.pdf", bbox_inches="tight")
plt.close()

sc.pl.umap(
    adata,
    color="palantir_entropy",
    show=False
)
plt.savefig(f"{args.out_prefix}.palantir_entropy.umap.pdf", bbox_inches="tight")
plt.close()

sc.pl.pca(
    adata,
    color="palantir_pseudotime",
    show=False
)
plt.savefig(f"{args.out_prefix}.palantir_pseudotime.pca.pdf", bbox_inches="tight")
plt.close()


# # Save processed AnnData
# adata.write_h5ad(f"{args.out_prefix}.palantir.processed.h5ad")
