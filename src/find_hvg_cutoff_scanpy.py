import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import seaborn as sns

# File paths
matrix_mtx = "/home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks/gex_lineage_tests/input_datasets/quinn_mouse_5k_all_raw_deposited_data/GSM4905335_matrix.5k.mtx"
barcodes_tsv = "/home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks/gex_lineage_tests/input_datasets/quinn_mouse_5k_all_raw_deposited_data/GSM4905335_barcodes.5k.tsv"
genes_tsv = "/home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks/gex_lineage_tests/input_datasets/quinn_mouse_5k_all_raw_deposited_data/GSM4905335_genes.5k.tsv"

outfile = "/home/staklins/projects/gex_lineage_project/gex_lineage_benchmarks/gex_lineage_tests/input_datasets/quinn_mouse_5k_all_raw_deposited_data/GSM4905335_matrix.5k.mtx.hvg_cutoff_scanpy.pdf"

# Load matrix
X = scipy.io.mmread(matrix_mtx).tocsr()

# Read in barcodes (cells)
barcodes = pd.read_csv(barcodes_tsv, sep="\t", header=None)
barcodes = barcodes.iloc[:, 0].astype(str).values

# Check gex matrix orientation. If not already cells x genes, then transpose
if X.shape[0] != len(barcodes):
    X = X.T.tocsr()

# Load genes
genes = pd.read_csv(genes_tsv, sep="\t", header=None)
if genes.shape[1] >= 2:
    gene_ids = genes.iloc[:, 0].astype(str).values
    gene_names = genes.iloc[:, 1].astype(str).values
else:
    gene_ids = genes.iloc[:, 0].astype(str).values
    gene_names = genes.iloc[:, 0].astype(str).values

# Build AnnData
adata = ad.AnnData(X=X)
adata.obs_names = barcodes
adata.var_names = gene_names
adata.var["gene_ids"] = gene_ids
adata.var_names_make_unique()

# Standard preprocessing
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Find HVGs
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# Compute per-gene variance of the log-normalized data
if sp.issparse(adata.X):
    mean = np.asarray(adata.X.mean(axis=0)).ravel()
    mean_sq = np.asarray(adata.X.power(2).mean(axis=0)).ravel()
    variance = mean_sq - mean**2
else:
    variance = np.var(adata.X, axis=0)

adata.var["variance"] = variance

# Rank genes by HVG score
ranked = adata.var.sort_values("dispersions_norm", ascending=False).copy()

n_hvg = int(adata.var["highly_variable"].sum())
cutoff_idx = n_hvg - 1

cutoff_gene = ranked.index[cutoff_idx]
cutoff_variance = ranked.iloc[cutoff_idx]["variance"]

print(f"Number of highly variable genes: {n_hvg}")
print(f"Last gene at cutoff: {cutoff_gene}")
print(f"Variance of that gene: {cutoff_variance}")

# Plot variance histogram of log-normalized data
ranked_variances = ranked["variance"].values

plt.figure(figsize=(8, 5))

sns.histplot(
    ranked_variances,
    bins=100,
    kde=False,
)

# Mark HVG cutoff variance
cutoff_variance = ranked_variances[cutoff_idx]
plt.axvline(cutoff_variance, linestyle="--", color="red", label=f"2000 HVG cutoff = {cutoff_variance:.4f}")

plt.xlabel("Variance")
plt.ylabel("Number of genes")
plt.title(f"Distribution of log-normalized gene variance for all genes")
plt.legend()
plt.tight_layout()
plt.savefig(outfile)
plt.close()
