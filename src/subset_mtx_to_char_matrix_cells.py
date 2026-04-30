#!/usr/bin/env python3

import sys
import pandas as pd
from scipy.io import mmread


char_matrix_path = sys.argv[1]
mtx_path = sys.argv[2]
barcodes_path = sys.argv[3]
features_path = sys.argv[4]
output_tsv = sys.argv[5]

# Read character matrix; first column is assumed to be cell/index names
char_df = pd.read_csv(char_matrix_path, sep="\t", index_col=0)
char_cells = set(char_df.index.astype(str))

# Read 10x-style expression inputs
X = mmread(mtx_path).tocsr()

barcodes = pd.read_csv(barcodes_path, sep="\t", header=None)[0].astype(str).tolist()
features = pd.read_csv(features_path, sep="\t", header=None)

gene_names = features.iloc[:, 1].astype(str).tolist()

# 10x mtx is usually genes x cells, so transpose to cells x genes
if X.shape[0] == len(gene_names) and X.shape[1] == len(barcodes):
    X = X.T
elif X.shape[0] == len(barcodes) and X.shape[1] == len(gene_names):
    pass
else:
    raise ValueError(
        f"Matrix shape {X.shape} does not match "
        f"{len(barcodes)} barcodes and {len(gene_names)} genes."
    )

expr_df = pd.DataFrame.sparse.from_spmatrix(
    X,
    index=barcodes,
    columns=gene_names,
)
initial_ncells, initial_ngenes = expr_df.shape
print(f"Initial expression matrix: {initial_ncells} cells x {initial_ngenes} genes")

# Keep only cells present in character matrix
expr_df = expr_df.loc[expr_df.index.isin(char_cells)]
filtered_ncells, filtered_ngenes = expr_df.shape
print(f"Filtered expression matrix after character matrix cell filtering: {filtered_ncells} cells x {filtered_ngenes} genes")

# Optional: match character matrix order where possible
expr_df = expr_df.reindex([c for c in char_df.index.astype(str) if c in expr_df.index])

# Remove any genes with zero expression across all remaining cells
expr_df = expr_df.loc[:, (expr_df != 0).any(axis=0)]
filtered_ncells, filtered_ngenes = expr_df.shape
print(f"Filtered expression matrix after removing zero-expression genes: {filtered_ncells} cells x {filtered_ngenes} genes")

# Add index row name for output TSV
expr_df.index.name = "cell"

expr_df.to_csv(output_tsv, sep="\t")
