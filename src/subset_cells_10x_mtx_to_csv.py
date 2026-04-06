#!/usr/bin/env python3

import argparse
import sys

import pandas as pd
from scipy.io import mmread
from scipy import sparse


def read_tsv_column(path, colnum=0):
    vals = []
    with open(path) as f:
        for line_num, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) <= colnum or parts[colnum] == "":
                raise ValueError(f"Empty column {colnum} in {path} at line {line_num}")
            vals.append(parts[colnum])
    return vals

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Subset a 10x-style barcodes.tsv + genes.tsv + matrix.mtx to a set of cells "
            "and write a dense TSV with cell names as rows and genes as columns."
        )
    )
    parser.add_argument("barcodes_tsv")
    parser.add_argument("genes_tsv")
    parser.add_argument("matrix_mtx")
    parser.add_argument("keep_barcodes_txt")
    parser.add_argument("output_tsv")
    args = parser.parse_args()

    barcodes = read_tsv_column(args.barcodes_tsv, colnum=0)
    genes = read_tsv_column(args.genes_tsv, colnum=1)
    
    num_input_cells = len(barcodes)
    num_input_genes = len(genes)
    print(f"Read {num_input_cells} barcodes and {num_input_genes} genes in the input matrix...", file=sys.stderr)

    keep_cells = []
    with open(args.keep_barcodes_txt) as f:
        for line in f:
            name = line.strip()
            if name:
                keep_cells.append(name)
    keep_set = set(keep_cells)
    
    print(f"Read {len(keep_cells)} barcodes designated to keep...", file=sys.stderr)

    X = mmread(args.matrix_mtx)
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    else:
        X = X.tocsr()

    n_genes, n_cells = X.shape
    if n_cells != num_input_cells:
        raise ValueError(f"Matrix column count ({n_cells}) does not match number of barcodes ({num_input_cells})")
    if n_genes != num_input_genes:
        raise ValueError(f"Matrix row count ({n_genes}) does not match number of genes ({num_input_genes})")

    barcode_to_idx = {bc: i for i, bc in enumerate(barcodes)}

    keep_indices = []
    keep_found_names = []
    missing = []
    for cell in keep_cells:
        if cell in barcode_to_idx:
            keep_indices.append(barcode_to_idx[cell])
            keep_found_names.append(cell)
        else:
            missing.append(cell)

    if len(keep_indices) == 0:
        raise ValueError("None of the requested cells were found in barcodes.tsv")

    if missing:
        print(f"Warning: {len(missing)} requested cells were not found and will be skipped: {' '.join(missing)}")

    print(f"Keeping {len(keep_indices)} cells...", file=sys.stderr)

    # Subset columns (cells)
    X_sub = X[:, keep_indices]  # genes x kept_cells

    # Find genes (colulmns) with any nonzero counts across kept cells
    gene_nonzero = (X_sub.sum(axis=1) > 0)
    gene_nonzero = pd.Series(gene_nonzero.A1 if hasattr(gene_nonzero, "A1") else gene_nonzero).astype(bool).values
    kept_gene_indices = [i for i, keep in enumerate(gene_nonzero) if keep]
    kept_genes = [genes[i] for i in kept_gene_indices]

    print(f"Keeping {len(kept_gene_indices)} genes with nonzero values...", file=sys.stderr)

    X_sub = X_sub[kept_gene_indices, :]   # kept_genes x kept_cells

    # Convert to dense cells x genes for TSV output
    dense = X_sub.T.toarray()  # kept_cells x kept_genes
    df = pd.DataFrame(dense, index=keep_found_names, columns=kept_genes)
    df.index.name = "barcode"
    df.to_csv(args.output_tsv, sep="\t")


if __name__ == "__main__":
    main()
    