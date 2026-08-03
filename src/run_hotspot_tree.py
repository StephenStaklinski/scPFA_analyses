#!/usr/bin/env python3
"""Identify heritable gene modules with Hotspot and a lineage tree."""

import argparse
from pathlib import Path

import anndata
import hotspot
import numpy as np
import pandas as pd
from ete3 import Tree
from scipy import sparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", required=True, type=Path)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--out-prefix", required=True, type=Path)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--min-cells", type=int, default=10)
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--gene-fdr", type=float, default=0.05)
    parser.add_argument("--module-fdr", type=float, default=0.05)
    parser.add_argument("--min-module-genes", type=int, default=50)
    parser.add_argument(
        "--max-genes",
        type=int,
        default=0,
        help="Maximum autocorrelation-ranked genes for pairwise analysis; 0 means all significant genes.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    prefix = str(args.out_prefix)

    counts = pd.read_table(args.counts, index_col=0)
    if not counts.index.is_unique:
        raise ValueError("Cell names in the count matrix are not unique")
    if not counts.columns.is_unique:
        raise ValueError("Gene names in the count matrix are not unique")
    if (counts.to_numpy() < 0).any():
        raise ValueError("The danb model requires nonnegative count data")

    tree = Tree(str(args.tree), format=1)
    tree_cells = {leaf.name for leaf in tree.iter_leaves()}
    common_cells = counts.index[counts.index.isin(tree_cells)]
    if len(common_cells) < 2:
        raise ValueError("Fewer than two count-matrix cells match tree leaves")

    missing_from_tree = counts.index.difference(pd.Index(tree_cells))
    missing_from_counts = pd.Index(tree_cells).difference(counts.index)
    print(
        f"Using {len(common_cells)} shared cells "
        f"({len(missing_from_tree)} count-only; {len(missing_from_counts)} tree-only)"
    )

    counts = counts.loc[common_cells]
    tree.prune(common_cells.tolist(), preserve_branch_length=True)

    detected = (counts.to_numpy() > 0).sum(axis=0) >= args.min_cells
    counts = counts.loc[:, detected]
    if counts.shape[1] == 0:
        raise ValueError("No genes remain after min-cells filtering")
    print(f"Using {counts.shape[1]} genes detected in at least {args.min_cells} cells")

    adata = anndata.AnnData(
        X=sparse.csr_matrix(counts.to_numpy()),
        obs=pd.DataFrame(index=counts.index.copy()),
        var=pd.DataFrame(index=counts.columns.copy()),
    )
    adata.obs["total_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()

    hs = hotspot.Hotspot(
        adata,
        model="danb",
        tree=tree,
        umi_counts_obs_key="total_counts",
    )
    hs.create_knn_graph(
        weighted_graph=False,
        n_neighbors=args.n_neighbors,
    )

    autocorrelations = hs.compute_autocorrelations(jobs=args.jobs)
    autocorrelations.to_csv(prefix + ".autocorrelations.tsv", sep="\t")

    selected = autocorrelations.loc[autocorrelations["FDR"] < args.gene_fdr]
    selected = selected.sort_values("Z", ascending=False)
    if args.max_genes > 0:
        selected = selected.head(args.max_genes)
    if len(selected) < 2:
        raise ValueError(
            f"Only {len(selected)} genes pass the autocorrelation selection; "
            "at least two are required for module analysis"
        )
    print(f"Computing local correlations for {len(selected)} heritable genes")

    local_correlations = hs.compute_local_correlations(
        selected.index,
        jobs=args.jobs,
    )
    local_correlations.to_csv(prefix + ".local_correlations.tsv", sep="\t")

    modules = hs.create_modules(
        min_gene_threshold=args.min_module_genes,
        core_only=True,
        fdr_threshold=args.module_fdr,
    )
    modules.rename("Module").to_csv(prefix + ".modules.tsv", sep="\t")

    module_scores = hs.calculate_module_scores()
    module_scores.to_csv(prefix + ".module_scores.tsv", sep="\t")

    parameters = pd.Series(
        {
            "cells": adata.n_obs,
            "genes_after_filtering": adata.n_vars,
            "selected_heritable_genes": len(selected),
            "assigned_module_genes": int((modules != -1).sum()),
            "modules": int(modules[modules != -1].nunique()),
            "model": "danb",
            "min_cells": args.min_cells,
            "n_neighbors": args.n_neighbors,
            "gene_fdr": args.gene_fdr,
            "module_fdr": args.module_fdr,
            "min_module_genes": args.min_module_genes,
            "max_genes": args.max_genes,
        },
        name="value",
    )
    parameters.to_csv(prefix + ".parameters.tsv", sep="\t", header=True)


if __name__ == "__main__":
    main()
