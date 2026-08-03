#!/usr/bin/env python3
"""Test lineage autocorrelation before and after regressing out factor 1."""

import argparse

import anndata
from ete3 import Tree
import hotspot
import numpy as np
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--module-scores", required=True)
    parser.add_argument("--factor-scores", required=True)
    parser.add_argument("--tree", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--jobs", type=int, default=4)
    return parser.parse_args()


def standardize(values):
    values = values - values.mean(axis=0, keepdims=True)
    scale = values.std(axis=0, ddof=0, keepdims=True)
    scale[scale == 0] = 1
    return values / scale


def hotspot_autocorrelations(values, cells, modules, tree, n_neighbors, jobs):
    adata = anndata.AnnData(
        X=standardize(values),
        obs=pd.DataFrame({"size_factor": 1.0}, index=cells),
        var=pd.DataFrame(index=modules),
    )
    hs = hotspot.Hotspot(
        adata,
        model="none",
        tree=tree,
        umi_counts_obs_key="size_factor",
    )
    hs.create_knn_graph(weighted_graph=False, n_neighbors=n_neighbors)
    return hs.compute_autocorrelations(jobs=jobs)


def main():
    args = parse_args()

    module_scores = pd.read_table(args.module_scores, index_col=0)
    factor_scores = pd.read_table(args.factor_scores, index_col=0)
    if "factor_1" not in factor_scores.columns:
        raise ValueError("Factor-score table does not contain factor_1")

    tree = Tree(args.tree, format=1)
    tree_cells = {leaf.name for leaf in tree.iter_leaves()}
    shared_cells = (
        module_scores.index.intersection(factor_scores.index)
        .intersection(pd.Index(tree_cells))
    )
    if len(shared_cells) < 3:
        raise ValueError("Fewer than three cells are shared across inputs")

    module_scores = module_scores.loc[shared_cells]
    factor_1 = factor_scores.loc[shared_cells, "factor_1"].to_numpy(float)
    tree.prune(shared_cells.tolist(), preserve_branch_length=True)

    design = np.column_stack([np.ones(len(shared_cells)), factor_1])
    values = module_scores.to_numpy(float)
    coefficients = np.linalg.lstsq(design, values, rcond=None)[0]
    fitted = design @ coefficients
    residuals = values - fitted

    total_ss = ((values - values.mean(axis=0)) ** 2).sum(axis=0)
    residual_ss = (residuals ** 2).sum(axis=0)
    r_squared = np.divide(
        total_ss - residual_ss,
        total_ss,
        out=np.zeros_like(total_ss),
        where=total_ss > 0,
    )

    modules = module_scores.columns.astype(str)
    original = hotspot_autocorrelations(
        values,
        shared_cells,
        modules,
        tree,
        args.n_neighbors,
        args.jobs,
    ).add_prefix("original_")
    residual = hotspot_autocorrelations(
        residuals,
        shared_cells,
        modules,
        tree,
        args.n_neighbors,
        args.jobs,
    ).add_prefix("residual_")

    regression = pd.DataFrame(
        {
            "factor_1_intercept": coefficients[0],
            "factor_1_slope": coefficients[1],
            "factor_1_r_squared": r_squared,
        },
        index=modules,
    )
    regression.index.name = "Module"
    regression.join(original).join(residual).to_csv(
        args.output,
        sep="\t",
    )


if __name__ == "__main__":
    main()
