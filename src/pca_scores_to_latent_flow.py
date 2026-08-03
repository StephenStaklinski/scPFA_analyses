#!/usr/bin/env python3
"""Project expression onto PCA loadings and reconstruct states on a cell tree."""

import argparse

import numpy as np
import pandas as pd
from ete3 import Tree


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("expression_tsv")
    parser.add_argument("eigenvectors_tsv")
    parser.add_argument("tree_nwk")
    parser.add_argument("output_tsv")
    parser.add_argument(
        "--max-components",
        type=int,
        default=10,
        help="Maximum number of leading PCs to include as latent factors.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if args.max_components < 2:
        raise ValueError("--max-components must be at least 2")

    expression = pd.read_table(args.expression_tsv, index_col=0)
    loadings = pd.read_table(args.eigenvectors_tsv, index_col=0)
    expression.index = expression.index.astype(str)
    expression.columns = expression.columns.astype(str)
    loadings.index = loadings.index.astype(str)

    if not expression.index.is_unique:
        raise ValueError("Expression-matrix cell names are not unique")
    if not expression.columns.is_unique:
        raise ValueError("Expression-matrix gene names are not unique")
    if not loadings.index.is_unique:
        raise ValueError("PCA-loading gene names are not unique")

    components = list(loadings.columns[: args.max_components])
    if len(components) < 2:
        raise ValueError("PCA-loading table has fewer than two components")

    missing_genes = expression.columns.difference(loadings.index)
    extra_genes = loadings.index.difference(expression.columns)
    if len(missing_genes) or len(extra_genes):
        raise ValueError(
            "Expression and PCA-loading genes differ "
            f"({len(missing_genes)} missing loadings; "
            f"{len(extra_genes)} loadings absent from expression)"
        )

    ordered_loadings = loadings.loc[expression.columns, components]
    tip_scores = pd.DataFrame(
        expression.to_numpy(float) @ ordered_loadings.to_numpy(float),
        index=expression.index,
        columns=components,
    )

    tree = Tree(args.tree_nwk, format=1)
    tree_cells = {leaf.name for leaf in tree.iter_leaves()}
    shared_cells = tip_scores.index[tip_scores.index.isin(tree_cells)]
    if len(shared_cells) < 2:
        raise ValueError("Fewer than two expression cells match tree leaves")

    tip_scores = tip_scores.loc[shared_cells]
    tree.prune(shared_cells.tolist(), preserve_branch_length=True)

    nodes = list(tree.traverse("preorder"))
    node_ids = {node: i for i, node in enumerate(nodes)}
    node_names = {
        node: node.name if node.is_leaf() else f"node_{node_ids[node]}"
        for node in nodes
    }

    node_scores = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf():
            node_scores[node] = tip_scores.loc[node.name].to_numpy(float)
        else:
            child_scores = np.vstack([node_scores[child] for child in node.children])
            node_scores[node] = child_scores.mean(axis=0)

    rows = []
    for node in nodes:
        parent = node.up
        values = node_scores[node]
        row = {
            "node_id": node_ids[node],
            "parent_id": -1 if parent is None else node_ids[parent],
            "node_name": node_names[node],
            "parent_name": "NA" if parent is None else node_names[parent],
            "is_tip": int(node.is_leaf()),
            "tree_depth": tree.get_distance(node),
            "branch_length": 0.0 if parent is None else float(node.dist),
            "pc1": values[0],
            "pc2": values[1],
        }
        for component, value in zip(components, values):
            row[f"factor_{component}"] = value
        rows.append(row)

    pd.DataFrame(rows).to_csv(args.output_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
