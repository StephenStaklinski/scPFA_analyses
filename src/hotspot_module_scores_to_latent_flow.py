#!/usr/bin/env python3
"""Convert terminal-cell Hotspot module scores to a latent tree table."""

import argparse

import numpy as np
import pandas as pd
from ete3 import Tree
from sklearn.decomposition import PCA


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("module_scores_tsv")
    parser.add_argument("tree_nwk")
    parser.add_argument("output_tsv")
    return parser.parse_args()


def main():
    args = parse_args()

    scores = pd.read_table(args.module_scores_tsv, index_col=0)
    scores.index = scores.index.astype(str)
    if not scores.index.is_unique:
        raise ValueError("Module-score cell names are not unique")
    if scores.shape[1] == 0:
        raise ValueError("Module-score table has no modules")

    tree = Tree(args.tree_nwk, format=1)
    tree_cells = {leaf.name for leaf in tree.iter_leaves()}
    shared_cells = scores.index[scores.index.isin(tree_cells)]
    if len(shared_cells) < 2:
        raise ValueError("Fewer than two module-score cells match tree leaves")

    scores = scores.loc[shared_cells]
    tree.prune(shared_cells.tolist(), preserve_branch_length=True)

    n_pcs = min(2, scores.shape[0], scores.shape[1])
    tip_pcs = PCA(n_components=n_pcs).fit_transform(scores.to_numpy(float))
    if n_pcs == 1:
        tip_pcs = np.column_stack([tip_pcs[:, 0], np.zeros(len(scores))])
    tip_pcs = pd.DataFrame(tip_pcs, index=scores.index, columns=["pc1", "pc2"])

    nodes = list(tree.traverse("preorder"))
    node_ids = {node: i for i, node in enumerate(nodes)}
    node_names = {
        node: node.name if node.is_leaf() else f"node_{node_ids[node]}"
        for node in nodes
    }

    rows = []
    for node in nodes:
        descendants = [leaf.name for leaf in node.iter_leaves()]
        mean_scores = scores.loc[descendants].mean(axis=0)
        mean_pcs = tip_pcs.loc[descendants].mean(axis=0)
        parent = node.up

        row = {
            "node_id": node_ids[node],
            "parent_id": -1 if parent is None else node_ids[parent],
            "node_name": node_names[node],
            "parent_name": "NA" if parent is None else node_names[parent],
            "is_tip": int(node.is_leaf()),
            "tree_depth": tree.get_distance(node),
            "branch_length": 0.0 if parent is None else float(node.dist),
            "pc1": mean_pcs["pc1"],
            "pc2": mean_pcs["pc2"],
        }
        for module, value in mean_scores.items():
            row[f"factor_{module}"] = value
        rows.append(row)

    pd.DataFrame(rows).to_csv(args.output_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
