#!/usr/bin/env python3

import sys
from itertools import combinations
from statistics import mean, pstdev
from Bio import Phylo


def collapse_unifurcations(tree):
    """Remove internal nodes with only one child by merging them with that child,
    accumulating branch lengths so total distances are preserved."""
    changed = True
    while changed:
        changed = False
        for clade in tree.find_clades(order="level"):
            if clade.is_terminal() or len(clade.clades) != 1:
                continue
            child = clade.clades[0]
            if clade.branch_length is not None and child.branch_length is not None:
                child.branch_length += clade.branch_length
            elif clade.branch_length is not None:
                child.branch_length = clade.branch_length
            clade.clades = child.clades
            clade.name = child.name
            changed = True
            break  # Restart after any modification


def cv(values):
    """Coefficient of variation: ratio of standard deviation to mean.
    Scale-free measure of dispersion; nan if empty or mean is zero."""
    if len(values) == 0:
        return float("nan")
    mu = mean(values)
    if mu == 0:
        return float("nan")
    return pstdev(values) / mu


def count_leaves(clade):
    """Return the number of terminal nodes (leaves) descending from a clade."""
    return len(clade.get_terminals())

def sackin_colless_single_pass(tree):
    """Sackin index: sum of the number of internal ancestors for each tip.
    Higher values indicate a more imbalanced (pectinate/ladder-like) tree.
    Lower values indicate a more balanced (symmetric) tree.
    
    Colless index: sum over all internal nodes of |L - R|, where L and R
    are the number of leaves in the left and right subtrees respectively.
    Only defined for strictly bifurcating trees; returns nan otherwise.
    Higher values = more imbalanced; 0 = perfectly balanced."""
    internal_depth = {id(tree.root): 0}
    sackin = 0
    colless = 0
    leaf_count = {}

    for clade in tree.find_clades(order="postorder"):
        if clade.is_terminal():
            leaf_count[id(clade)] = 1
        else:
            counts = [leaf_count[id(c)] for c in clade.clades]
            leaf_count[id(clade)] = sum(counts)
            if len(clade.clades) == 2:
                colless += abs(counts[0] - counts[1])

    for clade in tree.find_clades(order="preorder"):
        d = internal_depth.get(id(clade), 0)
        if clade.is_terminal():
            sackin += d
        else:
            for child in clade.clades:
                internal_depth[id(child)] = d + 1

    return sackin, colless


def cherry_count(tree):
    """Count cherries: internal nodes whose two children are both leaves.
    More cherries = more balanced / bushier tree.
    Fewer cherries = more ladder-like / imbalanced tree."""
    count = 0
    for clade in tree.find_clades(order="level"):
        if clade.is_terminal():
            continue
        if len(clade.clades) == 2 and all(child.is_terminal() for child in clade.clades):
            count += 1
    return count


def max_ladder_length_from_clade(clade):
    """Maximum ladder length: the longest uninterrupted chain of internal nodes
    where each node has exactly one internal child (the other child, if any, is a leaf).
    A long ladder indicates a highly asymmetric, pectinate subtree."""
    if clade.is_terminal() or len(clade.clades) != 2:
        return 0

    internal_children = [c for c in clade.clades if not c.is_terminal()]
    if len(internal_children) == 1:
        # One internal child: this node extends the ladder by 1
        return 1 + max_ladder_length_from_clade(internal_children[0])

    if len(internal_children) == 2:
        # Both children are internal: ladder resets; take the longer branch
        return max(max_ladder_length_from_clade(internal_children[0]),
                   max_ladder_length_from_clade(internal_children[1]))

    return 0


tree_path = sys.argv[1]
outpath = sys.argv[2]

tree = Phylo.read(tree_path, "newick")
collapse_unifurcations(tree)

# True if every internal node has exactly 2 children (fully resolved, no polytomies)
is_bifurcating = all(len(clade.clades) == 2 for clade in tree.find_clades(order="level") if not clade.is_terminal())

all_clades = list(tree.find_clades(order="level"))
tips = tree.get_terminals()
internal_nodes = [clade for clade in all_clades if not clade.is_terminal()]

# Collect all branch lengths across the whole tree (internal + terminal branches)
branch_lengths = []
for clade in all_clades:
    if clade.branch_length is not None:
        branch_lengths.append(clade.branch_length)

# Distance from the root to each tip (sum of branch lengths along the path)
root_to_tip_distances = []
for tip in tips:
    root_to_tip_distances.append(tree.distance(tree.root, tip))

# Patristic distance between every pair of tips (sum of branches on the path between them)
pairwise_tip_distances = []
for tip1, tip2 in combinations(tips, 2):
    pairwise_tip_distances.append(tree.distance(tip1, tip2))

sackin, colless = sackin_colless_single_pass(tree)
cherries = cherry_count(tree)
n_tips = len(tips)
# Fraction of tips that are part of a cherry (each cherry contributes 2 tips)
cherry_fraction_of_tips = (2.0 * cherries / n_tips) if n_tips > 0 else float("nan")
max_ladder_length = max_ladder_length_from_clade(tree.root)

# Summary statistics for individual branch lengths
branch_length_mean = mean(branch_lengths) if len(branch_lengths) > 0 else float("nan")
branch_length_sd = pstdev(branch_lengths) if len(branch_lengths) > 0 else float("nan")
branch_length_cv = cv(branch_lengths)   # Relative variability, independent of scale
branch_length_min = min(branch_lengths) if len(branch_lengths) > 0 else float("nan")
branch_length_max = max(branch_lengths) if len(branch_lengths) > 0 else float("nan")

# Summary statistics for tip-to-tip patristic distances
pairwise_distance_mean = mean(pairwise_tip_distances) if len(pairwise_tip_distances) > 0 else float("nan")
pairwise_distance_sd = pstdev(pairwise_tip_distances) if len(pairwise_tip_distances) > 0 else float("nan")
pairwise_distance_cv = cv(pairwise_tip_distances)   # Clock-like trees have low CV; rate variation increases it
pairwise_distance_min = min(pairwise_tip_distances) if len(pairwise_tip_distances) > 0 else float("nan")
pairwise_distance_max = max(pairwise_tip_distances) if len(pairwise_tip_distances) > 0 else float("nan")

header = [
    "is_bifurcating",
    "branch_length_mean",
    "branch_length_sd",
    "branch_length_cv",
    "branch_length_min",
    "branch_length_max",
    "pairwise_distance_mean",
    "pairwise_distance_sd",
    "pairwise_distance_cv",
    "pairwise_distance_min",
    "pairwise_distance_max",
    "sackin_index",
    "colless_index",
    "cherry_count",
    "cherry_fraction_of_tips",
    "max_ladder_length",
]

values = [
    str(is_bifurcating),
    str(branch_length_mean),
    str(branch_length_sd),
    str(branch_length_cv),
    str(branch_length_min),
    str(branch_length_max),
    str(pairwise_distance_mean),
    str(pairwise_distance_sd),
    str(pairwise_distance_cv),
    str(pairwise_distance_min),
    str(pairwise_distance_max),
    str(sackin),
    str(colless),
    str(cherries),
    str(cherry_fraction_of_tips),
    str(max_ladder_length),
]

with open(outpath, "w") as f:
    f.write("\t".join(header) + "\n")
    f.write("\t".join(values) + "\n")