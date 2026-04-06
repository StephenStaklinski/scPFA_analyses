#!/usr/bin/env python3

import argparse
import re
from Bio import Phylo


def read_mapping(mapping_file):
    mapping = {}

    with open(mapping_file) as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) != 2:
                raise ValueError(
                    f"Mapping file line {line_num} does not have exactly 2 tab-separated columns: {line}"
                )

            old_name, new_name = parts

            if old_name in mapping:
                raise ValueError(f"Duplicate source name in mapping file: {old_name}")

            mapping[old_name] = new_name

    return mapping


def rename_tree_terminals(tree, mapping, require_all=False):
    for terminal in tree.get_terminals():
        old_name = terminal.name

        if old_name in mapping:
            terminal.name = mapping[old_name]
        elif require_all:
            raise ValueError(f"No mapping found for tip name: {old_name}")


def extract_taxa_from_first_tree(tree):
    return [t.name for t in tree.get_terminals()]


def strip_root_length(newick):
    """
    Remove a trailing root branch length like:
        ((A:1,B:1):2,C:3):0.00000;
    ->  ((A:1,B:1):2,C:3);

    Only removes a top-level final ':number' immediately before ';'.
    """
    newick = newick.strip()
    if not newick.endswith(";"):
        return newick

    semi_idx = len(newick) - 1
    depth = 0
    colon_idx = None

    for i, ch in enumerate(newick[:-1]):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == ":" and depth == 0:
            colon_idx = i

    if colon_idx is not None:
        return newick[:colon_idx] + ";"

    return newick


def main():
    parser = argparse.ArgumentParser(
        description="Rename tip labels in a NEXUS tree file using a two-column mapping file."
    )
    parser.add_argument("input_nexus", help="Input NEXUS file")
    parser.add_argument("mapping_file", help="Tab-delimited mapping file: old_name<TAB>new_name")
    parser.add_argument("output_nexus", help="Output NEXUS file")
    parser.add_argument(
        "--require-all",
        action="store_true",
        help="Error if any tip name in the NEXUS file is missing from the mapping"
    )
    args = parser.parse_args()

    mapping = read_mapping(args.mapping_file)

    trees = list(Phylo.parse(args.input_nexus, "nexus"))
    if not trees:
        raise ValueError(f"No trees found in NEXUS file: {args.input_nexus}")

    for tree in trees:
        rename_tree_terminals(tree, mapping, require_all=args.require_all)

    for i, tree in enumerate(trees, start=1):
        if not getattr(tree, "name", None):
            tree.name = f"tree_{i}"

    taxa = extract_taxa_from_first_tree(trees[0])

    with open(args.output_nexus, "w") as fout:
        fout.write("#NEXUS\n")

        fout.write("Begin taxa;\n")
        fout.write(f"  Dimensions ntax={len(taxa)};\n")
        fout.write("  TaxLabels\n\t\t" + "\n\t\t".join(taxa) + ";\n")
        fout.write("End;\n\n")

        fout.write("Begin trees;\n")
        for i, tree in enumerate(trees, start=1):
            newick = tree.format("newick").strip()
            newick = strip_root_length(newick)
            fout.write(f"  Tree tree_{i} = {newick}\n")
        fout.write("End;\n")


if __name__ == "__main__":
    main()