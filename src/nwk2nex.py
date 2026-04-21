import argparse
from Bio import Phylo


def extract_taxa_from_first_tree(newick_file):
    """
    Read the first Newick tree and return the list of taxa names.
    """
    tree = next(Phylo.parse(newick_file, "newick"))
    return [t.name for t in tree.get_terminals()]

parser = argparse.ArgumentParser()
parser.add_argument("newick_file")
parser.add_argument("nexus_file")
args = parser.parse_args()

taxa = extract_taxa_from_first_tree(args.newick_file)

with open(args.newick_file) as fin, open(args.nexus_file, "w") as fout:

    fout.write("#NEXUS\n")

    fout.write("Begin taxa;\n")
    fout.write(f"  Dimensions ntax={len(taxa)};\n")
    fout.write("  TaxLabels\n\t\t" + "\n\t\t".join(taxa) + ";\n")
    fout.write("End;\n\n")

    fout.write("Begin trees;\n")

    tree_idx = 1
    for line in fin:
        line = line.strip()
        if not line:
            continue
        fout.write(f"  Tree tree_{tree_idx} = {line}\n")
        tree_idx += 1

    fout.write("End;\n")