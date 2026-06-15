#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize overlap of top up/down genes across learned L programs."
    )
    parser.add_argument("--runs-dir", default="runs")
    parser.add_argument("--base-summary", default=None)
    parser.add_argument("--output", default="overlap_summary_top_genes.tsv")
    parser.add_argument("--clone", default=None)
    parser.add_argument(
        "--ignore-direction",
        action="store_true",
        help="Compare gene symbols only instead of up/down-signed gene entries.",
    )
    return parser.parse_args()


def jaccard(a, b):
    union = a | b
    if not union:
        return None
    return len(a & b) / len(union)


def factor_gene_sets(path, use_direction):
    df = pd.read_csv(path, sep="\t")
    if "factor" not in df.columns:
        raise ValueError(f"{path} missing column: factor")

    up_cols = [col for col in df.columns if col.startswith("gene_up_")]
    down_cols = [col for col in df.columns if col.startswith("gene_down_")]
    if not up_cols and not down_cols:
        raise ValueError(f"{path} has no gene_up_* or gene_down_* columns")

    gene_sets = {}
    for _, row in df.iterrows():
        factor = str(row["factor"])
        values = set()
        for col in up_cols:
            gene = row[col]
            if pd.notna(gene) and str(gene) != "":
                values.add(str(gene) if not use_direction else f"up:{gene}")
        for col in down_cols:
            gene = row[col]
            if pd.notna(gene) and str(gene) != "":
                values.add(str(gene) if not use_direction else f"down:{gene}")
        gene_sets[factor] = values
    return gene_sets


def summarize_top_genes(path, use_direction):
    gene_sets = factor_gene_sets(path, use_direction)
    factors = sorted(gene_sets)
    similarities = []
    for i, factor_a in enumerate(factors):
        for factor_b in factors[i + 1:]:
            sim = jaccard(gene_sets[factor_a], gene_sets[factor_b])
            if sim is not None:
                similarities.append(sim)

    if similarities:
        mean_jaccard = sum(similarities) / len(similarities)
        uniqueness = 1.0 - mean_jaccard
    else:
        mean_jaccard = float("nan")
        uniqueness = float("nan")

    set_sizes = [len(gene_sets[factor]) for factor in factors]
    return {
        "top_gene_mean_pairwise_jaccard": mean_jaccard,
        "top_gene_uniqueness": uniqueness,
        "top_gene_n_pairwise_comparisons": len(similarities),
        "top_gene_n_factors": len(factors),
        "top_gene_mean_genes_per_factor": (
            sum(set_sizes) / len(set_sizes) if set_sizes else float("nan")
        ),
    }


def parse_clone_condition(path):
    prefix = path.with_suffix("")
    prefix = Path(str(prefix).removesuffix(".L.top_genes"))
    return prefix.parent.name, prefix.name


def main():
    args = parse_args()
    runs_dir = Path(args.runs_dir)
    rows = []

    for top_genes_path in sorted(runs_dir.glob("*/*.L.top_genes.tsv")):
        clone, condition = parse_clone_condition(top_genes_path)
        if args.clone and clone != args.clone:
            continue
        try:
            metrics = summarize_top_genes(top_genes_path, not args.ignore_direction)
        except ValueError as exc:
            sys.stderr.write(f"Skipping {top_genes_path}: {exc}\n")
            continue
        rows.append({"clone": clone, "condition": condition, **metrics})

    if not rows:
        sys.stderr.write(f"No top-gene files found under {runs_dir}\n")
        return 1

    top_gene_df = pd.DataFrame(rows)
    if args.base_summary:
        base_df = pd.read_csv(args.base_summary, sep="\t")
        out = base_df.merge(top_gene_df, on=["clone", "condition"], how="left")
    else:
        out = top_gene_df

    keep_cols = [
        "clone",
        "condition",
        "mean_offdiag_abs_L_correlation",
        "mean_offdiag_L_correlation",
        "final_objective",
        "brownian_negll",
        "observation_negll",
        "top_gene_mean_pairwise_jaccard",
        "top_gene_uniqueness",
        "top_gene_n_factors",
        "top_gene_mean_genes_per_factor",
    ]
    keep_cols = [col for col in keep_cols if col in out.columns]
    out = out[keep_cols]

    out.to_csv(args.output, sep="\t", index=False)
    print(f"Saved top-gene Jaccard summary to: {args.output}")
    print(f"Summarized {len(top_gene_df)} top-gene run(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
