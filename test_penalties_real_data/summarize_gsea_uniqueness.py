#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize biological uniqueness of L programs from factor GSEA results."
    )
    parser.add_argument("--runs-dir", default="runs")
    parser.add_argument("--base-summary", default=None)
    parser.add_argument("--output", default="overlap_summary_gsea.tsv")
    parser.add_argument("--clone", default=None)
    parser.add_argument("--padj-cutoff", type=float, default=0.05)
    parser.add_argument("--top-n", type=int, default=20)
    parser.add_argument(
        "--ignore-direction",
        action="store_true",
        help="Compare GO term IDs only instead of direction:GO term pairs.",
    )
    return parser.parse_args()


def jaccard(a, b):
    union = a | b
    if not union:
        return None
    return len(a & b) / len(union)


def factor_term_sets(path, padj_cutoff, top_n, use_direction):
    df = pd.read_csv(path, sep="\t")
    required = {"ID", "factor", "direction", "p.adjust", "NES"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")

    df["p.adjust"] = pd.to_numeric(df["p.adjust"], errors="coerce")
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    df = df.dropna(subset=["ID", "factor", "direction", "p.adjust", "NES"])
    df = df[df["p.adjust"] <= padj_cutoff].copy()
    if df.empty:
        return {}

    df["abs_NES"] = df["NES"].abs()
    df = df.sort_values(["factor", "p.adjust", "abs_NES", "ID"],
                        ascending=[True, True, False, True])
    if top_n > 0:
        df = df.groupby("factor", sort=True).head(top_n)

    terms = {}
    for factor, factor_df in df.groupby("factor", sort=True):
        if use_direction:
            term_set = set(
                factor_df["direction"].astype(str) + ":" + factor_df["ID"].astype(str)
            )
        else:
            term_set = set(factor_df["ID"].astype(str))
        terms[str(factor)] = term_set
    return terms


def summarize_gsea(path, padj_cutoff, top_n, use_direction):
    terms = factor_term_sets(path, padj_cutoff, top_n, use_direction)
    factors = sorted(terms)
    similarities = []
    for i, factor_a in enumerate(factors):
        for factor_b in factors[i + 1:]:
            sim = jaccard(terms[factor_a], terms[factor_b])
            if sim is not None:
                similarities.append(sim)

    if similarities:
        mean_jaccard = sum(similarities) / len(similarities)
        uniqueness = 1.0 - mean_jaccard
    else:
        mean_jaccard = float("nan")
        uniqueness = float("nan")

    term_counts = [len(terms[factor]) for factor in factors]
    return {
        "gsea_mean_pairwise_jaccard": mean_jaccard,
        "gsea_uniqueness": uniqueness,
        "gsea_n_pairwise_comparisons": len(similarities),
        "gsea_n_enriched_factors": len(factors),
        "gsea_mean_terms_per_enriched_factor": (
            sum(term_counts) / len(term_counts) if term_counts else float("nan")
        ),
    }


def parse_clone_condition(path):
    prefix = path.with_suffix("")
    prefix = Path(str(prefix).removesuffix(".L.gsea"))
    return prefix.parent.name, prefix.name


def main():
    args = parse_args()
    runs_dir = Path(args.runs_dir)
    rows = []

    for gsea_path in sorted(runs_dir.glob("*/*.L.gsea.tsv")):
        clone, condition = parse_clone_condition(gsea_path)
        if args.clone and clone != args.clone:
            continue
        try:
            metrics = summarize_gsea(
                gsea_path,
                args.padj_cutoff,
                args.top_n,
                not args.ignore_direction,
            )
        except ValueError as exc:
            sys.stderr.write(f"Skipping {gsea_path}: {exc}\n")
            continue
        rows.append({"clone": clone, "condition": condition, **metrics})

    if not rows:
        sys.stderr.write(f"No GSEA result files found under {runs_dir}\n")
        return 1

    gsea_df = pd.DataFrame(rows)
    if args.base_summary:
        base_df = pd.read_csv(args.base_summary, sep="\t")
        out = base_df.merge(gsea_df, on=["clone", "condition"], how="left")
    else:
        out = gsea_df

    keep_cols = [
        "clone",
        "condition",
        "mean_offdiag_abs_L_correlation",
        "mean_offdiag_L_correlation",
        "final_objective",
        "brownian_negll",
        "observation_negll",
        "gsea_uniqueness",
        "gsea_mean_pairwise_jaccard",
        "gsea_n_enriched_factors",
        "gsea_mean_terms_per_enriched_factor",
    ]
    keep_cols = [col for col in keep_cols if col in out.columns]
    out = out[keep_cols]

    out.to_csv(args.output, sep="\t", index=False)
    print(f"Saved GSEA uniqueness summary to: {args.output}")
    print(f"Summarized {len(gsea_df)} GSEA run(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
