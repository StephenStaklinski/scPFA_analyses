#!/usr/bin/env python3

import argparse
import math
import os
from pathlib import Path
import re
import sys

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import pandas as pd


CONDITION_PREFIXES = ["Ll", "Lo", "Rm", "Fa", "K", "T"]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Summarize likelihood, top-gene Jaccard, and L-row Pearson metrics "
            "from completed penalty sweep runs."
        )
    )
    parser.add_argument("--runs-dir", default="runs")
    parser.add_argument("--output", required=True)
    parser.add_argument("--clone", default=None)
    parser.add_argument("--conditions", nargs="*", default=None)
    parser.add_argument("--top-n", type=int, default=10)
    return parser.parse_args()


def parse_condition(condition):
    parts = {}
    for token in condition.split("_"):
        for prefix in CONDITION_PREFIXES:
            if token.startswith(prefix):
                parts[prefix] = token[len(prefix):]
                break
    return parts


def condition_columns(condition):
    parts = parse_condition(condition)
    return {
        "tree_mode": parts.get("T", ""),
        "remove_ribo_mito": parts.get("Rm", ""),
        "final_absorbing_factor": parts.get("Fa", ""),
        "dim": parts.get("K", ""),
        "L_l1_strength": parts.get("Ll", ""),
        "L_loading_overlap_strength": parts.get("Lo", ""),
    }


def best_fit_values_from_log(path):
    header = None
    best_row = None
    best_objective = None

    def row_values(row, source):
        brownian_negll = parse_float(row.get("brownian_neglprior"))
        observation_negll = parse_float(row.get("observation_negll"))
        combined_log_likelihood = float("nan")
        if not math.isnan(brownian_negll) and not math.isnan(observation_negll):
            combined_log_likelihood = -(brownian_negll + observation_negll)
        return {
            "final_objective": parse_float(row.get("objective")),
            "brownian_negll": brownian_negll,
            "observation_negll": observation_negll,
            "brownian_log_likelihood": -brownian_negll,
            "observation_log_likelihood": -observation_negll,
            "combined_log_likelihood": combined_log_likelihood,
            "fit_metric_source": source,
        }

    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if fields[0] == "step":
                header = fields
                continue
            if fields[0] == "# best_state" and header is not None:
                return row_values(dict(zip(header, fields[1:])), "best_state")
            if line.startswith("#") or header is None or not fields[0].isdigit():
                continue
            row = dict(zip(header, fields))
            objective = parse_float(row.get("objective"))
            if math.isnan(objective):
                continue
            if best_objective is None or objective < best_objective:
                best_objective = objective
                best_row = row

    if best_row is None:
        return {
            "final_objective": float("nan"),
            "brownian_negll": float("nan"),
            "observation_negll": float("nan"),
            "brownian_log_likelihood": float("nan"),
            "observation_log_likelihood": float("nan"),
            "combined_log_likelihood": float("nan"),
            "fit_metric_source": "NA",
        }
    return row_values(best_row, "best_sample")


def parse_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def pearson_metrics(path):
    corr = pd.read_csv(path, sep="\t", index_col=0)
    values = []
    for row_idx in range(corr.shape[0]):
        for col_idx in range(corr.shape[1]):
            if row_idx != col_idx:
                values.append(float(corr.iloc[row_idx, col_idx]))
    if not values:
        return {
            "mean_offdiag_abs_L_pearson": float("nan"),
            "mean_offdiag_L_pearson": float("nan"),
            "max_offdiag_abs_L_pearson": float("nan"),
            "n_offdiag_L_pearson": 0,
        }
    abs_values = [abs(v) for v in values]
    return {
        "mean_offdiag_abs_L_pearson": sum(abs_values) / len(abs_values),
        "mean_offdiag_L_pearson": sum(values) / len(values),
        "max_offdiag_abs_L_pearson": max(abs_values),
        "n_offdiag_L_pearson": len(values),
    }


def top_gene_columns(columns, direction, top_n):
    matches = []
    pattern = re.compile(rf"^gene_{direction}_(\d+)$")
    for col in columns:
        match = pattern.match(col)
        if match:
            rank = int(match.group(1))
            if rank <= top_n:
                matches.append((rank, col))
    return [col for _, col in sorted(matches)]


def factor_gene_sets(path, top_n):
    df = pd.read_csv(path, sep="\t")
    if "factor" not in df.columns:
        raise ValueError(f"{path} missing column: factor")

    up_cols = top_gene_columns(df.columns, "up", top_n)
    down_cols = top_gene_columns(df.columns, "down", top_n)
    if not up_cols and not down_cols:
        raise ValueError(f"{path} has no gene_up_* or gene_down_* columns")

    up_sets = {}
    down_sets = {}
    signed_sets = {}
    unsigned_sets = {}
    for _, row in df.iterrows():
        factor = str(row["factor"])
        up = {
            str(row[col])
            for col in up_cols
            if pd.notna(row[col]) and str(row[col]) != ""
        }
        down = {
            str(row[col])
            for col in down_cols
            if pd.notna(row[col]) and str(row[col]) != ""
        }
        up_sets[factor] = up
        down_sets[factor] = down
        signed_sets[factor] = {f"up:{gene}" for gene in up} | {
            f"down:{gene}" for gene in down
        }
        unsigned_sets[factor] = up | down
    return up_sets, down_sets, signed_sets, unsigned_sets


def jaccard(a, b):
    union = a | b
    if not union:
        return float("nan")
    return len(a & b) / len(union)


def pairwise_jaccard_values(gene_sets):
    factors = sorted(gene_sets)
    values = []
    for i, factor_a in enumerate(factors):
        for factor_b in factors[i + 1 :]:
            value = jaccard(gene_sets[factor_a], gene_sets[factor_b])
            if not math.isnan(value):
                values.append(value)
    return values


def mean_pairwise_jaccard(gene_sets):
    values = pairwise_jaccard_values(gene_sets)
    if not values:
        return float("nan"), 0
    return sum(values) / len(values), len(values)


def top_gene_metrics(path, top_n):
    up_sets, down_sets, signed_sets, unsigned_sets = factor_gene_sets(path, top_n)
    up_jaccard, n_pairs = mean_pairwise_jaccard(up_sets)
    down_jaccard, _ = mean_pairwise_jaccard(down_sets)
    signed_jaccard, _ = mean_pairwise_jaccard(signed_sets)
    unsigned_jaccard, _ = mean_pairwise_jaccard(unsigned_sets)
    signed_values = pairwise_jaccard_values(signed_sets)
    max_signed_jaccard = max(signed_values) if signed_values else float("nan")
    return {
        "top_gene_top_n": top_n,
        "top_gene_mean_pairwise_signed_jaccard": signed_jaccard,
        "top_gene_max_pairwise_signed_jaccard": max_signed_jaccard,
        "top_gene_mean_pairwise_unsigned_jaccard": unsigned_jaccard,
        "top_gene_mean_pairwise_up_jaccard": up_jaccard,
        "top_gene_mean_pairwise_down_jaccard": down_jaccard,
        "top_gene_n_pairwise_comparisons": n_pairs,
        "top_gene_n_factors": len(signed_sets),
    }


def main():
    args = parse_args()
    if args.top_n <= 0:
        sys.stderr.write("--top-n must be positive.\n")
        return 1

    runs_dir = Path(args.runs_dir)
    conditions = set(args.conditions) if args.conditions else None
    rows = []
    for log_path in sorted(runs_dir.glob("*/*.log")):
        prefix = log_path.with_suffix("")
        clone = prefix.parent.name
        if args.clone and clone != args.clone:
            continue
        if conditions is not None and prefix.name not in conditions:
            continue

        corr_path = Path(f"{prefix}.L.pearson_correlation.tsv")
        top_genes_path = Path(f"{prefix}.L.top_genes.tsv")
        if not corr_path.exists() or not top_genes_path.exists():
            sys.stderr.write(f"Skipping incomplete run: {prefix}\n")
            continue

        try:
            gene_metrics = top_gene_metrics(top_genes_path, args.top_n)
        except ValueError as exc:
            sys.stderr.write(f"Skipping {top_genes_path}: {exc}\n")
            continue

        rows.append(
            {
                "clone": clone,
                "condition": prefix.name,
                **condition_columns(prefix.name),
                **best_fit_values_from_log(log_path),
                **pearson_metrics(corr_path),
                **gene_metrics,
            }
        )

    if not rows:
        sys.stderr.write(f"No complete runs found under {runs_dir}\n")
        return 1

    out = pd.DataFrame(rows)
    out = out.sort_values(["clone", "combined_log_likelihood", "condition"])
    out.to_csv(args.output, sep="\t", index=False)
    print(f"Saved selection metrics to: {args.output}")
    print(f"Summarized {len(out)} run(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
