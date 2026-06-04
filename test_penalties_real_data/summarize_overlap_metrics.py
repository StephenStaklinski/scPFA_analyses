#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize L-correlation overlap and fit metrics from completed runs."
    )
    parser.add_argument("--runs-dir", default="runs")
    parser.add_argument("--output", default="overlap_summary.tsv")
    return parser.parse_args()


def best_fit_values_from_log(path):
    values = {
        "objective": "NA",
        "brownian_negll": "NA",
        "observation_negll": "NA",
        "source": "NA",
    }
    header = None
    best_row = None
    best_objective = None

    def row_values(row, source):
        return {
            "objective": row.get("objective", "NA"),
            "brownian_negll": row.get("brownian_neglprior", "NA"),
            "observation_negll": row.get("observation_negll", "NA"),
            "source": source,
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
                row = dict(zip(header, fields[1:]))
                return row_values(row, "best_state")
            if line.startswith("#") or header is None:
                continue
            if fields and fields[0].isdigit():
                row = dict(zip(header, fields))
                try:
                    objective = float(row["objective"])
                except (KeyError, TypeError, ValueError):
                    continue
                if best_objective is None or objective < best_objective:
                    best_objective = objective
                    best_row = row

    if best_row is not None:
        values = row_values(best_row, "best_sample")
    return values


def correlation_metrics(path):
    corr = pd.read_csv(path, sep="\t", index_col=0)
    if corr.empty:
        return "NA", "NA", 0

    values = []
    for row_idx in range(corr.shape[0]):
        for col_idx in range(corr.shape[1]):
            if row_idx != col_idx:
                values.append(float(corr.iloc[row_idx, col_idx]))

    if not values:
        return "NA", "NA", 0
    mean_abs = sum(abs(v) for v in values) / len(values)
    mean_signed = sum(values) / len(values)
    return f"{mean_abs:.17g}", f"{mean_signed:.17g}", len(values)


def main():
    args = parse_args()
    runs_dir = Path(args.runs_dir)
    rows = []

    for corr_path in sorted(runs_dir.glob("*/*.L.pearson_correlation.tsv")):
        if ".pca." in corr_path.name:
            continue
        prefix = corr_path.with_suffix("")
        prefix = Path(str(prefix).removesuffix(".L.pearson_correlation"))
        log_path = Path(f"{prefix}.log")
        if not log_path.exists():
            sys.stderr.write(f"Skipping incomplete run: {prefix}\n")
            continue

        mean_abs, mean_signed, n = correlation_metrics(corr_path)
        values = best_fit_values_from_log(log_path)
        rows.append(
            {
                "clone": prefix.parent.name,
                "condition": prefix.name,
                "mean_offdiag_abs_L_correlation": mean_abs,
                "mean_offdiag_L_correlation": mean_signed,
                "n_offdiag_entries": n,
                "final_objective": values["objective"],
                "brownian_negll": values["brownian_negll"],
                "observation_negll": values["observation_negll"],
                "fit_metric_source": values["source"],
            }
        )

    if not rows:
        sys.stderr.write(f"No complete runs found under {runs_dir}\n")
        return 1

    out = pd.DataFrame(rows)
    out.to_csv(args.output, sep="\t", index=False)

    for _, row in out.iterrows():
        metric_path = runs_dir / row["clone"] / f"{row['condition']}.overlap.tsv"
        row.to_frame().T.to_csv(metric_path, sep="\t", index=False, header=False)

    print(f"Saved summary to: {args.output}")
    print(f"Summarized {len(out)} run(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
