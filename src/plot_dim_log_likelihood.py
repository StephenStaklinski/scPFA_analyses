#!/usr/bin/env python3

import argparse
import math
from pathlib import Path
import re
import sys

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import seaborn as sns


sns.set_theme(style="ticks", context="paper")
plt.rcParams.update(
    {
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "legend.fontsize": 8,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "text.color": "black",
        "axes.edgecolor": "black",
        "axes.labelcolor": "black",
        "xtick.color": "black",
        "ytick.color": "black",
    }
)

PENALTY_COLORS = ["#F0E442", "#CC79A7", "#009E73", "#56B4E9", "#E69F00"]
#                   yellow     pink       teal      sky-blue   amber
OBS_COLOR      = PENALTY_COLORS[2]   # teal     — observation log-lik
COMBINED_COLOR = PENALTY_COLORS[3]   # sky-blue — combined log-lik
COUNT_COLOR    = PENALTY_COLORS[4]   # amber    — GO factor count
NEW_GO_COLOR   = PENALTY_COLORS[1]   # pink     — new GO terms per factor
RUNTIME_COLOR  = PENALTY_COLORS[0]   # yellow   — runtime

BROWNIAN_PLOT_THRESHOLDS = [0.25, 0.50, 0.75, 0.99]
BROWNIAN_THRESHOLD_COLORS = [
    PENALTY_COLORS[3],  # sky-blue — 25 %
    PENALTY_COLORS[2],  # teal     — 50 %
    PENALTY_COLORS[4],  # amber    — 75 %
    PENALTY_COLORS[1],  # pink     — 99 %
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot dimensionality sweep diagnostics."
    )
    parser.add_argument(
        "--gsea-padj-cutoff",
        type=float,
        default=0.05,
        help="Adjusted p-value cutoff for counting GO-enriched factors.",
    )
    parser.add_argument("output_pdf")
    parser.add_argument("logs", nargs="+")
    return parser.parse_args()


def parse_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def dim_from_path(path):
    match = re.search(r"\.K([0-9]+)\.fit\.log$", path.name)
    if not match:
        raise ValueError(f"Could not parse dimension from log filename: {path}")
    return int(match.group(1))


def row_values(row, source):
    observation_negll = parse_float(row.get("observation_negll"))
    brownian_negll = parse_float(row.get("brownian_neglprior"))
    combined_log_likelihood = float("nan")
    if not math.isnan(observation_negll) and not math.isnan(brownian_negll):
        combined_log_likelihood = -(observation_negll + brownian_negll)
    return {
        "observation_negll": observation_negll,
        "brownian_negll": brownian_negll,
        "observation_log_likelihood": -observation_negll,
        "brownian_log_likelihood": -brownian_negll,
        "combined_log_likelihood": combined_log_likelihood,
        "objective": parse_float(row.get("objective")),
        "source": source,
    }


def best_fit_values_from_log(path):
    header = None
    best_row = None
    best_objective = None

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
        raise ValueError(f"No fit rows found in {path}")
    return row_values(best_row, "best_sample")


def n_factors_for_threshold(sigma2_latent, total_sigma2_latent, threshold):
    if total_sigma2_latent <= 0:
        return float("nan")
    cumulative = 0.0
    for i, value in enumerate(sigma2_latent, start=1):
        cumulative += value
        if cumulative / total_sigma2_latent >= threshold:
            return i
    return float("nan")


def values_from_summary(path):
    df = pd.read_csv(path, sep="\t")
    if not {"parameter", "value"}.issubset(df.columns):
        raise ValueError(f"Summary file has unexpected columns: {path}")
    values = dict(zip(df["parameter"], df["value"]))
    sigma2_latent = sorted(
        [
            parse_float(v)
            for p, v in values.items()
            if re.fullmatch(r"sigma2_latent_LF[0-9]+", p)
            and math.isfinite(parse_float(v))
            and parse_float(v) > 0
        ],
        reverse=True,
    )
    total_sigma2_latent = sum(sigma2_latent)
    result = {
        "n_cells": int(float(values["n_cells"])),
        "n_genes": int(float(values["n_genes"])),
        "summary_k": int(float(values["k"])),
        "total_sigma2_latent": total_sigma2_latent,
    }
    for t in BROWNIAN_PLOT_THRESHOLDS:
        pct = int(round(t * 100))
        result[f"n_factors_bv_{pct}pct"] = n_factors_for_threshold(
            sigma2_latent, total_sigma2_latent, t
        )
    return result


def summary_path_from_log(path):
    return path.with_suffix(".summary.tsv")


def gsea_path_from_log(path):
    return path.with_suffix(".L.gsea.tsv")


def time_path_from_log(path):
    return path.with_suffix(".time")


def runtime_seconds_from_time_file(path):
    if not path.exists():
        return float("nan")
    with path.open() as fh:
        content = fh.read()
    match = re.search(r"(\d+):(\d+(?:\.\d+)?)elapsed", content)
    if not match:
        return float("nan")
    return int(match.group(1)) * 60 + float(match.group(2))


def gsea_values(path, k, padj_cutoff):
    if not path.exists():
        return {
            "gsea_padj_cutoff": padj_cutoff,
            "n_go_enriched_factors": float("nan"),
            "pct_go_enriched_factors": float("nan"),
        }

    df = pd.read_csv(path, sep="\t")
    required = {"factor", "p.adjust"}
    if not required.issubset(df.columns):
        raise ValueError(f"GSEA file has unexpected columns: {path}")

    sig = df[pd.to_numeric(df["p.adjust"], errors="coerce") < padj_cutoff]
    n_enriched = sig["factor"].nunique()
    return {
        "gsea_padj_cutoff": padj_cutoff,
        "n_go_enriched_factors": n_enriched,
        "pct_go_enriched_factors": 100.0 * n_enriched / k if k > 0 else float("nan"),
    }


def compute_new_go_terms(sorted_df, padj_cutoff):
    """Return a list of new-GO-term counts in the same row order as sorted_df."""
    seen = set()
    result = []
    for _, row in sorted_df.iterrows():
        gsea_path = gsea_path_from_log(Path(row["log_path"]))
        if not gsea_path.exists():
            result.append(float("nan"))
            continue
        gdf = pd.read_csv(gsea_path, sep="\t")
        if not {"ID", "p.adjust"}.issubset(gdf.columns):
            result.append(float("nan"))
            continue
        sig_ids = set(
            gdf.loc[pd.to_numeric(gdf["p.adjust"], errors="coerce") < padj_cutoff, "ID"]
        )
        result.append(len(sig_ids - seen))
        seen |= sig_ids
    return result


def set_integer_k_axis(ax, values):
    values = pd.Series(values).dropna()
    if values.empty:
        return

    span = values.max() - values.min()
    if span >= 20:
        spacing = 5
    elif span >= 10:
        spacing = 2
    else:
        spacing = 1

    ax.xaxis.set_major_locator(mticker.MultipleLocator(spacing))
    ax.xaxis.set_major_formatter(mticker.StrMethodFormatter("{x:.0f}"))


def main():
    args = parse_args()
    rows = []
    for log in args.logs:
        path = Path(log)
        if not path.exists():
            sys.stderr.write(f"WARNING: skipping missing log: {path}\n")
            continue
        try:
            summary = summary_path_from_log(path)
            if not summary.exists():
                raise ValueError(f"Missing summary file for {path}: {summary}")
            rows.append(
                {
                    "k": dim_from_path(path),
                    "log": path.name,
                    "log_path": str(path),
                    **best_fit_values_from_log(path),
                    **values_from_summary(summary),
                    **gsea_values(
                        gsea_path_from_log(path),
                        dim_from_path(path),
                        args.gsea_padj_cutoff,
                    ),
                    "runtime_s": runtime_seconds_from_time_file(time_path_from_log(path)),
                }
            )
        except ValueError as exc:
            sys.stderr.write(f"WARNING: {exc}\n")

    if not rows:
        sys.stderr.write("No usable log files found.\n")
        return 1

    df = pd.DataFrame(rows).sort_values("k")
    output_pdf = Path(args.output_pdf)

    df["delta_k"] = df["k"].diff()
    df["delta_observation_log_likelihood"] = (
        df["observation_log_likelihood"].diff() / df["delta_k"]
    )
    df["delta_combined_log_likelihood"] = (
        df["combined_log_likelihood"].diff() / df["delta_k"]
    )
    reference = df.loc[df["k"] == 2, "delta_observation_log_likelihood"]
    if reference.empty or not math.isfinite(reference.iloc[0]) or reference.iloc[0] == 0:
        reference = df["delta_observation_log_likelihood"].dropna()
    observation_reference_gain = (
        reference.iloc[0] if not reference.empty and reference.iloc[0] != 0 else float("nan")
    )

    reference = df.loc[df["k"] == 2, "delta_combined_log_likelihood"]
    if reference.empty or not math.isfinite(reference.iloc[0]) or reference.iloc[0] == 0:
        reference = df["delta_combined_log_likelihood"].dropna()
    combined_reference_gain = (
        reference.iloc[0] if not reference.empty and reference.iloc[0] != 0 else float("nan")
    )

    df["relative_observation_log_likelihood_gain"] = (
        df["delta_observation_log_likelihood"] / observation_reference_gain
    )
    df["relative_combined_log_likelihood_gain"] = (
        df["delta_combined_log_likelihood"] / combined_reference_gain
    )
    df.to_csv(output_pdf.with_suffix(".tsv"), sep="\t", index=False)

    delta_df = df.dropna(
        subset=["delta_observation_log_likelihood", "delta_combined_log_likelihood"]
    )
    go_df = df.dropna(subset=["n_go_enriched_factors"])
    df["n_new_go_terms"] = compute_new_go_terms(df, args.gsea_padj_cutoff)
    df["new_go_terms_per_factor"] = df["n_new_go_terms"] / df["delta_k"].fillna(df["k"])
    new_go_df = df.dropna(subset=["new_go_terms_per_factor"])
    bv_cols = [f"n_factors_bv_{int(round(t * 100))}pct" for t in BROWNIAN_PLOT_THRESHOLDS]
    brownian_df = df.dropna(subset=bv_cols, how="all")
    runtime_df = df.dropna(subset=["runtime_s"])

    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5, figsize=(17.0, 3.2))

    if not delta_df.empty:
        ax1.axhline(0, linewidth=0.8, color="0.45", linestyle=":")
        observation_line, = ax1.plot(
            delta_df["k"],
            delta_df["delta_observation_log_likelihood"],
            marker="o",
            linewidth=1.4,
            color=OBS_COLOR,
            label="Observation log-likelihood",
        )
        combined_line, = ax1.plot(
            delta_df["k"],
            delta_df["delta_combined_log_likelihood"],
            marker="o",
            linewidth=1.0,
            color=COMBINED_COLOR,
            label="Observation + Brownian\nlog-likelihoods",
        )
        ax1.set_xlabel("Latent dimension (k)")
        ax1.set_ylabel("Log-likelihood gain\nper dimension")
        set_integer_k_axis(ax1, delta_df["k"])
        ax1.legend(
            [observation_line, combined_line],
            [observation_line.get_label(), combined_line.get_label()],
            frameon=False,
            fontsize=7,
        )
    else:
        ax1.set_visible(False)
        print("Skipped log-likelihood gain panel: fewer than two usable k values.")

    if not go_df.empty:
        ax2.plot(
            go_df["k"],
            go_df["n_go_enriched_factors"],
            marker="s",
            linewidth=1.4,
            color=COUNT_COLOR,
            label=f"Factors with >=1 GO term\n(padj < {args.gsea_padj_cutoff:g})",
        )
        ax2.set_xlabel("Latent dimension (k)")
        ax2.set_ylabel(f"Factors with significant\nGO enrichment (padj < {args.gsea_padj_cutoff:g})")
        set_integer_k_axis(ax2, go_df["k"])
        ax2.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
    else:
        ax2.set_visible(False)
        print("Skipped GO enrichment panel: no usable GSEA files found.")

    if not new_go_df.empty:
        ax3.plot(
            new_go_df["k"],
            new_go_df["new_go_terms_per_factor"],
            marker="s",
            linewidth=1.4,
            color=NEW_GO_COLOR,
        )
        ax3.set_xlabel("Latent dimension (k)")
        ax3.set_ylabel(f"New GO terms per\nadded factor (padj < {args.gsea_padj_cutoff:g})")
        set_integer_k_axis(ax3, new_go_df["k"])
    else:
        ax3.set_visible(False)
        print("Skipped new GO terms panel: no usable GSEA files found.")

    if not brownian_df.empty:
        for t, col, color in zip(BROWNIAN_PLOT_THRESHOLDS, bv_cols, BROWNIAN_THRESHOLD_COLORS):
            pct = int(round(t * 100))
            sub = brownian_df.dropna(subset=[col])
            if sub.empty:
                continue
            ax4.plot(
                sub["k"],
                sub[col],
                marker="o",
                linewidth=1.2,
                color=color,
                label=f"{pct}%",
            )
        ax4.set_xlabel("Latent dimension (k)")
        ax4.set_ylabel("Factors explaining X%\nof Brownian variance")
        set_integer_k_axis(ax4, brownian_df["k"])
        ax4.yaxis.set_major_locator(mticker.MaxNLocator(integer=True))
        ax4.legend(frameon=False, fontsize=7)
    else:
        ax4.set_visible(False)
        print("Skipped Brownian variance panel: no usable sigma2_latent values found.")

    if not runtime_df.empty:
        ax5.plot(
            runtime_df["k"],
            runtime_df["runtime_s"],
            marker="o",
            linewidth=1.4,
            color=RUNTIME_COLOR,
        )
        ax5.set_xlabel("Latent dimension (k)")
        ax5.set_ylabel("Runtime (s)")
        set_integer_k_axis(ax5, runtime_df["k"])
    else:
        ax5.set_visible(False)
        print("Skipped runtime panel: no usable time files found.")

    for ax in (ax1, ax2, ax3, ax4, ax5):
        if ax.get_visible():
            sns.despine(ax=ax)

    fig.tight_layout()
    fig.savefig(output_pdf)
    plt.close(fig)

    print(f"Saved dim summary plot to: {output_pdf}")
    print(f"Saved summary to: {output_pdf.with_suffix('.tsv')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
