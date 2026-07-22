#!/usr/bin/env python3

import argparse
import math
import os
from pathlib import Path
import sys

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib.pyplot as plt
import numpy as np
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

PENALTY_VALUES = [0.0, 0.01, 0.1, 1.0, 10.0]
DISPLAY_L1_VALUES = PENALTY_VALUES
RECOMMENDED_L1 = 0.1
RECOMMENDED_OVERLAP = 1.0
PENALTY_COLORS = ["#F0E442", "#CC79A7", "#009E73", "#56B4E9", "#E69F00"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Make cross-dataset paper figures for the penalty sweep."
    )
    parser.add_argument(
        "--figure",
        choices=["tradeoff", "normalization"],
        required=True,
    )
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--selection-summaries",
        nargs="+",
        required=True,
        help="Per-clone selection_metrics TSV files.",
    )
    parser.add_argument(
        "--runs-dir",
        default="runs",
        help="Directory containing per-clone run outputs and summary TSVs.",
    )
    parser.add_argument(
        "--gsea-padj-cutoff",
        type=float,
        default=0.05,
        help="Adjusted-p-value cutoff used to define significant GO terms.",
    )
    return parser.parse_args()


def format_penalty(value):
    value = float(value)
    return "0" if value == 0 else f"{value:g}"


def read_run_size(path):
    summary = pd.read_csv(path, sep="\t")
    values = dict(zip(summary["parameter"], summary["value"]))
    return int(float(values["n_cells"])), int(float(values["n_genes"]))


def find_run_size(runs_dir, clone):
    paths = sorted((runs_dir / clone).glob("*.summary.tsv"))
    if not paths:
        raise FileNotFoundError(f"No run summary TSV found for {clone}")
    return read_run_size(paths[0])


def load_data(selection_paths, runs_dir):
    frames = [pd.read_csv(path, sep="\t") for path in selection_paths]
    df = pd.concat(frames, ignore_index=True)

    numeric_columns = [
        "L_l1_strength",
        "L_loading_overlap_strength",
        "combined_log_likelihood",
        "mean_offdiag_abs_L_pearson",
        "top_gene_mean_pairwise_signed_jaccard",
        "top_gene_mean_pairwise_unsigned_jaccard",
    ]
    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    df = df.dropna(subset=numeric_columns).copy()

    sizes = {}
    for clone in sorted(df["clone"].unique()):
        sizes[clone] = find_run_size(runs_dir, clone)
    df["n_cells"] = df["clone"].map(lambda clone: sizes[clone][0])
    df["n_genes"] = df["clone"].map(lambda clone: sizes[clone][1])
    df["n_observations"] = df["n_cells"] * df["n_genes"]

    baseline = (
        df[
            (df["L_l1_strength"] == 0)
            & (df["L_loading_overlap_strength"] == 0)
        ]
        .set_index("clone")["combined_log_likelihood"]
    )
    missing = sorted(set(df["clone"]) - set(baseline.index))
    if missing:
        raise ValueError(f"Missing unpenalized baseline for: {', '.join(missing)}")

    df["log_likelihood_change_per_observation"] = df.apply(
        lambda row: (
            row["combined_log_likelihood"] - baseline.loc[row["clone"]]
        )
        / row["n_observations"],
        axis=1,
    )
    return df


def mean_pairwise_jaccard(factor_term_sets):
    enriched = {factor: terms for factor, terms in factor_term_sets.items() if terms}
    factors = sorted(enriched)
    values = []
    for index, factor_a in enumerate(factors):
        for factor_b in factors[index + 1:]:
            terms_a = enriched[factor_a]
            terms_b = enriched[factor_b]
            union = terms_a | terms_b
            values.append(len(terms_a & terms_b) / len(union) if union else float("nan"))
    valid = [value for value in values if not math.isnan(value)]
    return sum(valid) / len(valid) if valid else float("nan")


def load_gsea_data(selection_df, runs_dir, padj_cutoff):
    condition_penalties = (
        selection_df[
            [
                "clone",
                "condition",
                "dim",
                "L_l1_strength",
                "L_loading_overlap_strength",
            ]
        ]
        .drop_duplicates()
        .set_index(["clone", "condition"])
    )
    rows = []
    for clone in sorted(selection_df["clone"].unique()):
        for path in sorted((runs_dir / clone).glob("*.L.gsea.tsv")):
            condition = path.name.removesuffix(".L.gsea.tsv")
            key = (clone, condition)
            if key not in condition_penalties.index:
                continue
            gsea = pd.read_csv(path, sep="\t", usecols=["ID", "factor", "p.adjust"])
            gsea["p.adjust"] = pd.to_numeric(gsea["p.adjust"], errors="coerce")
            significant = gsea[gsea["p.adjust"] < padj_cutoff]
            factor_term_sets = {
                factor: set(group["ID"])
                for factor, group in significant.groupby("factor")
            }
            penalties = condition_penalties.loc[key]
            rows.append(
                {
                    "clone": clone,
                    "condition": condition,
                    "L_l1_strength": penalties["L_l1_strength"],
                    "L_loading_overlap_strength": penalties[
                        "L_loading_overlap_strength"
                    ],
                    "percent_sig_factors": (
                        100.0
                        * significant["factor"].nunique()
                        / float(penalties["dim"])
                    ),
                    "mean_pairwise_gsea_jaccard": mean_pairwise_jaccard(
                        factor_term_sets
                    ),
                }
            )
    if not rows:
        raise ValueError(f"No GSEA result files found under {runs_dir}")
    return pd.DataFrame(rows)


def penalty_positions(values):
    positions = {value: i for i, value in enumerate(PENALTY_VALUES)}
    return np.array([positions[float(value)] for value in values])


def draw_response_panel(ax, df, metric, ylabel):
    for color, l1_value in zip(PENALTY_COLORS, DISPLAY_L1_VALUES):
        sub = df[df["L_l1_strength"] == l1_value]
        grouped = sub.groupby("L_loading_overlap_strength")[metric]
        median = grouped.median().reindex(PENALTY_VALUES)
        lower = grouped.quantile(0.25).reindex(PENALTY_VALUES)
        upper = grouped.quantile(0.75).reindex(PENALTY_VALUES)
        x = penalty_positions(median.index)
        label = rf"$\lambda_{{L1}}={format_penalty(l1_value)}$"
        ax.plot(
            x,
            median,
            marker="o",
            linewidth=1.8,
            markersize=4,
            color=color,
            label=label,
        )
        ax.fill_between(x, lower, upper, color=color, alpha=0.13, linewidth=0)

    ax.set_xticks(range(len(PENALTY_VALUES)))
    ax.set_xticklabels([format_penalty(value) for value in PENALTY_VALUES])
    ax.set_xlabel(r"Loading-overlap penalty $\lambda_{\mathrm{overlap}}$")
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", color="black", alpha=0.18, linewidth=0.6)
    sns.despine(ax=ax)


def plot_tradeoff(df, gsea_df, output, gsea_padj_cutoff):
    fig, axes = plt.subplots(2, 3, figsize=(10.3, 6.2), constrained_layout=True)
    draw_response_panel(
        axes[0, 0],
        df,
        "top_gene_mean_pairwise_unsigned_jaccard",
        "Mean pairwise top-gene Jaccard",
    )
    draw_response_panel(
        axes[0, 1],
        df,
        "mean_offdiag_abs_L_pearson",
        "Mean pairwise L row\nPearson correlation ($|r|$)",
    )
    draw_response_panel(
        axes[0, 2],
        df,
        "log_likelihood_change_per_observation",
        "Change in log likelihood\nper cell–gene observation",
    )
    draw_response_panel(
        axes[1, 0],
        gsea_df,
        "mean_pairwise_gsea_jaccard",
        "Mean pairwise GO-term Jaccard",
    )
    draw_response_panel(
        axes[1, 1],
        gsea_df,
        "percent_sig_factors",
        f"Factors with significant GO terms (%)\n(adjusted p < {gsea_padj_cutoff:g})",
    )
    axes[1, 2].set_visible(False)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        frameon=False,
        title="L1 penalty",
        loc="center left",
        bbox_to_anchor=(1.0, 0.5),
    )
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def clone_label(clone):
    if clone.startswith("quinn"):
        return f"Quinn {clone.removeprefix('quinn')}"
    if clone.startswith("tls"):
        return f"TLS {clone.removeprefix('tls')}"
    return clone


def annotate_points(ax, sub, metric):
    x = np.log10(sub["n_observations"].to_numpy())
    y = sub[metric].to_numpy()
    median_x = float(np.median(x))
    for xi, yi, clone in zip(x, y, sub["clone"]):
        offset = (5, 4) if xi <= median_x else (-5, 4)
        ha = "left" if xi <= median_x else "right"
        ax.annotate(
            clone_label(clone),
            (10**xi, yi),
            xytext=offset,
            textcoords="offset points",
            fontsize=7.5,
            ha=ha,
            va="bottom",
        )


def draw_size_panel(ax, sub, metric, ylabel, letter):
    x = sub["n_observations"].to_numpy()
    y = sub[metric].to_numpy()
    ax.scatter(
        x,
        y,
        s=38,
        color="#4C72B0",
        edgecolor="white",
        linewidth=0.7,
        zorder=3,
    )
    ax.axhline(np.median(y), color="#555555", linestyle="--", linewidth=1)
    annotate_points(ax, sub, metric)
    ax.set_xscale("log")
    ax.set_xlabel("Modeled cell–gene observations")
    ax.set_ylabel(ylabel)
    ax.set_title(letter, loc="left", fontweight="bold")
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.6)
    sns.despine(ax=ax)


def plot_normalization(df, output):
    sub = df[
        np.isclose(df["L_l1_strength"], RECOMMENDED_L1)
        & np.isclose(df["L_loading_overlap_strength"], RECOMMENDED_OVERLAP)
    ].copy()
    if sub.empty:
        raise ValueError(
            "Recommended setting is absent: "
            f"L1={RECOMMENDED_L1:g}, overlap={RECOMMENDED_OVERLAP:g}"
        )
    size_fold = sub["n_observations"].max() / sub["n_observations"].min()

    fig, axes = plt.subplots(1, 3, figsize=(10.3, 3.25), constrained_layout=True)
    draw_size_panel(
        axes[0],
        sub,
        "top_gene_mean_pairwise_signed_jaccard",
        "Mean top-gene overlap\n(signed Jaccard)",
        "A",
    )
    draw_size_panel(
        axes[1],
        sub,
        "mean_offdiag_abs_L_pearson",
        r"Mean loading correlation ($|r|$)",
        "B",
    )
    draw_size_panel(
        axes[2],
        sub,
        "log_likelihood_change_per_observation",
        "Change in log likelihood\nper cell–gene observation",
        "C",
    )
    fig.suptitle(
        rf"$\lambda_{{L1}}={RECOMMENDED_L1:g}$, "
        rf"$\lambda_{{\mathrm{{overlap}}}}={RECOMMENDED_OVERLAP:g}$ "
        f"across a {size_fold:.0f}-fold matrix-size range",
        fontsize=10,
    )
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    try:
        df = load_data(
            [Path(path) for path in args.selection_summaries],
            Path(args.runs_dir),
        )
        if args.figure == "tradeoff":
            gsea_df = load_gsea_data(
                df,
                Path(args.runs_dir),
                args.gsea_padj_cutoff,
            )
            plot_tradeoff(df, gsea_df, output, args.gsea_padj_cutoff)
        else:
            plot_normalization(df, output)
    except (FileNotFoundError, KeyError, ValueError) as exc:
        sys.stderr.write(f"ERROR: {exc}\n")
        return 1
    print(f"Saved figure to: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
