#!/usr/bin/env python3

import argparse
import os
import sys

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import seaborn as sns

sns.set_theme(style="ticks", context="paper")

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.labelsize": 10,
        "axes.titlesize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot penalty sweep overlap against unpenalized log-likelihood."
    )
    parser.add_argument("summary_tsv")
    parser.add_argument("output_pdf")
    parser.add_argument("--metric-col", default=None)
    parser.add_argument("--metric-label", default=None)
    parser.add_argument("--clone", default=None)
    return parser.parse_args()


def parse_condition(condition):
    prefixes = ["Fl2", "Ll", "Lc", "Lo", "Fc", "Fo", "Fa", "K", "V"]
    parts = {}
    for token in condition.split("_"):
        for prefix in prefixes:
            if token.startswith(prefix):
                parts[prefix] = token[len(prefix):]
                break
    return parts


def is_zero(value):
    try:
        return float(value) == 0.0
    except (TypeError, ValueError):
        return False


def is_false(value):
    return str(value).lower() in {"0", "0.0", "false", "no", "none"}


def is_no_penalty_condition(condition):
    parts = parse_condition(condition)
    return (
        is_zero(parts.get("Ll"))
        and is_zero(parts.get("Lc"))
        and is_zero(parts.get("Lo", "0"))
        and is_zero(parts.get("Fc"))
        and is_zero(parts.get("Fo"))
        and is_false(parts.get("Fa"))
        and is_false(parts.get("V"))
    )


BASELINE_ORDER_KEYS = ["K", "Ll", "Lc", "Lo", "Fc", "Fo", "Fa", "Fl2", "V"]


def baseline_sort_value(value):
    if str(value).lower() in {"none", "0", "0.0", "false", "no", ""}:
        return 0.0
    if str(value).lower() == "tree":
        return 1.0
    if str(value).lower() == "iid":
        return 2.0
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("inf")


def mark_simplest_baseline(df):
    df = df.copy()
    df["is_no_penalty"] = False
    parts = df["condition"].map(parse_condition)
    for key in BASELINE_ORDER_KEYS:
        df[key] = parts.map(lambda p, k=key: p.get(k, ""))

    for _, clone_df in df.groupby("clone", sort=True):
        baseline_idx = min(
            clone_df.index,
            key=lambda idx: tuple(
                baseline_sort_value(df.at[idx, key]) for key in BASELINE_ORDER_KEYS
            )
            + (str(df.at[idx, "condition"]),),
        )
        df.at[baseline_idx, "is_no_penalty"] = True

    return df


def condition_short_label(condition):
    parts = parse_condition(condition)
    keys = ["K", "Ll", "Lc", "Lo", "Fc", "Fo", "Fa", "V"]
    return " ".join(f"{key}={parts[key]}" for key in keys if key in parts)


def add_likelihood_component(df, component):
    ll_candidates = [
        f"{component}_ll",
        f"{component}_log_likelihood",
        f"{component}_likelihood",
    ]
    for col in ll_candidates:
        if col in df.columns:
            return pd.to_numeric(df[col], errors="coerce")

    negll_col = f"{component}_negll"
    if negll_col in df.columns:
        return -pd.to_numeric(df[negll_col], errors="coerce")

    sys.stderr.write(
        f"Missing {component} likelihood column. Expected one of "
        f"{ll_candidates} or {negll_col}.\n"
    )
    raise SystemExit(1)


def main():
    args = parse_args()

    df = pd.read_csv(args.summary_tsv, sep="\t")
    if args.clone:
        df = df[df["clone"] == args.clone].copy()
        if df.empty:
            sys.stderr.write(f"No rows found for clone '{args.clone}'.\n")
            return 1

    if args.metric_col:
        plot_specs = [(args.metric_col, args.metric_label or args.metric_col)]
    else:
        plot_specs = [
            (
                "mean_offdiag_abs_L_correlation",
                "Mean off-diagonal |L correlation|",
            ),
            (
                "mean_offdiag_L_correlation",
                "Mean off-diagonal L correlation",
            ),
        ]

    required_cols = {"clone", "condition"} | {col for col, _ in plot_specs}
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        sys.stderr.write(f"Missing required columns: {sorted(missing_cols)}\n")
        return 1

    numeric_cols = [col for col, _ in plot_specs]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["brownian_log_likelihood"] = add_likelihood_component(df, "brownian")
    df["observation_log_likelihood"] = add_likelihood_component(df, "observation")
    df["total_log_likelihood"] = (
        df["brownian_log_likelihood"] + df["observation_log_likelihood"]
    )
    numeric_cols.extend(
        ["brownian_log_likelihood", "observation_log_likelihood", "total_log_likelihood"]
    )

    df = df.dropna(subset=numeric_cols).copy()
    if df.empty:
        sys.stderr.write("No complete rows to plot.\n")
        return 1

    df = mark_simplest_baseline(df)
    clone_order = sorted(df["clone"].unique())

    nrows = len(clone_order)
    fig, axes = plt.subplots(
        nrows,
        len(plot_specs),
        figsize=(4.9 * len(plot_specs), 3.3 * nrows),
        sharex="col",
        sharey=False,
        squeeze=False,
    )

    for row_idx, clone in enumerate(clone_order):
        clone_df = df[df["clone"] == clone]
        default_df = clone_df[clone_df["is_no_penalty"]]
        for col_idx, (x_col, xlabel) in enumerate(plot_specs):
            ax = axes[row_idx, col_idx]
            ax.scatter(
                clone_df[x_col],
                clone_df["total_log_likelihood"],
                s=50,
                color="#F0E442",
                edgecolor="white",
                linewidth=0.5,
                alpha=0.9,
            )

            if not default_df.empty:
                default_row = default_df.iloc[0]
                ax.scatter(
                    [default_row[x_col]],
                    [default_row["total_log_likelihood"]],
                    s=120,
                    marker="o",
                    color="#E69F00",
                    edgecolor="white",
                    linewidth=0.7,
                    zorder=5,
                )

            if x_col == "mean_offdiag_L_correlation":
                ax.axvline(0, color="0.2", linewidth=0.8, linestyle="--", alpha=0.7)

            ax.set_title(clone)
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Brownian + observation log-likelihood")
            ax.tick_params(
                axis="both",
                which="both",
                labelbottom=True,
                labelleft=True,
            )
            ax.grid(True, axis="both", linewidth=0.5, alpha=0.35)
            ax.set_axisbelow(True)
            sns.despine(ax=ax)

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#F0E442",
            markeredgecolor="white",
            markersize=7,
            label="Altered model",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#E69F00",
            markeredgecolor="white",
            markersize=10,
            label="Default model",
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="center left",
        frameon=False,
        bbox_to_anchor=(0.92, 0.5),
    )

    if not df["is_no_penalty"].any():
        sys.stderr.write("Warning: no simplest baseline condition found.\n")

    fig.subplots_adjust(hspace=0.48, wspace=0.32)
    fig.savefig(args.output_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure to: {args.output_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
