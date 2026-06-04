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
    return parser.parse_args()


def parse_condition(condition):
    prefixes = ["Fl2", "Ll", "Lc", "Fc", "Fo", "Fa", "V"]
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
    return str(value).lower() in {"0", "0.0", "false", "no"}


def is_no_penalty_condition(condition):
    parts = parse_condition(condition)
    return (
        is_zero(parts.get("Ll"))
        and is_zero(parts.get("Lc"))
        and is_zero(parts.get("Fc"))
        and is_zero(parts.get("Fo"))
        and is_false(parts.get("Fa"))
        and is_false(parts.get("V"))
    )


def condition_short_label(condition):
    parts = parse_condition(condition)
    keys = ["Ll", "Lc", "Fc", "Fo", "Fa", "V"]
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
    required_cols = {
        "clone",
        "condition",
        "mean_offdiag_abs_L_correlation",
        "mean_offdiag_L_correlation",
    }
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        sys.stderr.write(f"Missing required columns: {sorted(missing_cols)}\n")
        return 1

    numeric_cols = [
        "mean_offdiag_abs_L_correlation",
        "mean_offdiag_L_correlation",
    ]
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

    df["abs_signed_overlap"] = df["mean_offdiag_L_correlation"].abs()
    df["is_no_penalty"] = df["condition"].map(is_no_penalty_condition)
    clone_order = sorted(df["clone"].unique())

    nrows = len(clone_order)
    fig, axes = plt.subplots(
        nrows,
        2,
        figsize=(9.4, 3.3 * nrows),
        sharex="col",
        sharey=False,
        squeeze=False,
    )
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

    for row_idx, clone in enumerate(clone_order):
        clone_df = df[df["clone"] == clone]
        default_df = clone_df[clone_df["is_no_penalty"]]
        for col_idx, (x_col, xlabel) in enumerate(plot_specs):
            ax = axes[row_idx, col_idx]
            ax.scatter(
                clone_df[x_col],
                clone_df["total_log_likelihood"],
                s=50,
                color="#56B4E9",
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
                    color="#D55E00",
                    edgecolor="white",
                    linewidth=0.7,
                    zorder=5,
                )

            if x_col == "mean_offdiag_L_correlation":
                ax.axvline(0, color="0.2", linewidth=0.8, linestyle="--", alpha=0.7)

            title_suffix = "absolute overlap" if col_idx == 0 else "signed overlap"
            ax.set_title(f"{clone}: {title_suffix}")
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
            markerfacecolor="#56B4E9",
            markeredgecolor="white",
            markersize=7,
            label="Various penalties applied",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#D55E00",
            markeredgecolor="white",
            markersize=10,
            label="Default no-penalty model",
        ),
    ]
    fig.legend(
        handles=legend_handles,
        loc="center left",
        frameon=False,
        bbox_to_anchor=(0.92, 0.5),
    )

    if not df["is_no_penalty"].any():
        sys.stderr.write("Warning: no all-penalties-off baseline condition found.\n")

    fig.subplots_adjust(hspace=0.48, wspace=0.32)
    fig.savefig(args.output_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure to: {args.output_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
