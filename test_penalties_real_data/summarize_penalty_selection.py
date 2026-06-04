#!/usr/bin/env python3

import argparse
import math
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
from matplotlib.gridspec import GridSpecFromSubplotSpec
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
        description="Select and plot penalty settings with good fit/overlap tradeoffs."
    )
    parser.add_argument("summary_tsv")
    parser.add_argument("selected_tsv")
    parser.add_argument("output_pdf")
    parser.add_argument("--worst-tsv", default="worst_penalty_runs.tsv")
    parser.add_argument(
        "--top-frac",
        type=float,
        default=0.25,
        help="Fraction of runs to keep per clone by combined rank.",
    )
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


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def is_zero(value):
    value = as_float(value)
    return not math.isnan(value) and value == 0.0


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


def add_likelihood_component(df, component):
    for col in [
        f"{component}_ll",
        f"{component}_log_likelihood",
        f"{component}_likelihood",
    ]:
        if col in df.columns:
            return pd.to_numeric(df[col], errors="coerce")

    negll_col = f"{component}_negll"
    if negll_col in df.columns:
        return -pd.to_numeric(df[negll_col], errors="coerce")

    sys.stderr.write(f"Missing likelihood column for {component}.\n")
    raise SystemExit(1)


def load_summary(path):
    df = pd.read_csv(path, sep="\t")
    required_cols = {
        "clone",
        "condition",
        "mean_offdiag_abs_L_correlation",
        "mean_offdiag_L_correlation",
    }
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        sys.stderr.write(f"Missing required columns: {sorted(missing_cols)}\n")
        raise SystemExit(1)

    for col in ["mean_offdiag_abs_L_correlation", "mean_offdiag_L_correlation"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["brownian_log_likelihood"] = add_likelihood_component(df, "brownian")
    df["observation_log_likelihood"] = add_likelihood_component(df, "observation")
    df["total_log_likelihood"] = (
        df["brownian_log_likelihood"] + df["observation_log_likelihood"]
    )
    df["abs_signed_overlap"] = df["mean_offdiag_L_correlation"].abs()
    df["is_no_penalty"] = df["condition"].map(is_no_penalty_condition)

    parts = df["condition"].map(parse_condition)
    for key in ["Ll", "Lc", "Fc", "Fo", "Fa", "V"]:
        df[key] = parts.map(lambda p, k=key: p.get(k, ""))
    return df.dropna(
        subset=[
            "mean_offdiag_abs_L_correlation",
            "mean_offdiag_L_correlation",
            "total_log_likelihood",
        ]
    ).copy()


def add_baseline_deltas(df):
    out = []
    for clone, clone_df in df.groupby("clone", sort=True):
        default_df = clone_df[clone_df["is_no_penalty"]]
        if default_df.empty:
            sys.stderr.write(f"Warning: no no-penalty baseline for {clone}; skipping.\n")
            continue

        default = default_df.iloc[0]
        clone_df = clone_df.copy()
        clone_df["default_total_log_likelihood"] = default["total_log_likelihood"]
        clone_df["default_abs_overlap"] = default["mean_offdiag_abs_L_correlation"]
        clone_df["default_abs_signed_overlap"] = default["abs_signed_overlap"]
        clone_df["log_likelihood_delta"] = (
            clone_df["total_log_likelihood"] - default["total_log_likelihood"]
        )
        clone_df["log_likelihood_loss"] = (
            default["total_log_likelihood"] - clone_df["total_log_likelihood"]
        ).clip(lower=0)
        clone_df["abs_overlap_improvement"] = (
            default["mean_offdiag_abs_L_correlation"]
            - clone_df["mean_offdiag_abs_L_correlation"]
        )
        clone_df["signed_overlap_improvement"] = (
            default["abs_signed_overlap"] - clone_df["abs_signed_overlap"]
        )
        out.append(clone_df)

    if not out:
        sys.stderr.write("No clone had a no-penalty baseline.\n")
        raise SystemExit(1)
    return pd.concat(out, ignore_index=True)


def score_tradeoffs(df):
    scored = []
    for clone, clone_df in df.groupby("clone", sort=True):
        clone_df = clone_df.copy()
        clone_df["log_likelihood_loss_rank"] = clone_df["log_likelihood_loss"].rank(
            method="min", ascending=True
        )
        clone_df["abs_overlap_improvement_rank"] = clone_df[
            "abs_overlap_improvement"
        ].rank(method="min", ascending=False)
        clone_df["tradeoff_rank"] = (
            clone_df["log_likelihood_loss_rank"]
            + clone_df["abs_overlap_improvement_rank"]
        )
        scored.append(clone_df)
    return pd.concat(scored, ignore_index=True)


def select_top_runs(df, top_frac):
    selected = []
    for clone, clone_df in df.groupby("clone", sort=True):
        n_keep = max(1, math.ceil(len(clone_df) * top_frac))
        clone_df = clone_df.sort_values(
            [
                "tradeoff_rank",
                "log_likelihood_loss_rank",
                "abs_overlap_improvement_rank",
                "condition",
            ]
        )
        selected.append(clone_df.head(n_keep))
    return pd.concat(selected, ignore_index=True)


def select_worst_runs(df, top_frac):
    selected = []
    for clone, clone_df in df.groupby("clone", sort=True):
        n_keep = max(1, math.ceil(len(clone_df) * top_frac))
        clone_df = clone_df.sort_values(
            [
                "tradeoff_rank",
                "log_likelihood_loss_rank",
                "abs_overlap_improvement_rank",
                "condition",
            ],
            ascending=[False, False, False, True],
        )
        selected.append(clone_df.head(n_keep))
    return pd.concat(selected, ignore_index=True)


def penalty_active(row, key):
    if key in {"Fa", "V"}:
        return not is_false(row[key])
    return not is_zero(row[key])


PENALTY_LABELS = {
    "Ll": "L L1",
    "Lc": "L corr.",
    "Fc": "F corr.",
    "Fo": "F orth.",
    "Fa": "Final absorbing",
    "V": "Varimax",
}


def active_penalty_keys(row):
    return tuple(key for key in PENALTY_LABELS if penalty_active(row, key))


def make_combination_data(selected_df):
    rows = []
    for _, row in selected_df.iterrows():
        keys = active_penalty_keys(row)
        rows.append({"active_keys": keys, "combination": "+".join(keys) or "none"})
    combo_df = pd.DataFrame(rows)
    if combo_df.empty:
        return pd.DataFrame(columns=["active_keys", "combination", "n"])

    counts = (
        combo_df.groupby(["active_keys", "combination"], sort=False)
        .size()
        .reset_index(name="n")
        .sort_values(["n", "combination"], ascending=[False, True])
        .reset_index(drop=True)
    )
    return counts


def plot_combination_panel(subplot_spec, selected_clone, title, color):
    combo_df = make_combination_data(selected_clone)
    inner = GridSpecFromSubplotSpec(
        2,
        2,
        subplot_spec=subplot_spec,
        height_ratios=[1.25, 1.55],
        width_ratios=[1.0, 0.55],
        hspace=0.05,
        wspace=0.08,
    )
    ax_bar = plt.subplot(inner[0, 0])
    ax_blank = plt.subplot(inner[0, 1])
    ax_matrix = plt.subplot(inner[1, 0], sharex=ax_bar)
    ax_marginal = plt.subplot(inner[1, 1], sharey=ax_matrix)
    ax_blank.axis("off")

    if combo_df.empty:
        ax_bar.set_title(title)
        ax_bar.text(0.5, 0.5, "No selected runs", ha="center", va="center")
        ax_bar.axis("off")
        ax_marginal.axis("off")
        ax_matrix.axis("off")
        return ax_bar, ax_marginal, ax_matrix

    x = list(range(len(combo_df)))
    ax_bar.bar(x, combo_df["n"], color=color, width=0.72)
    ax_bar.set_title(title)
    ax_bar.set_ylabel("Count")
    ax_bar.tick_params(axis="x", labelbottom=False, bottom=False)
    ax_bar.grid(True, axis="y", linewidth=0.5, alpha=0.35)
    ax_bar.set_axisbelow(True)
    sns.despine(ax=ax_bar, bottom=True)

    penalty_keys = list(PENALTY_LABELS.keys())
    y_positions = list(range(len(penalty_keys)))
    row_colors = ["#F5F5F5", "white"]
    for y in y_positions:
        ax_matrix.axhspan(y - 0.5, y + 0.5, color=row_colors[y % 2], zorder=0)
        ax_marginal.axhspan(y - 0.5, y + 0.5, color=row_colors[y % 2], zorder=0)

    marginal_counts = [
        selected_clone.apply(lambda row, k=key: penalty_active(row, k), axis=1).sum()
        for key in penalty_keys
    ]
    ax_marginal.barh(y_positions, marginal_counts, color=color, height=0.58)
    max_marginal = max(marginal_counts) if marginal_counts else 0
    ax_marginal.set_xlim(0, max(max_marginal * 1.25, 1))
    ax_marginal.set_yticks(y_positions)
    ax_marginal.set_yticklabels([])
    ax_marginal.tick_params(axis="y", left=False, labelleft=False)
    ax_marginal.invert_yaxis()
    ax_marginal.set_xlabel("Count")
    ax_marginal.grid(True, axis="x", linewidth=0.5, alpha=0.35)
    ax_marginal.set_axisbelow(True)
    sns.despine(ax=ax_marginal, left=True)

    for xpos, active_keys in zip(x, combo_df["active_keys"]):
        active_y = [i for i, key in enumerate(penalty_keys) if key in active_keys]
        if active_y:
            ax_matrix.plot(
                [xpos, xpos],
                [min(active_y), max(active_y)],
                color=color,
                linewidth=1.0,
                zorder=2,
            )
        for y, key in enumerate(penalty_keys):
            is_active = key in active_keys
            ax_matrix.scatter(
                xpos,
                y,
                s=46 if is_active else 38,
                color=color if is_active else "0.88",
                edgecolor=color if is_active else "0.85",
                linewidth=0.4,
                zorder=3,
            )

    ax_matrix.set_yticks(y_positions)
    ax_matrix.set_yticklabels([PENALTY_LABELS[key] for key in penalty_keys])
    ax_matrix.invert_yaxis()
    ax_matrix.set_xlabel("Penalty combination")
    ax_matrix.set_xlim(-0.6, len(combo_df) - 0.4)
    ax_matrix.set_xticks(x)
    ax_matrix.set_xticklabels([])
    ax_matrix.tick_params(axis="x", bottom=False)
    ax_matrix.tick_params(axis="y", length=0)
    sns.despine(ax=ax_matrix, bottom=True, left=True)
    return ax_bar, ax_marginal, ax_matrix


def make_plot(scored_df, selected_df, worst_df, output_pdf):
    clone_order = sorted(scored_df["clone"].unique())
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#F0E442",
            markeredgecolor="white",
            markersize=7,
            label="All runs",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#56B4E9",
            markeredgecolor="white",
            markersize=8,
            label="Best quartile",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#009E73",
            markeredgecolor="white",
            markersize=8,
            label="Worst quartile",
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
    fig = plt.figure(figsize=(17.0, 3.8 * len(clone_order)))
    outer = fig.add_gridspec(
        len(clone_order),
        2,
        width_ratios=[1.0, 3.0],
        hspace=0.48,
        wspace=0.62,
    )

    for row_idx, clone in enumerate(clone_order):
        clone_df = scored_df[scored_df["clone"] == clone]
        selected_clone = selected_df[selected_df["clone"] == clone]
        worst_clone = worst_df[worst_df["clone"] == clone]

        ax = fig.add_subplot(outer[row_idx, 0])
        ax.scatter(
            clone_df["abs_overlap_improvement"],
            clone_df["log_likelihood_delta"],
            s=42,
            color="#F0E442",
            edgecolor="white",
            linewidth=0.4,
            alpha=0.85,
        )
        ax.scatter(
            selected_clone["abs_overlap_improvement"],
            selected_clone["log_likelihood_delta"],
            s=56,
            color="#56B4E9",
            edgecolor="white",
            linewidth=0.5,
            alpha=0.95,
        )
        ax.scatter(
            worst_clone["abs_overlap_improvement"],
            worst_clone["log_likelihood_delta"],
            s=56,
            color="#009E73",
            edgecolor="white",
            linewidth=0.5,
            alpha=0.95,
        )
        default_df = clone_df[clone_df["is_no_penalty"]]
        if not default_df.empty:
            default = default_df.iloc[0]
            ax.scatter(
                [default["abs_overlap_improvement"]],
                [default["log_likelihood_delta"]],
                s=110,
                color="#E69F00",
                edgecolor="white",
                linewidth=0.7,
                zorder=5,
            )

        ax.axhline(0, color="0.2", linewidth=0.8, linestyle=":", alpha=0.8)
        ax.axvline(0, color="0.2", linewidth=0.8, linestyle=":", alpha=0.8)
        ax.set_title(clone)
        ax.set_xlabel("Reduction in mean |L correlation| vs. default")
        ax.set_ylabel("Log-likelihood change vs. default")
        ax.grid(True, axis="both", linewidth=0.5, alpha=0.35)
        ax.set_axisbelow(True)
        ax.legend(
            handles=legend_handles,
            loc="center left",
            frameon=False,
            bbox_to_anchor=(1.02, 0.5),
        )
        sns.despine(ax=ax)

        upset_cols = GridSpecFromSubplotSpec(
            1,
            2,
            subplot_spec=outer[row_idx, 1],
            wspace=0.38,
        )

        plot_combination_panel(
            upset_cols[0],
            selected_clone,
            clone,
            "#56B4E9",
        )
        plot_combination_panel(
            upset_cols[1],
            worst_clone,
            clone,
            "#009E73",
        )

    fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    if args.top_frac <= 0 or args.top_frac > 1:
        sys.stderr.write("--top-frac must be in (0, 1].\n")
        return 1

    df = load_summary(args.summary_tsv)
    scored_df = add_baseline_deltas(df)
    scored_df = score_tradeoffs(scored_df)
    selected_df = select_top_runs(scored_df, args.top_frac)
    worst_df = select_worst_runs(scored_df, args.top_frac)

    output_cols = [
        "clone",
        "condition",
        "total_log_likelihood",
        "default_total_log_likelihood",
        "log_likelihood_delta",
        "log_likelihood_loss",
        "mean_offdiag_abs_L_correlation",
        "default_abs_overlap",
        "abs_overlap_improvement",
        "mean_offdiag_L_correlation",
        "signed_overlap_improvement",
        "tradeoff_rank",
        "Ll",
        "Lc",
        "Fc",
        "Fo",
        "Fa",
        "V",
    ]
    selected_df[output_cols].to_csv(args.selected_tsv, sep="\t", index=False)
    worst_df[output_cols].to_csv(args.worst_tsv, sep="\t", index=False)
    make_plot(scored_df, selected_df, worst_df, args.output_pdf)
    print(f"Saved selected runs to: {args.selected_tsv}")
    print(f"Saved worst runs to: {args.worst_tsv}")
    print(f"Saved figure to: {args.output_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
