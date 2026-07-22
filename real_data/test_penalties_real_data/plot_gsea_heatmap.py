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
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import seaborn as sns


sns.set_theme(style="white", context="paper")

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

CONDITION_PREFIXES = ["Ll", "Lo", "Rm", "Fa", "K", "T"]

RM_LABELS = {
    0: "Ribo/mito genes included",
    1: "Ribo/mito genes removed",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Heatmap of GSEA biological enrichment across penalty sweep."
    )
    parser.add_argument("runs_dir")
    parser.add_argument("output_pdf")
    parser.add_argument("--clone", required=True)
    parser.add_argument("--padj-cutoff", type=float, default=0.05)
    parser.add_argument(
        "--no-annot",
        action="store_true",
        help="Do not print metric values inside heatmap cells.",
    )
    return parser.parse_args()


def parse_condition(condition):
    parts = {}
    for token in condition.split("_"):
        for prefix in CONDITION_PREFIXES:
            if token.startswith(prefix):
                parts[prefix] = token[len(prefix):]
                break
    return parts


def format_penalty(value):
    value = float(value)
    return "0" if value == 0 else f"{value:g}"


def mean_pairwise_jaccard(factor_term_sets):
    """Mean pairwise Jaccard of significant GO term sets.

    Only considers pairs where both factors have ≥1 significant term.
    Returns NaN if fewer than 2 enriched factors exist.
    """
    enriched = {f: terms for f, terms in factor_term_sets.items() if terms}
    factors = sorted(enriched)
    values = []
    for i, fa in enumerate(factors):
        for fb in factors[i + 1:]:
            a, b = enriched[fa], enriched[fb]
            union = a | b
            values.append(len(a & b) / len(union) if union else float("nan"))
    valid = [v for v in values if not math.isnan(v)]
    return sum(valid) / len(valid) if valid else float("nan")


def load_gsea_metrics(runs_dir, clone, padj_cutoff):
    rows = []
    for path in sorted(runs_dir.glob(f"{clone}/*.L.gsea.tsv")):
        condition = path.name.removesuffix(".L.gsea.tsv")
        parts = parse_condition(condition)
        try:
            ll = float(parts["Ll"])
            lo = float(parts["Lo"])
            rm = int(parts.get("Rm", 0))
        except (KeyError, ValueError):
            continue

        try:
            df = pd.read_csv(path, sep="\t", usecols=["ID", "factor", "p.adjust"])
        except Exception as exc:
            sys.stderr.write(f"Skipping {path}: {exc}\n")
            continue

        sig = df[df["p.adjust"] < padj_cutoff]
        n_sig_factors = sig["factor"].nunique()

        factor_term_sets = {
            factor: set(grp["ID"])
            for factor, grp in sig.groupby("factor")
        }
        all_factors = df["factor"].unique()
        for f in all_factors:
            if f not in factor_term_sets:
                factor_term_sets[f] = set()

        jaccard = mean_pairwise_jaccard(factor_term_sets)

        rows.append(
            {
                "condition": condition,
                "L_l1_strength": ll,
                "L_loading_overlap_strength": lo,
                "remove_ribo_mito": rm,
                "n_sig_factors": n_sig_factors,
                "mean_pairwise_gsea_jaccard": jaccard,
            }
        )

    return pd.DataFrame(rows)


def make_pivot(df, rm_value, value_col):
    sub = df[df["remove_ribo_mito"] == rm_value]
    pivot = sub.pivot_table(
        index="L_loading_overlap_strength",
        columns="L_l1_strength",
        values=value_col,
        aggfunc="mean",
    ).sort_index(ascending=False)
    return pivot.reindex(sorted(pivot.columns), axis=1)


def draw_heatmap(
    ax,
    pivot,
    cmap,
    vmin,
    vmax,
    fmt,
    cbar_label,
    title=None,
    nan_label="—",
    show_annotations=True,
):
    # Build a string annotation array to label NaN cells cleanly
    annot = pivot.copy().astype(object)
    for r in annot.index:
        for c in annot.columns:
            v = pivot.loc[r, c]
            if pd.isna(v):
                annot.loc[r, c] = nan_label
            else:
                annot.loc[r, c] = format(v, fmt)

    cmap_obj = plt.get_cmap(cmap) if isinstance(cmap, str) else cmap
    cmap_obj = cmap_obj.copy()
    cmap_obj.set_bad(color="#DDDDDD")

    sns.heatmap(
        pivot,
        ax=ax,
        cmap=cmap_obj,
        vmin=vmin,
        vmax=vmax,
        annot=annot if show_annotations else False,
        annot_kws={"fontsize": 10},
        fmt="",
        linewidths=0.6,
        linecolor="white",
        cbar_kws={"label": cbar_label},
    )
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold", pad=8)
    ax.set_xlabel("L1 penalty")
    ax.set_ylabel("L-overlap penalty")
    ax.set_xticklabels(
        [format_penalty(x.get_text()) for x in ax.get_xticklabels()]
    )
    ax.set_yticklabels(
        [format_penalty(y.get_text()) for y in ax.get_yticklabels()], rotation=0
    )


def main():
    args = parse_args()
    runs_dir = Path(args.runs_dir)

    df = load_gsea_metrics(runs_dir, args.clone, args.padj_cutoff)
    if df.empty:
        sys.stderr.write(f"No GSEA files found for clone {args.clone} under {runs_dir}\n")
        return 1

    rm_values_present = [rm for rm in [0, 1] if rm in df["remove_ribo_mito"].values]

    try:
        n_factors_total = int(parse_condition(df["condition"].iloc[0]).get("K", 5))
    except (ValueError, IndexError):
        n_factors_total = 5

    cmap_enrichment = LinearSegmentedColormap.from_list(
        "gsea_enrichment",
        ["#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"],
    )
    cmap_jaccard = LinearSegmentedColormap.from_list(
        "gsea_jaccard",
        ["#08306B", "#6BAED6", "#F7FBFF", "#FDD0A2", "#99000D"],
    )

    jaccard_max = df["mean_pairwise_gsea_jaccard"].max()
    jaccard_vmax = jaccard_max if not pd.isna(jaccard_max) else 1.0

    n_cols = len(rm_values_present)
    fig, axes = plt.subplots(
        2, n_cols,
        figsize=(7.0 * n_cols, 11.2),
        constrained_layout=True,
    )
    # Normalise axes to always be 2-D indexable
    if n_cols == 1:
        axes = axes.reshape(2, 1)

    for col, rm in enumerate(rm_values_present):
        title = RM_LABELS[rm]

        pivot_n = make_pivot(df, rm, "n_sig_factors")
        draw_heatmap(
            axes[0, col],
            pivot_n,
            cmap=cmap_enrichment,
            vmin=0,
            vmax=n_factors_total,
            fmt=".0f",
            cbar_label=f"Factors with\nsignificant GO term (of {n_factors_total})",
            title=title,
            show_annotations=not args.no_annot,
        )

        pivot_j = make_pivot(df, rm, "mean_pairwise_gsea_jaccard")
        draw_heatmap(
            axes[1, col],
            pivot_j,
            cmap=cmap_jaccard,
            vmin=0.0,
            vmax=jaccard_vmax,
            fmt=".2f",
            cbar_label="Mean pairwise Jaccard\nsignificant GO terms",
            title=None,
            show_annotations=not args.no_annot,
        )

    # Row labels on the left of the figure
    row_labels = [
        "Enriched factors",
        "GO term redundancy\nbetween factors",
    ]
    for row, label in enumerate(row_labels):
        fig.text(
            -0.01,
            1 - (row + 0.5) / 2,
            label,
            va="center",
            ha="right",
            fontsize=11,
            fontweight="bold",
            rotation=90,
        )

    fig.savefig(args.output_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure to: {args.output_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
