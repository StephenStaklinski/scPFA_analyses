#!/usr/bin/env python3

import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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

input_tsv = sys.argv[1]
output_pdf = sys.argv[2]

df = pd.read_csv(input_tsv, sep="\t")

fit_cols = [c for c in df.columns if c.endswith("_fit")]
pairs = []

for fit_col in fit_cols:
    base = fit_col[:-4]
    sim_col = f"{base}_simulated"
    if sim_col in df.columns:
        pairs.append((base, fit_col, sim_col))

long_rows = []

for _, row in df.iterrows():
    condition = row["condition"] if "condition" in df.columns else "all"
    sim_num = row["sim_num"] if "sim_num" in df.columns else np.nan

    for base, fit_col, sim_col in pairs:
        long_rows.append(
            {
                "condition": str(condition),
                "sim_num": sim_num,
                "metric": base,
                "source": "fit",
                "value": row[fit_col],
            }
        )
        long_rows.append(
            {
                "condition": str(condition),
                "sim_num": sim_num,
                "metric": base,
                "source": "simulated",
                "value": row[sim_col],
            }
        )

long_df = pd.DataFrame(long_rows)

long_df["source"] = pd.Categorical(
    long_df["source"],
    categories=["simulated", "fit"],
    ordered=True
)

metrics = list(dict.fromkeys(long_df["metric"]))
n = len(metrics)

ncols = int(math.ceil(math.sqrt(n)))
nrows = int(math.ceil(n / ncols))

fig, axes = plt.subplots(
    nrows,
    ncols,
    figsize=(3.0 * ncols, 3.0 * nrows),
    sharex=False
)
axes = np.array(axes).reshape(-1)

for i, metric in enumerate(metrics):
    ax = axes[i]
    sub = long_df[long_df["metric"] == metric].copy()

    sns.boxplot(
        data=sub,
        x="condition",
        y="value",
        hue="source",
        ax=ax,
        dodge=True,
        fliersize=2,
        linewidth=1,
    )

    sns.stripplot(
        data=sub,
        x="condition",
        y="value",
        hue="source",
        ax=ax,
        dodge=True,
        alpha=0.7,
        size=3,
        legend=False,
    )

    y = pd.to_numeric(sub["value"], errors="coerce").dropna().to_numpy()
    y_positive = y[y > 0]

    if len(y_positive) == len(y) and len(y_positive) > 0:
        if y_positive.max() / y_positive.min() > 1000:
            ax.set_yscale("log")

    ax.set_title(metric, pad=6)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", rotation=0)

    if i != 0 and ax.get_legend() is not None:
        ax.get_legend().remove()

    sns.despine(ax=ax)
    ax.grid(True, axis="y", alpha=1.0)
    ax.grid(False, axis="x")

for j in range(i + 1, len(axes)):
    axes[j].axis("off")

if len(metrics) > 0 and axes[0].get_legend() is not None:
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles[:2], labels[:2], loc="upper center", ncol=2, frameon=False)
    axes[0].get_legend().remove()

fig.supxlabel("")
fig.supylabel("")

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(output_pdf, bbox_inches="tight")
plt.close()
