#!/usr/bin/env python3

import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math

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

log_file = sys.argv[1]
n_cols = int(sys.argv[2])
skip_steps = int(sys.argv[3])
outprefix = sys.argv[4]

df = pd.read_csv(log_file, sep="\t", index_col=0, dtype=float, comment="#")
x = df.index.values

# skip steps
if skip_steps > 0:
    mask = x >= skip_steps
    x = x[mask]
    df = df[mask]

# ticks
tick_idx = np.linspace(0, len(x) - 1, 5, dtype=int)
tick_pos = x[tick_idx]

if n_cols == -1:
    n_cols = len(df.columns)

cols = df.columns[:n_cols]
n = len(cols)

# choose grid shape (roughly square)
ncols = int(math.ceil(math.sqrt(n)))
nrows = int(math.ceil(n / ncols))

fig, axes = plt.subplots(
    nrows,
    ncols,
    figsize=(4.2 * ncols, 2.8 * nrows),
    sharex=False
)
axes = np.array(axes).reshape(-1)

for i, col in enumerate(cols):
    ax = axes[i]
    y = df[col].to_numpy()

    sns.lineplot(
        x=x,
        y=y,
        ax=ax,
        linewidth=1.8
    )

    # conditional log scale
    y_positive = y[y > 0]
    if len(y_positive) == len(y) and len(y_positive) > 0:
        if y_positive.max() / y_positive.min() > 1000:
            ax.set_yscale("log")

    ax.set_title(col, pad=6)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels([str(int(v)) for v in tick_pos], rotation=0)

    # cleaner axes
    sns.despine(ax=ax)
    ax.grid(True, axis="y", alpha=1.0)
    ax.grid(True, axis="x", alpha=1.0)
    ax.margins(x=0.03)
    ax.set_xlabel("")
    ax.set_ylabel("")

# hide unused axes
for j in range(i + 1, len(axes)):
    axes[j].axis("off")

fig.supxlabel("Optimization step")

plt.tight_layout()
plt.savefig(f"{outprefix}.pdf", bbox_inches="tight")
plt.close()