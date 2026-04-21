#!/usr/bin/env python3

import sys
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

# Read TSV
df = pd.read_csv(input_tsv, sep="\t", index_col=0)

# Plot style
sns.set_theme(style="white")

# Figure size scales with matrix shape
nrows, ncols = df.shape
if nrows > ncols:
    figsize = (4, 6)
else:
    figsize = (6, 4)

plt.figure(figsize=figsize)

# Heatmap
sns.heatmap(
    df,
    cmap="coolwarm",
    linewidths=0,
    cbar=True,
    # xticklabels=False,
    # yticklabels=False
)

plt.tight_layout()
plt.savefig(output_pdf)
plt.close()