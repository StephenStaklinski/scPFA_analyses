#!/usr/bin/env python3

import argparse
import os
import re
import sys
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib.pyplot as plt
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

PENALTIES = [
    ("Ll", "L l1"),
    ("Lc", "L corr."),
    ("Fc", "F corr."),
    ("Fo", "F orth."),
    ("Fa", "Final absorbing"),
    ("V", "Varimax"),
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot marginal runtime by penalty setting from GNU time files."
    )
    parser.add_argument("output_tsv")
    parser.add_argument("output_pdf")
    parser.add_argument("time_files", nargs="+")
    return parser.parse_args()


def parse_elapsed_to_seconds(text):
    match = re.search(r"(\S+)elapsed", text)
    if not match:
        raise ValueError("missing elapsed field")

    elapsed = match.group(1)
    days = 0
    if "-" in elapsed:
        day_text, elapsed = elapsed.split("-", 1)
        days = int(day_text)

    parts = elapsed.split(":")
    if len(parts) == 3:
        hours, minutes, seconds = parts
    elif len(parts) == 2:
        hours = 0
        minutes, seconds = parts
    elif len(parts) == 1:
        hours = 0
        minutes = 0
        seconds = parts[0]
    else:
        raise ValueError(f"unrecognized elapsed time: {elapsed}")

    return (
        days * 86400
        + int(hours) * 3600
        + int(minutes) * 60
        + float(seconds)
    )


def parse_condition(condition):
    prefixes = ["Fl2", "Ll", "Lc", "Fc", "Fo", "Fa", "V"]
    parts = {}
    for token in condition.split("_"):
        for prefix in prefixes:
            if token.startswith(prefix):
                parts[prefix] = token[len(prefix):]
                break
    return parts


def setting_label(key, value):
    return "Off" if str(value).lower() in {"0", "0.0", "false", "no"} else "On"


def load_runtime_rows(paths):
    rows = []
    for path_text in paths:
        path = Path(path_text)
        try:
            sec = parse_elapsed_to_seconds(path.read_text())
        except OSError as exc:
            sys.stderr.write(f"Could not read {path}: {exc}\n")
            return None
        except ValueError as exc:
            sys.stderr.write(f"Could not parse {path}: {exc}\n")
            return None

        prefix = path.with_suffix("")
        rows.append(
            {
                "clone": prefix.parent.name,
                "condition": prefix.name,
                "sec": sec,
            }
        )
    return rows


def build_long_table(df):
    long_rows = []
    for _, row in df.iterrows():
        parts = parse_condition(row["condition"])
        for key, label in PENALTIES:
            if key not in parts:
                continue
            long_rows.append(
                {
                    "clone": row["clone"],
                    "condition": row["condition"],
                    "sec": row["sec"],
                    "penalty": label,
                    "setting": setting_label(key, parts[key]),
                }
            )
    return pd.DataFrame(long_rows)


def main():
    args = parse_args()
    rows = load_runtime_rows(args.time_files)
    if rows is None:
        return 1
    if not rows:
        sys.stderr.write("No runtime files provided.\n")
        return 1

    df = pd.DataFrame(rows)
    df.to_csv(args.output_tsv, sep="\t", index=False)

    long_df = build_long_table(df)
    if long_df.empty:
        sys.stderr.write("No penalty settings found in runtime file names.\n")
        return 1

    penalty_order = [label for _, label in PENALTIES]
    setting_order = [value for value in ["Off", "On"] if value in set(long_df["setting"])]
    palette = {
        "Off": "#E69F00",
        "On": "#56B4E9",
    }

    fig, ax = plt.subplots(1, 1, figsize=(7.0, 3.5))
    sns.barplot(
        data=long_df,
        x="penalty",
        y="sec",
        hue="setting",
        order=penalty_order,
        hue_order=setting_order,
        palette={key: palette.get(key, "#009E73") for key in setting_order},
        errorbar="sd",
        capsize=0.08,
        err_kws={"linewidth": 1},
        ax=ax,
    )

    ax.set_xlabel("")
    ax.set_ylabel("Runtime (s)")
    ax.set_title("")
    ax.tick_params(axis="x", rotation=25)
    for label in ax.get_xticklabels():
        label.set_horizontalalignment("right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linewidth=0.5, alpha=0.3)
    ax.set_axisbelow(True)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend_.remove()
    fig.legend(
        handles,
        labels,
        title="Setting",
        loc="center left",
        frameon=False,
        bbox_to_anchor=(0.92, 0.5),
    )

    fig.savefig(args.output_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved runtimes to: {args.output_tsv}")
    print(f"Saved figure to: {args.output_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
