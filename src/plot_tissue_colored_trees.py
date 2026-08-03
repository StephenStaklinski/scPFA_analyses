#!/usr/bin/env python3

import argparse
from collections import defaultdict

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from Bio import Phylo


TISSUE_COLORS = {
    "M1": "#4c78a8",
    "M2": "#f58518",
    "RW": "#54a24b",
    "RE": "#e45756",
    "LL": "#b279a2",
    "Liv": "#ffbf79",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot phylogenies with terminal branches colored by tissue."
    )
    parser.add_argument(
        "--tree",
        nargs=2,
        action="append",
        required=True,
        metavar=("LABEL", "NEXUS"),
        help="Panel label and NEXUS tree path; repeat once per panel.",
    )
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def tissue_from_tip(name):
    return str(name).split(".", 1)[0]


def tree_segments(tree):
    depths = tree.depths()
    if not depths or max(depths.values()) == 0:
        depths = tree.depths(unit_branch_lengths=True)

    terminals = tree.get_terminals()
    theta = {
        tip: 2.0 * np.pi * (index + 0.5) / len(terminals)
        for index, tip in enumerate(terminals)
    }

    def assign_internal_angle(clade):
        if clade in theta:
            return theta[clade]
        child_angles = [assign_internal_angle(child) for child in clade.clades]
        theta[clade] = float(np.mean(child_angles))
        return theta[clade]

    assign_internal_angle(tree.root)
    max_depth = max(depths.values())
    radius = {clade: depths[clade] / max_depth for clade in depths}

    internal_segments = []
    terminal_segments = defaultdict(list)
    tissues = set()

    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            segments = []
            if abs(theta[child] - theta[parent]) > 1e-12 and radius[parent] > 0:
                arc_theta = np.linspace(theta[parent], theta[child], 12)
                segments.append(
                    np.column_stack(
                        [arc_theta, np.full_like(arc_theta, radius[parent])]
                    )
                )
            segments.append(
                np.array([
                    [theta[child], radius[parent]],
                    [theta[child], radius[child]],
                ])
            )
            if child.is_terminal():
                tissue = tissue_from_tip(child.name)
                terminal_segments[tissue].extend(segments)
                tissues.add(tissue)
            else:
                internal_segments.extend(segments)

    return internal_segments, terminal_segments, tissues


def main():
    args = parse_args()
    mpl.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
    })

    figure, axes = plt.subplots(
        1, len(args.tree),
        figsize=(4.0 * len(args.tree), 4.5),
        squeeze=False,
        subplot_kw={"projection": "polar"},
    )
    all_tissues = set()

    for axis, (label, tree_path) in zip(axes[0], args.tree):
        tree = Phylo.read(tree_path, "nexus")
        internal, terminal, tissues = tree_segments(tree)
        all_tissues.update(tissues)

        axis.add_collection(
            LineCollection(
                internal, colors="0.72", linewidths=0.22, alpha=0.55,
                rasterized=True, transform=axis.transData,
            )
        )
        for tissue, segments in terminal.items():
            axis.add_collection(
                LineCollection(
                    segments,
                    colors=TISSUE_COLORS.get(tissue, "0.25"),
                    linewidths=0.48,
                    alpha=0.9,
                    rasterized=True,
                    transform=axis.transData,
                )
            )

        axis.set_theta_direction(-1)
        axis.set_theta_offset(np.pi / 2.0)
        axis.set_ylim(0, 1.03)
        axis.set_xticks([])
        axis.set_yticks([])
        axis.grid(False)
        axis.spines["polar"].set_visible(False)
        axis.set_title(label, pad=8, fontweight="normal")

    legend_order = [
        tissue for tissue in TISSUE_COLORS if tissue in all_tissues
    ] + sorted(all_tissues - set(TISSUE_COLORS))
    handles = [
        Line2D(
            [0], [0], color=TISSUE_COLORS.get(tissue, "0.25"),
            linewidth=2.2, label=tissue,
        )
        for tissue in legend_order
    ]
    figure.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.03),
        ncol=len(handles),
        frameon=False,
        title="Tissue",
    )
    figure.tight_layout(rect=(0, 0.11, 1, 1))
    figure.savefig(args.output, bbox_inches="tight")
    plt.close(figure)


if __name__ == "__main__":
    main()
