#!/usr/bin/env python3

"""Plot reconstructed active-factor values on four circular clone trees."""

import argparse
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import pandas as pd
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clone", action="append", nargs=3, required=True,
                        metavar=("LABEL", "LATENT_FLOW", "SUMMARY"))
    parser.add_argument("--active-variance-threshold", type=float, default=1e-5)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def active_factor_columns(summary_path, threshold):
    summary = pd.read_table(summary_path)
    factors = []
    for row in summary.itertuples(index=False):
        match = re.fullmatch(r"sigma2_latent_LF([0-9]+)", str(row.parameter))
        if match and float(row.value) > threshold:
            number = int(match.group(1))
            factors.append((number, f"factor_{number}"))
    if not factors:
        raise ValueError(f"No active factors found in {summary_path}")
    return factors


def circular_layout(frame):
    node_ids = frame["node_id"].astype(int).tolist()
    node_set = set(node_ids)
    children = {node_id: [] for node_id in node_ids}
    parent = {}
    for row in frame[["node_id", "parent_id"]].itertuples(index=False):
        node_id, parent_id = int(row.node_id), int(row.parent_id)
        parent[node_id] = parent_id
        if parent_id >= 0:
            if parent_id not in node_set:
                raise ValueError(f"Parent {parent_id} is absent")
            children[parent_id].append(node_id)
    for child_list in children.values():
        child_list.sort()

    roots = [node_id for node_id in node_ids if parent[node_id] < 0]
    if len(roots) != 1:
        raise ValueError("Expected exactly one root")
    root_id = roots[0]

    ordered_tips, postorder = [], []
    stack = [(root_id, False)]
    while stack:
        node_id, visited = stack.pop()
        if visited:
            postorder.append(node_id)
            continue
        stack.append((node_id, True))
        for child_id in reversed(children[node_id]):
            stack.append((child_id, False))
        if not children[node_id]:
            ordered_tips.append(node_id)
    if not ordered_tips:
        raise ValueError("No terminal nodes found")

    tip_position = {tip: i for i, tip in enumerate(ordered_tips)}
    descendant_span = {tip: (tip_position[tip], tip_position[tip])
                       for tip in ordered_tips}
    for node_id in postorder:
        if node_id not in descendant_span:
            descendant_span[node_id] = (
                min(descendant_span[child][0] for child in children[node_id]),
                max(descendant_span[child][1] for child in children[node_id]),
            )
    theta = {
        node_id: 2 * np.pi * (0.5 * (span[0] + span[1]) + 0.5) / len(ordered_tips)
        for node_id, span in descendant_span.items()
    }

    depth = dict(zip(frame["node_id"].astype(int), frame["tree_depth"].astype(float)))
    max_depth = max(depth.values())
    if max_depth <= 0:
        raise ValueError("Tree depth must be positive")
    radius = {node_id: depth[node_id] / max_depth for node_id in node_ids}

    segments, segment_children = [], []
    for child_id in node_ids:
        parent_id = parent[child_id]
        if parent_id < 0:
            continue
        # Draw the orthogonal circular-tree edge as an arc at the parent's
        # radius followed by a radial daughter branch at the child's angle.
        angle_delta = theta[child_id] - theta[parent_id]
        n_arc_segments = max(1, min(16, int(np.ceil(abs(angle_delta) / (np.pi / 64)))))
        arc_angles = np.linspace(theta[parent_id], theta[child_id], n_arc_segments + 1)
        arc_points = np.column_stack([
            radius[parent_id] * np.sin(arc_angles),
            radius[parent_id] * np.cos(arc_angles),
        ])
        child_arc = arc_points[-1]
        child_tip = np.array([radius[child_id] * np.sin(theta[child_id]),
                              radius[child_id] * np.cos(theta[child_id])])
        if radius[parent_id] > 0 and abs(angle_delta) > 1e-12:
            for index in range(n_arc_segments):
                segments.append([arc_points[index], arc_points[index + 1]])
                segment_children.append(child_id)
        segments.append([child_arc, child_tip])
        segment_children.append(child_id)
    return root_id, parent, np.asarray(segments), segment_children


def read_clone(label, flow_path, summary_path, threshold):
    frame = pd.read_table(flow_path)
    factors = active_factor_columns(summary_path, threshold)
    missing = [column for _number, column in factors if column not in frame.columns]
    if missing:
        raise ValueError(f"Active factors absent for {label}: {', '.join(missing)}")
    root_id, parent, segments, segment_children = circular_layout(frame)
    values = frame.copy()
    values.index = frame["node_id"].astype(int)
    factor_matrix = values[[column for _number, column in factors]].to_numpy(float)
    root_position = values.index.get_loc(root_id)
    values["total_distance"] = np.linalg.norm(
        factor_matrix - factor_matrix[root_position], axis=1
    )
    return dict(label=label, factors=factors, values=values, parent=parent,
                segments=segments, segment_children=segment_children)


def edge_values(clone, column):
    values = clone["values"][column].to_dict()
    return np.asarray([
        0.5 * (float(values[clone["parent"][child]]) + float(values[child]))
        for child in clone["segment_children"]
    ])


def draw_tree(axis, clone, column, title, sequential=False, sequential_vmax=None):
    values = clone["values"][column].to_numpy(float)
    if sequential:
        vmax = (
            float(sequential_vmax)
            if sequential_vmax is not None
            else float(np.nanmax(values))
        )
        norm = mpl.colors.Normalize(0, max(vmax, 1e-12))
        cmap = mpl.colormaps["YlOrBr"]
    else:
        bound = max(float(np.nanmax(np.abs(values))), 1e-12)
        norm = mpl.colors.TwoSlopeNorm(vmin=-bound, vcenter=0, vmax=bound)
        cmap = mpl.colormaps["RdBu_r"]
    collection = LineCollection(clone["segments"], cmap=cmap, norm=norm,
                                linewidths=0.45, alpha=0.95, rasterized=True)
    collection.set_array(edge_values(clone, column))
    axis.add_collection(collection)
    axis.set(xlim=(-1.04, 1.04), ylim=(-1.04, 1.04), aspect="equal")
    axis.axis("off")
    axis.set_title(title, fontsize=13, pad=5)
    color_axis = axis.inset_axes([0.20, -0.035, 0.60, 0.025])
    colorbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                            cax=color_axis, orientation="horizontal")
    colorbar.ax.tick_params(labelsize=9, length=3, pad=2)
    colorbar.outline.set_linewidth(0.4)


def main():
    args = parse_args()
    clones = [read_clone(label, flow, summary, args.active_variance_threshold)
              for label, flow, summary in args.clone]
    max_factor = max(number for clone in clones for number, _ in clone["factors"])
    distance_vmax = max(
        float(clone["values"]["total_distance"].max()) for clone in clones
    )
    sns.set_theme(style="white", context="paper")
    mpl.rcParams.update({"figure.dpi": 140, "savefig.dpi": 300})
    figure, axes = plt.subplots(len(clones), max_factor + 1,
                                figsize=(2.55 * (max_factor + 1), 2.75 * len(clones)),
                                squeeze=False)
    for row, clone in enumerate(clones):
        draw_tree(
            axes[row, 0], clone, "total_distance", "Euclidean distance",
            sequential=True, sequential_vmax=distance_vmax,
        )
        factor_by_number = dict(clone["factors"])
        for number in range(1, max_factor + 1):
            if number in factor_by_number:
                draw_tree(axes[row, number], clone, factor_by_number[number],
                          f"Factor {number}")
            else:
                axes[row, number].axis("off")
        axes[row, 0].text(-0.20, 0.50, clone["label"],
                          transform=axes[row, 0].transAxes, ha="right", va="center",
                          fontsize=15, fontweight="normal")
    figure.subplots_adjust(left=0.075, right=0.995, bottom=0.045, top=0.975,
                           wspace=0.18, hspace=0.30)
    figure.savefig(args.output, bbox_inches="tight", facecolor="white")
    plt.close(figure)


if __name__ == "__main__":
    main()
