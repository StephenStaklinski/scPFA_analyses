#!/usr/bin/env python3

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D


def branch_times(branches):
    """Recover root-to-node time from the parent/child edge table."""
    parent = dict(zip(branches.child_id, branches.parent_id))
    length = dict(zip(branches.child_id, branches.branch_length))
    children = set(parent)
    roots = set(branches.parent_id) - children
    if len(roots) != 1:
        raise ValueError("Expected exactly one root in each module's branch table")
    root = roots.pop()
    cache = {root: 0.0}

    def depth(node):
        if node not in cache:
            cache[node] = depth(parent[node]) + length[node]
        return cache[node]

    parent_depth = np.array([depth(node) for node in branches.parent_id])
    child_depth = np.array([depth(node) for node in branches.child_id])
    scale = max(child_depth.max(), 1e-12)
    return (parent_depth + child_depth) / (2.0 * scale)


def tree_layout(branches):
    children = {}
    length = {}
    for row in branches.itertuples():
        children.setdefault(row.parent_id, []).append(row.child_id)
        length[row.child_id] = row.branch_length
    child_nodes = set(length)
    roots = set(children) - child_nodes
    if len(roots) != 1:
        raise ValueError("Expected exactly one root in each module's branch table")
    root = roots.pop()
    depth = {root: 0.0}

    def set_depths(node):
        for child in children.get(node, []):
            depth[child] = depth[node] + length[child]
            set_depths(child)

    set_depths(root)
    scale = max(depth.values()) or 1.0
    depth = {node: value / scale for node, value in depth.items()}

    y = {}
    next_tip = [0]

    def set_y(node):
        descendants = sorted(children.get(node, []))
        if not descendants:
            y[node] = next_tip[0]
            next_tip[0] += 1
        else:
            for child in descendants:
                set_y(child)
            y[node] = np.mean([y[child] for child in descendants])

    set_y(root)
    return children, depth, y, next_tip[0]


def summarize_timing(branches, bins):
    rows = []
    profiles = {}
    edges = np.linspace(0.0, 1.0, bins + 1)
    for module, frame in branches.groupby("module", sort=True):
        frame = frame.loc[frame.shift_eligible.astype(bool)].copy()
        frame["time"] = branch_times(branches.loc[branches.module == module])[
            branches.loc[branches.module == module, "shift_eligible"].astype(bool).to_numpy()
        ]
        frame["signed_probability"] = np.where(
            frame.reconstructed_delta >= 0,
            frame.shift_probability,
            -frame.shift_probability,
        )
        positive, _ = np.histogram(
            frame.time, edges,
            weights=np.where(frame.signed_probability > 0, frame.shift_probability, 0.0),
        )
        negative, _ = np.histogram(
            frame.time, edges,
            weights=np.where(frame.signed_probability < 0, frame.shift_probability, 0.0),
        )
        total = positive + negative
        opportunities, _ = np.histogram(frame.time, edges)
        mass = total.sum()
        if mass > 0:
            posterior_timing = total / mass
            diffuse_timing = opportunities / opportunities.sum()
            burstiness = 0.5 * np.abs(posterior_timing - diffuse_timing).sum()
        else:
            burstiness = np.nan
        rows.append({
            "module": module,
            "posterior_expected_shifts": mass,
            "positive_shift_mass": positive.sum(),
            "negative_shift_mass": negative.sum(),
            "temporal_burstiness": burstiness,
            "branches_probability_ge_0.5": int((frame.shift_probability >= 0.5).sum()),
        })
        profiles[module] = (frame, edges, positive, negative, opportunities)
    return pd.DataFrame(rows), profiles


def make_plots(branches, summary, timing, profiles, output):
    plt.rcParams.update({"font.size": 11, "axes.titlesize": 13, "axes.labelsize": 11})
    merged = timing.merge(
        summary[["module", "posterior_expected_shifts"]],
        on="module", how="left", suffixes=("_binned", "_model"),
    )
    with PdfPages(output) as pdf:
        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        ax.scatter(
            merged.posterior_expected_shifts_binned,
            merged.temporal_burstiness,
            s=65, color="#5B4B8A", alpha=0.85,
        )
        for row in merged.itertuples():
            ax.annotate(f"Module {row.module}",
                        (row.posterior_expected_shifts_binned, row.temporal_burstiness),
                        xytext=(5, 4), textcoords="offset points")
        ax.set_xlabel("Posterior expected number of shifts")
        ax.set_ylabel("Timing departure from diffuse (0–1)")
        ax.set_title("Module shift burden and timing")
        ax.set_ylim(bottom=0)
        ax.grid(alpha=0.2)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        for module in sorted(profiles):
            frame, edges, positive, negative, opportunities = profiles[module]
            centers = (edges[:-1] + edges[1:]) / 2
            width = np.diff(edges) * 0.88
            metric = timing.loc[timing.module == module].iloc[0]
            diffuse = opportunities / opportunities.sum() * metric.posterior_expected_shifts
            fig, (ax_time, ax_event) = plt.subplots(
                2, 1, figsize=(8.2, 7.0), sharex=True,
                gridspec_kw={"height_ratios": [1.05, 1.0]},
            )
            ax_time.bar(centers, positive, width=width, color="#D55E00",
                        label="Positive change")
            ax_time.bar(centers, -negative, width=width, color="#0072B2",
                        label="Negative change")
            ax_time.plot(centers, diffuse, color="0.3", linestyle="--", marker="o",
                         markersize=3, linewidth=1.2, label="Diffuse baseline")
            ax_time.axhline(0, color="0.25", linewidth=0.8)
            ax_time.set_ylabel("Posterior expected shifts\nper time bin")
            ax_time.set_title(
                f"Module {module}: {metric.posterior_expected_shifts:.2f} expected shifts; "
                f"timing departure = {metric.temporal_burstiness:.2f}"
            )
            ax_time.legend(frameon=False, ncol=2)

            sizes = 18 + 150 * frame.shift_probability
            colors = np.where(frame.reconstructed_delta >= 0, "#D55E00", "#0072B2")
            ax_event.scatter(frame.time, frame.reconstructed_delta, s=sizes, c=colors,
                             alpha=np.clip(0.2 + 0.8 * frame.shift_probability, 0, 1),
                             edgecolors="white", linewidths=0.4)
            ax_event.axhline(0, color="0.25", linewidth=0.8)
            ax_event.set_xlabel("Relative time from root (root = 0, tips = 1)")
            ax_event.set_ylabel("Reconstructed module change")
            ax_event.set_xlim(0, 1)
            ax_event.grid(axis="x", alpha=0.2)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def make_tree_plots(branches, output, minimum_probability):
    with PdfPages(output) as pdf:
        for module, frame in branches.groupby("module", sort=True):
            frame = frame.copy()
            children, depth, y, tip_count = tree_layout(frame)
            horizontal = []
            vertical = []
            for parent, child_list in children.items():
                for child in child_list:
                    horizontal.append([(depth[parent], y[child]), (depth[child], y[child])])
                if len(child_list) > 1:
                    child_y = [y[child] for child in child_list]
                    vertical.append([(depth[parent], min(child_y)),
                                     (depth[parent], max(child_y))])

            fig, ax = plt.subplots(figsize=(10.0, 7.0))
            ax.add_collection(LineCollection(horizontal, colors="0.72", linewidths=0.35,
                                             rasterized=True, zorder=1))
            ax.add_collection(LineCollection(vertical, colors="0.72", linewidths=0.35,
                                             rasterized=True, zorder=1))

            shown = frame.loc[
                frame.shift_eligible.astype(bool)
                & (frame.shift_probability >= minimum_probability)
            ].copy()
            shown["x"] = [
                (depth[parent] + depth[child]) / 2
                for parent, child in zip(shown.parent_id, shown.child_id)
            ]
            shown["y"] = [y[child] for child in shown.child_id]
            colors = np.where(shown.reconstructed_delta >= 0, "#D55E00", "#0072B2")
            ax.scatter(
                shown.x, shown.y,
                s=12 + 130 * shown.shift_probability,
                c=colors,
                alpha=0.85,
                edgecolors="white",
                linewidths=0.35,
                rasterized=True,
                zorder=3,
            )

            direction_legend = [
                Line2D([], [], marker="o", linestyle="", color="#D55E00",
                       markersize=7, label="Positive change"),
                Line2D([], [], marker="o", linestyle="", color="#0072B2",
                       markersize=7, label="Negative change"),
            ]
            size_legend = [
                ax.scatter([], [], s=12 + 130 * probability, facecolor="0.45",
                           edgecolor="white", label=f"P(shift) = {probability:g}")
                for probability in (0.25, 0.5, 0.9)
                if probability >= minimum_probability
            ]
            first = ax.legend(handles=direction_legend, loc="upper left",
                              frameon=False, title="Direction")
            ax.add_artist(first)
            ax.legend(handles=size_legend, loc="lower left",
                      frameon=False, title="Posterior probability")
            ax.set_xlim(-0.01, 1.01)
            ax.set_ylim(-1, tip_count)
            ax.set_yticks([])
            ax.set_xlabel("Relative time from root (root = 0, tips = 1)")
            ax.set_ylabel(f"{tip_count:,} terminal lineages")
            ax.set_title(
                f"Module {module}: posterior shift locations "
                f"(P ≥ {minimum_probability:g})"
            )
            ax.spines[["left", "right", "top"]].set_visible(False)
            fig.tight_layout()
            pdf.savefig(fig, dpi=250)
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--branches", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--tree-output", required=True)
    parser.add_argument("--timing-output", required=True)
    parser.add_argument("--bins", type=int, default=12)
    parser.add_argument("--tree-min-probability", type=float, default=0.1)
    args = parser.parse_args()
    branches = pd.read_csv(args.branches, sep="\t")
    summary = pd.read_csv(args.summary, sep="\t")
    timing, profiles = summarize_timing(branches, args.bins)
    timing.to_csv(args.timing_output, sep="\t", index=False, float_format="%.8g")
    make_plots(branches, summary, timing, profiles, args.output)
    make_tree_plots(branches, args.tree_output, args.tree_min_probability)


if __name__ == "__main__":
    main()
