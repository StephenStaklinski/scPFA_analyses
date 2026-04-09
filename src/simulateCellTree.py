
import argparse
import cassiopeia as cas
import numpy as np


def positive_normal(scale):
    """Sample from a normal distribution and return the absolute value."""
    return abs(np.random.normal(loc=scale, scale=scale / 5))


parser = argparse.ArgumentParser(description="Simulate CRISPR barcodes and tree.")
parser.add_argument("--out_tree", type=str, required=True, help="Output path for ground truth tree")
parser.add_argument("--birth_rate", type=float, default=0.4, help="Birth rate for tree process")
parser.add_argument("--death_rate", type=float, default=0.05, help="Death rate for tree process")
parser.add_argument("--num_tips", type=int, default=50, help="Number of tips in the tree")
parser.add_argument("--desired_time", type=int, default=50, help="Desired time for the tree height. Note: this will lead to uniform scaling that alters the scale of birth and death rates.")
args = parser.parse_args()

out_tree = args.out_tree
birth_rate = args.birth_rate
death_rate = args.death_rate
num_tips = args.num_tips
desired_time = args.desired_time

# Commenting this out in favor of using a fixed input tree instead of simulating the tree here
# Simulate tree using birth-death process
simulator = cas.sim.BirthDeathFitnessSimulator(
    birth_waiting_distribution=lambda scale: positive_normal(scale),
    initial_birth_scale=1 / birth_rate,
    death_waiting_distribution=lambda: positive_normal(1 / death_rate),
    num_extant=num_tips
)
ground_truth_tree = simulator.simulate_tree()

# Uniformly scale the tree to the desired height
current_heights = ground_truth_tree.get_times().values()
scaling_factor = desired_time / max(current_heights)
new_times = {node: ground_truth_tree.get_time(node) * scaling_factor for node in ground_truth_tree.nodes}
ground_truth_tree.set_times(new_times)

# Write ground truth tree to file
nwk = ground_truth_tree.get_newick(record_branch_lengths=True, record_node_names=False)
with open(out_tree, "w") as f:
    f.write(nwk)
