from argparse import ArgumentParser

import cassiopeia as cas
from cassiopeia.tl.branch_length_estimator import IIDExponentialMLE
import ete3
import networkx as nx
import numpy as np

parser = ArgumentParser()

parser.add_argument("--input_prefix", type=str, help="filepath prefix for the input")
parser.add_argument("--output_prefix", type=str, help="filepath prefix for the output")

args = parser.parse_args()

def to_newick(tree: nx.DiGraph, record_branch_lengths: bool = False) -> str:
    """Converts a networkx graph to a newick string.

    Args:
        tree: A networkx tree
        record_branch_lengths: Whether to record branch lengths on the tree in
            the newick string

    Returns:
        A newick string representing the topology of the tree
    """

    def _to_newick_str(g, node):
        is_leaf = g.out_degree(node) == 0
        weight_string = ""

        if record_branch_lengths and g.in_degree(node) > 0:
            parent = list(g.predecessors(node))[0]
            weight_string = ":" + str(g[parent][node]["length"])

        _name = str(node)
        return (
            "%s" % (_name,) + weight_string
            if is_leaf
            else (
                "("
                + ",".join(
                    _to_newick_str(g, child) for child in g.successors(node)
                )
                + ")"
                + "%s" % (_name,)
                + weight_string
            )
        )

    root = [node for node in tree if tree.in_degree(node) == 0][0]
    return _to_newick_str(tree, root) + ";"

def generate_cm_and_time_estimated_newick(newick):
   
    # round 1
    #total_time = 11
    #num_sites = 50
    #num_states = 1000
    #percent_mutated = 0.5
    #percent_dropped_h = 0.1
    #percent_dropped_s = 0.11111

    # round 2
    #total_time = 11
    #num_sites = 30
    #num_states = 32
    #percent_mutated = 0.1
    #percent_dropped_h = 0.2
    #percent_dropped_s = 0.15
    
    # round 2
    total_time = 11
    num_sites = 30
    num_states = 32
    percent_mutated = 0.1
    percent_dropped_h = 0.25
    percent_dropped_s = 0.2

    tree = ete3.Tree(newick, 1)
    g = nx.DiGraph()
    internal_node_iter = 0
    for n in tree.traverse():
        if n.name == "":
            n.name = f"cassiopeia_internal_node{internal_node_iter}"
            internal_node_iter += 1

        if n.is_root():
            continue

        g.add_edge(n.up.name, n.name, length = n.dist/total_time)
    tree = cas.data.CassiopeiaTree(tree = g)
    
    sim = cas.simulator.Cas9LineageTracingDataSimulator(
        number_of_cassettes = num_sites,
        size_of_cassette = 1,
        mutation_rate = -np.log(1 - percent_mutated),
        number_of_states = num_states,
        heritable_silencing_rate = -np.log(1 - percent_dropped_heritable),
        stochastic_silencing_rate = (total_percent_dropped - percent_dropped_heritable)/(1 - percent_dropped_heritable)
    )
    sim.overlay_data(tree)

    ble = IIDExponentialMLE(minimum_branch_length = 0.01)
    tree.reconstruct_ancestral_characters()
    ble.estimate_branch_lengths(tree)

    branch_length_dict = {}
    for e1, e2 in tree.edges:
        branch_length_dict[(e1, e2)] = tree.get_branch_length(e1, e2) * total_time
    tree.set_branch_lengths(branch_length_dict)

    return to_newick(tree._CassiopeiaTree__network, record_branch_lengths = True), tree.character_matrix


def main():
    file_loc = f"{args.input_prefix}_tree.nwk"
    with open(file_loc, 'r') as f:
        for line in f:
            newick = line
            
    out_newick, cm = generate_cm_and_time_estimated_newick(newick)
    out_file_loc = f"{args.output_prefix}_tree.nwk"
    with open(out_file_loc, 'w') as f:
        f.write(out_newick)
    
    cm.to_csv(f"{args.output_prefix}_cm.txt", sep='\t', header = True)

if __name__ == "__main__":
    main()
