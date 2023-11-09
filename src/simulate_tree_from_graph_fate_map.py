from argparse import ArgumentParser

import json
import pandas as pd
import networkx as nx

import graph_generation_utils as utils

parser = ArgumentParser()

parser.add_argument("--prefix", type=str, help="filepath prefix for this run")
parser.add_argument("--location_file_prefix", type=str, help="filepath prefix for this run")
parser.add_argument("--graph_file_location", type=str, help="")
parser.add_argument("--num_sampled_per_cell_type", type=int, help="")
parser.add_argument("--subsample_rate", type=float, help="")
parser.add_argument("--seed", type=int, help="")

args = parser.parse_args()

def to_newick(tree, record_branch_lengths: bool = False) -> str:
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

def main():
    graph_dicts = []
    with open(args.graph_file_location) as graph_file:
        for line in graph_file:
            graph_dicts.append(json.loads(line))

    edge_dict = {}
    for key, value in graph_dicts[0].items():
        if key != "root":
            edge_dict[int(key)] = value
        else:
            edge_dict[key] = value
            
    total_times = {}
    for key, value in graph_dicts[1].items():
        if key != "root":
            total_times[int(key)] = value
        else:
            total_times[key] = value
            
    id_to_progen = {}
    for key, value in graph_dicts[2].items():
        if key != "root":
            id_to_progen[int(key)] = value
        else:
            id_to_progen[key] = value

    states = list(max(id_to_progen.values(), key=len))
    num_extant = int(args.num_sampled_per_cell_type/args.subsample_rate * len(states))
    
    tree, used_seed1 = utils.build_tree(num_extant, args.seed)
    tree, used_seed2 = utils.overlay_fate_map_on_tree(edge_dict, tree, states, total_times, args.num_sampled_per_cell_type, args.seed)
    utils.subsample_tree(tree, states, args.num_sampled_per_cell_type, args.seed)
    u = utils.count_unrealizations(tree, id_to_progen)

    with open(f"{args.prefix}_tree.txt", "w") as f:
        f.write(to_newick(tree._CassiopeiaTree__network, record_branch_lengths = True))

    meta = pd.DataFrame({"cellBC": tree.leaves, 
             "cell_state": [tree.get_attribute(l, "state_labels")[0] for l in tree.leaves]})
    meta.to_csv(f"{args.prefix}_meta.txt", sep='\t', index = False)

    with open(f"{args.prefix}_unrealizations.txt", "w") as f:
        f.write(f"{u}\n")

    with open(f"{args.prefix}_usedseeds.txt", "w") as f:
        f.write(f"{used_seed1}\t{used_seed2}\n")

    with open(f"{args.location_file_prefix}_locations.txt", "w") as f:
        f.write(f"{args.prefix}_tree.txt\t{args.prefix}_meta.txt\n")

if __name__ == "__main__":
    main()
