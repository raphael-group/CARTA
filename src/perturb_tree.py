from argparse import ArgumentParser

import cassiopeia as cas
import ete3
import networkx as nx
import numpy as np

parser = ArgumentParser()

parser.add_argument("--prefix", type=str, help="")
parser.add_argument("--tree_file_prefix", type=str, help="")
parser.add_argument("--location_file_prefix", type=str, help="")
parser.add_argument("--scaling_factor", type=float, help="")
parser.add_argument("--seed", type=int, help="")

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

def main():
    newick = f"{args.tree_file_prefix}_tree.txt"
    tree = ete3.Tree(newick, 1)

    g = nx.DiGraph()
    internal_node_iter = 0
    for n in tree.traverse():
        if n.name == "":
            n.name = f"cassiopeia_internal_node{internal_node_iter}"
            internal_node_iter += 1

        if n.is_root():
            continue

        g.add_edge(n.up.name, n.name, length = n.dist)
    tree = cas.data.CassiopeiaTree(tree = g)

    np.random.seed(args.seed)

    branch_length_mapping = {}
    for u, v in tree.edges:
        l = tree.get_branch_length(u, v)
        l_ = np.random.normal(loc=0.0, scale=args.scaling_factor * l) + l
        if l_ < 0:
            l_ = 0.01 * l
        branch_length_mapping[(u, v)] = l_
    tree.set_branch_lengths(branch_length_mapping)

    longest_time = 0
    for l in tree.leaves:
        leaf_time = tree.get_time(l)
        if leaf_time > longest_time:
            longest_time = leaf_time
            
    times_dict = {}
    for n in tree.nodes:
        if tree.is_leaf(n):
            times_dict[n] = 1
        else:
            times_dict[n] = tree.get_time(n)/longest_time
    tree.set_times(times_dict)

    out_newick = to_newick(tree._CassiopeiaTree__network, record_branch_lengths = True)

    out_file_loc = f"{args.prefix}_tree.txt"
    with open(out_file_loc, 'w') as f:
        f.write(out_newick)
    
    with open(f"{args.location_file_prefix}_locations.txt", "w") as f:
        f.write(f"{args.prefix}_tree.txt\t{args.tree_file_prefix}_meta.txt\n")


if __name__ == "__main__":
    main()
