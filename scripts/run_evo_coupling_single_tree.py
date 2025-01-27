from argparse import ArgumentParser

import cassiopeia as cas
import ete3
import networkx as nx
import pandas as pd

import sys
import os
sys.path.insert(0, os.path.abspath('/n/fs/ragr-data/users/palash/carta-rebuttal/src'))
import evo_coupling

parser = ArgumentParser()

parser.add_argument("--prefix", type=str, help="filepath prefix for this run")
parser.add_argument("--tree_file", type=str, help="")
parser.add_argument("--label_file", type=str, help="")
parser.add_argument(
    "--use_branch_lengths",
    default=False, action='store_true',
    help="whether to use branch lengths in computing cell distances",
)

args = parser.parse_args()

def main():

    df_cellTypes = pd.read_csv(args.label_file, sep = "\t")
    df_cellTypes.set_index("cellBC", inplace = True)

    states = list(set(df_cellTypes["cell_state"]))

    if args.use_branch_lengths:
        tree = ete3.Tree(args.tree_file, 1)

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
    else:
        tree = cas.data.CassiopeiaTree(tree = args.tree_file)

    out_newick = evo_coupling.avg_evolutionary_coupling(tree, states, df_cellTypes, args.use_branch_lengths)
    
    with open(f"{args.prefix}_EvoCGraph.nwk", "w") as f:
        f.write(out_newick)

if __name__ == "__main__":
    main()
