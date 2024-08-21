from argparse import ArgumentParser

import cassiopeia as cas
import ete3
import networkx as nx
import pandas as pd

import cell_fate_mapping_evo_coupling as evo_coupling

parser = ArgumentParser()

parser.add_argument("--prefix", type=str, help="filepath prefix for this run")
parser.add_argument("--file_locations", type=str, help="a file that on each line specifies the filepath of a newick file and a corresponding state metadata file in that order, seperated by a tab")
parser.add_argument("--states_file", type=str, help="a file containing the cell states present in this dataset")
parser.add_argument(
    "--use_branch_lengths",
    default=False, action='store_true',
    help="whether to use branch lengths in computing cell distances",
)

args = parser.parse_args()

def main():
    states = []
    with open(args.states_file, "r") as f:
        for line in f:
            states.append(line.rstrip())

    with open(args.file_locations) as file_locations:
        for line in file_locations:
            nw, metadata = line.rstrip().split("\t")
            df_cellTypes = pd.read_csv(metadata, sep = "\t")
            df_cellTypes.set_index("cellBC", inplace = True)

    if args.use_branch_lengths:
        tree = ete3.Tree(nw, 1)

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
        tree = cas.data.CassiopeiaTree(tree = nw)

    out_newick = evo_coupling.avg_evolutionary_coupling(tree, states, df_cellTypes, args.use_branch_lengths)
    
    with open(f"{args.prefix}_EvoCGraph.nwk", "w") as f:
        f.write(out_newick)

if __name__ == "__main__":
    main()
