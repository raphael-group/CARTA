from argparse import ArgumentParser

import cassiopeia as cas
import networkx as nx
import pandas as pd
import pickle as pic

import utils
import ilp

parser = ArgumentParser()

parser.add_argument("--prefix", type=str, help="filepath prefix for this run")
parser.add_argument("--file_locations", type=str, help="a file that on each line specifies the filepath of a newick file and a corresponding state metadata file in that order, seperated by a tab")
parser.add_argument("--states_file", type=str, help="a file containing the cell states present in this dataset")
parser.add_argument("--node_labels_file", type=str, help="a file where each line provides a tree index, a node name, and the states that node is labeled with")
parser.add_argument(
    "--edge_weight_scheme",
    default="raw_frequency",
    choices=["raw_frequency", "weighted_by_subtree_size"],
    help="the scheme by which to weight frequencies/edge weights in the output graph",
)

args = parser.parse_args()

def main():
    states = []
    with open(args.states_file, "r") as f:
        for line in f:
            states.append(line.rstrip())

    node_labels_df = pd.read_table(args.node_labels_file)

    line_ind = 0
    labeled_trees = []
    with open(args.file_locations) as file_locations:
        for line in file_locations:
            nw, metadata = line.rstrip().split("\t")
            tree = cas.data.CassiopeiaTree(tree = nw)
            utils.label_tree_with_leaf_states(tree, metadata)
            utils.prune_unwanted_states(tree, states)
            utils.impute_states_from_children(tree)
            
            for index, row in node_labels_df[node_labels_df["tree_index"] == line_ind].iterrows():
                states_at_node = row["label"][1:-1].replace(" ", "").replace("'", "").split(",")
                tree.set_attribute(row["node_name"], "state_labels", states_at_node)
            
            labeled_trees.append(tree)
            line_ind += 1
            
    if args.edge_weight_scheme == "raw_frequency":
        transition_counts = utils.get_transition_frequencies_over_trees(labeled_trees, states)
    elif args.edge_weight_scheme == "weighted_by_subtree_size":
        transition_counts = utils.get_transition_frequencies_weighted_by_subtree_size(labeled_trees, states)

    G = utils.generate_hasse_diagram(states)
    G = utils.annotate_hasse_with_transitions(G, states, transition_counts, remove_zero_weight_edges=True)

    pic.dump(G, open(f"{args.prefix}_cell_fate_map.pkl", "wb"))

    with open(f"{args.prefix}_cell_fate_map.cyjs", "w") as f:
        out = str(nx.cytoscape_data(G)).replace('\'', '\"').replace("True", "true").replace("False", "false")
        f.write(out)

if __name__ == "__main__":
    main()
