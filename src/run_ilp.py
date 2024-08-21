from argparse import ArgumentParser

import cassiopeia as cas
import cell_fate_mapping_utils as utils
import cell_fate_mapping_ilp as ilp

parser = ArgumentParser()

parser.add_argument("--prefix", type=str, help="filepath prefix for this run")
parser.add_argument("-k", type=int, help="number of a prior progenitors")
parser.add_argument("--file_locations", type=str, help="a file that on each line specifies the filepath of a newick file and a corresponding state metadata file in that order, seperated by a tab")
parser.add_argument("--states_file", type=str, help="a file containing the cell states present in this dataset")
parser.add_argument(
    "--normalize_method",
    default="no_normalization",
    choices=["no_normalization", "cell_proportion_before_pruning", "cell_proportion_after_pruning"],
    help="the scheme by which to normalize the objective function",
)
parser.add_argument(
    "--time_limit_sec",
    default=2880000,
    help="time limit on ilp runs",
)
parser.add_argument("--enforce_tree", default=False, action='store_true', help="whether or not to enforce that the progenitors can create a tree"
)

args = parser.parse_args()

def main():
    states = []
    with open(args.states_file, "r") as f:
        for line in f:
            states.append(line.rstrip())

    labeled_trees = []
    with open(args.file_locations) as file_locations:
        for line in file_locations:
            nw, metadata = line.rstrip().split("\t")
            print(nw, metadata)
            tree = cas.data.CassiopeiaTree(tree = nw)
            utils.label_tree_with_leaf_states(tree, metadata)
            utils.prune_unwanted_states(tree, states)
            # utils.impute_states_from_children(tree)
            labeled_trees.append(tree)

    if args.normalize_method == "no_normalization":
        weights = None
    elif args.normalize_method == "cell_proportion_before_pruning":
        weights = ilp.get_inverse_state_proportions_from_trees(labeled_trees)
    elif args.normalize_method == "cell_proportion_after_pruning":
        weights = ilp.get_inverse_state_proportions_from_trees_post_pruning(labeled_trees)

    if args.enforce_tree:
        # solved_model, observed_potencies = ilp.solve_large_k_problem_tree(
        #     labeled_trees, states, args.k, weights, 
        # )
        # out = ilp.post_process_solution_tree(solved_model, observed_potencies, states, labeled_trees)
        solved_model = ilp.solve_large_k_problem_tree(
            labeled_trees, states, args.k, weights, args.time_limit_sec
        )
        out = ilp.post_process_solution_tree(solved_model, states)

    else:
        solved_model = ilp.solve_large_k_problem(
            labeled_trees, states, args.k, weights, args.time_limit_sec
        )
        out = ilp.post_process_solution_from_node_state_labels(solved_model)
    
    with open(f"{args.prefix}_results.txt", "w") as f:
        f.write("objective_score\tmodel_runtime\n")
        f.write(f"{out[0]}\t{solved_model.runtime}\n")

    with open(f"{args.prefix}_progenitors.txt", "w") as f:
        for progen in out[1]:
            f.write(f"{progen}\n")

    # with open(f"{args.prefix}_nodeLabels.txt", "w") as f:
    #     f.write("tree_index\tnode_name\tlabel\n")
    #     for node, label in out[2].items():
    #         tree_index = node.split("-")[-1]
    #         node_name = "-".join(node.split("-")[:-1])
    #         f.write(f"{tree_index}\t{node_name}\t{label}\n")


if __name__ == "__main__":
    main()
