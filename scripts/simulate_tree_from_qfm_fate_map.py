from argparse import ArgumentParser

import json
import pandas as pd
import networkx as nx
import numpy as np

import graph_generation_utils as utils
import MLE_time_estimate_trees as MLE

parser = ArgumentParser()

parser.add_argument("--prefix", type=str, help="filepath prefix for this run")
parser.add_argument("--tree_ind", type=int, help="")
parser.add_argument("--graph_file_location", type=str, help="")
parser.add_argument("--num_sampled_per_cell_type", type=int, help="", default = 100)
parser.add_argument("--min_sampled_per_cell_type", type=int, help="", default = 100)
parser.add_argument("--subsample_rate", type=float, help="", default = 0.1)
parser.add_argument("--seed", type=int, help="", default = None)
parser.add_argument("--graph_format", choices=['qfm', 'carta'], help="", default = 'qfm')

args = parser.parse_args()

def annotate_state(node_id, transitions, potencies):
    if len(transitions[node_id]) == 0:
        potencies[node_id] = node_id
        return [node_id]
    else:
        all_states = []
        for i in transitions[node_id]:
            all_states.extend(annotate_state(i, transitions, potencies))
        potencies[node_id] = all_states
        return all_states

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

    if args.graph_format == 'qfm':
        # Open the JSON file
        with open(args.graph_file_location) as f:
            # Load the JSON data as a dictionary
            data = json.load(f)

        edge_dict = data["merge"]
        edge_dict['root'] = [data['root_id']]
        for i in data["tip_id"]:
            edge_dict[i] = []

        id_to_progen = {}
        annotate_state(data['root_id'], edge_dict, id_to_progen)
        tot_time = data["target_time"]

        total_times = dict([(u, t/data['target_time']) for u, t in data["diff_time"].items() if t != 'Inf'])
        for i in data["tip_id"]:
            total_times[i] = 1.0

        prob_dict = dict(zip(data['diff_mode_probs'].keys(),[i[:2] for i in data['diff_mode_probs'].values()]))
        
        root_state = data["root_id"]

        incidence_dict = dict(zip(data['tip_id'], [1] * len(data['tip_id'])))

    else:
        graph_dicts = []
        with open(args.graph_file_location) as graph_file:
            for line in graph_file:
                graph_dicts.append(json.loads(line))

        edge_dict = {}
        for key, value in graph_dicts[0].items():
            edge_dict[str(key)] = value
                
        total_times = {}
        for key, value in graph_dicts[1].items():
            total_times[str(key)] = value
                
        id_to_progen = {}
        for key, value in graph_dicts[2].items():
            id_to_progen[str(key)] = value
        
        root_state = str(graph_dicts[3])

        incidence_dict = {}
        for key, value in graph_dicts[4].items():
            incidence_dict[str(key)] = value

        prob_dict = {}
        for key, value in graph_dicts[5].items():
            prob_dict[str(key)] = value


    states = list(max(id_to_progen.values(), key=len))
    num_extant = int((args.num_sampled_per_cell_type * len(states))/args.subsample_rate)

    if args.seed is None:
        seed = np.random.choice(100000)
    
    tree, used_seed1 = utils.build_tree(num_extant, seed)

    # if args.num_sampled_per_cell_type is None:
    #     tree, used_seed2 = utils.overlay_fate_map_on_tree_random_sample_cells(edge_dict, tree, states, data["root_id"], prob_dict, total_times, seed)
    #     utils.subsample_tree_random_sample_cells(tree, states, int(100 * len(states)), id_to_progen, seed)


    # else:
    #     tree, used_seed2 = utils.overlay_fate_map_on_tree_min_cells_per_type(edge_dict, tree, states, data["root_id"], prob_dict, total_times, args.num_leaves_per_cell_type, seed)
    #     if tree is not None:
    #         utils.subsample_tree_num_cells_per_type(tree, states, args.num_leaves_per_cell_type, seed)
    
    tree, used_seed2 = utils.overlay_fate_map_on_tree_num_cells_per_type(edge_dict, tree, states, root_state, prob_dict, incidence_dict, total_times, args.min_sampled_per_cell_type, seed)
    
    if tree is not None:
        full_tree_nwk, cm = MLE.generate_cm(tree)

        tot_sample = args.num_sampled_per_cell_type * sum(list(incidence_dict.values())) 
        utils.subsample_tree_min_cells_per_type(tree, states, tot_sample, id_to_progen, incidence_dict, args.min_sampled_per_cell_type, seed)
        subsample_cm = cm.loc[tree.leaves]
        subsample_cm.to_csv(f"{args.prefix}_cm_{args.tree_ind}.txt", sep='\t', header = True)

        with open(f"{args.prefix}_tree_{args.tree_ind}.txt", "w") as f:
            f.write(to_newick(tree._CassiopeiaTree__network, record_branch_lengths = True))

        # utils.subsample_tree_num_cells_per_type(tree, states, args.num_sampled_per_cell_type, seed)
        u = utils.count_unrealizations(tree, id_to_progen)

        meta = pd.DataFrame({"cellBC": tree.leaves, 
                "cell_state": [tree.get_attribute(l, "state_labels")[0] for l in tree.leaves]})
        meta.to_csv(f"{args.prefix}_meta_{args.tree_ind}.txt", sep='\t', index = False)

        with open(f"{args.prefix}_unrealizations_{args.tree_ind}.txt", "w") as f:
            f.write(f"{u}\n")

        with open(f"{args.prefix}_usedseeds_{args.tree_ind}.txt", "w") as f:
            f.write(f"{used_seed1}\t{used_seed2}\n")

        out_newick = MLE.estimate_times(tree)
        out_file_loc = f"{args.prefix}_time_estimate_tree_{args.tree_ind}.nwk"
        with open(out_file_loc, 'w') as f:
            f.write(out_newick)

    else:
        print("Could not find sampling")
        with open(f"{args.prefix}_tree_{args.tree_ind}.txt", "w") as f:
            f.write("Could not find sampling")

if __name__ == "__main__":
    main()
