import cassiopeia as cas
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

from collections import defaultdict

import os
import sys
cwd = os.getcwd()
sys.path.append("/".join(cwd.split("/")[:-1]) + "/scripts")
import utils as utils
import algs as algs

def gen_trees(folder_prefix):
    tree_stub = "_newick_noMutationlessEdges_Labeled.nwk"
    state_label_stub = "_metadata.txt"

    tree97 = f"{folder_prefix}input_trees/AM-DNA-097_hybrid{tree_stub}"
    tree98 = f"{folder_prefix}input_trees/AM-DNA-098_hybrid{tree_stub}"
    label97 = f"{folder_prefix}formatted_and_reduced_labels/AM-DNA-097{state_label_stub}"
    label98 = f"{folder_prefix}formatted_and_reduced_labels/AM-DNA-098{state_label_stub}"

    states = ['Endoderm', 'Endothelial', 'NMPs', 'NeuralTube', 'PCGLC', 'Somite']
    states_ = ['Endoderm', 'NMPs', 'NeuralTube', 'PCGLC', 'Somite']

    state_transform_mapping = {
        "Somite": ['Somite-1', 'Somite0', 'SomiteDermo', 'SomiteSclero'],
        "NeuralTube": ["NeuralTube1", "NeuralTube2"]
    }

    trees = []

    tree = cas.data.CassiopeiaTree(tree = tree97)
    trees.append(tree)

    tree = cas.data.CassiopeiaTree(tree = tree98)
    trees.append(tree)

    for i in range(1, 25):
        if i == 17:
            continue
        tree = cas.data.CassiopeiaTree(tree = f"{folder_prefix}input_trees/Bar{i}{tree_stub}")
        trees.append(tree)

    TLS_trees = []
    TLSCL_trees = []

    experiment_labels_df = pd.read_table(f"{folder_prefix}multiseq_barcodes.txt")
    experiment_labels_dict = dict(zip(experiment_labels_df["ID"], experiment_labels_df["TLS ID"]))

    tree = cas.data.CassiopeiaTree(tree = tree97)
    utils.label_tree_with_leaf_states(tree, label97, state_transform_mapping)
    utils.prune_unwanted_states(tree, states)
    utils.impute_states_from_children(tree)
    TLS_trees.append(tree)
    # pic.dump(tree, open(f"TLS_inputs/TLS1.pkl", "wb"))

    tree = cas.data.CassiopeiaTree(tree = tree98)
    utils.label_tree_with_leaf_states(tree, label98, state_transform_mapping)
    utils.prune_unwanted_states(tree, states)
    utils.impute_states_from_children(tree)
    TLS_trees.append(tree)
    # pic.dump(tree, open(f"TLS_inputs/TLS2.pkl", "wb"))

    for i in range(1, 25):
        if i == 17:
            continue
        tree = cas.data.CassiopeiaTree(tree = f"{folder_prefix}input_trees/Bar{i}{tree_stub}")
        utils.label_tree_with_leaf_states(tree, f"{folder_prefix}formatted_and_reduced_labels/Bar{i}{state_label_stub}", state_transform_mapping)

        if "CL" in experiment_labels_dict[f"Bar{i}"]:
            utils.prune_unwanted_states(tree, states_)
            utils.impute_states_from_children(tree)
            TLSCL_trees.append(tree)
    #         pic.dump(tree, open(f"_inputs/multi_seq{i}.pkl", "wb"))
        else:
            utils.prune_unwanted_states(tree, states)
            utils.impute_states_from_children(tree)
            TLS_trees.append(tree)
    #         pic.dump(tree, open(f"TLS_inputs/multi_seq{i}.pkl", "wb"))
    return TLS_trees, states

def annotate_trees_from_progens(trees, states, progens):
    for tree_ind, tree in enumerate(trees):
        unrealizations, potency_sets_at_nodes, ind_to_potency_set = algs.min_unrealizations_for_set(tree, states, progens)

        for n in tree.depth_first_traverse_nodes(postorder = False):
            if tree.is_leaf(n):
                continue
            if n == tree.root:
                progen_ind = [i[0] for i in sorted(potency_sets_at_nodes[n].items(), key=lambda item: item[1][1])][0]
                progen_set = ind_to_potency_set[progen_ind]
                tree.set_attribute(n, "state_labels", list(sorted(progen_set)))
                continue
            parent = tree.parent(n)
            progen_inds = [i[0] for i in sorted(potency_sets_at_nodes[n].items(), key=lambda item: item[1][1])]
            for progen_ind in progen_inds:
                progen_set = ind_to_potency_set[progen_ind]
                if progen_set.issubset(set(tree.get_attribute(parent, "state_labels"))):
                    tree.set_attribute(n, "state_labels", list(sorted(progen_set)))
                    break
            leaf_set = set([tree.get_attribute(l, "state_labels")[0] for l in tree.leaves_in_subtree(n)])
            unrealizations += len(list(sorted(progen_set))) - len(leaf_set)
            
def make_heatmap_from_progens(folder_prefix, progens, singleton_state_order, index_map):
    TLS_trees, states = gen_trees(folder_prefix)
    annotate_trees_from_progens(TLS_trees, states, progens)

    transition_counts = utils.get_transition_frequencies_weighted_by_subtree_size(TLS_trees, states)
    trans_mat = pd.DataFrame(transition_counts).T
    
    warped_ind = ["|".join(sorted(list(p))) for p in progens]
    for i in warped_ind:
        if i not in trans_mat.index:
            trans_mat.loc[i,] = np.NaN
    
    state_order = sorted(trans_mat.index, key = len, reverse = True)
    state_order_progen = [s for s in state_order if s.count("|") > 0]
    state_order = state_order_progen + singleton_state_order

    trans_mat = trans_mat.reindex(state_order, axis=1)
    trans_mat = trans_mat.loc[state_order]
    trans_mat = trans_mat.loc[:,singleton_state_order]
#     np.fill_diagonal(trans_mat.values, np.NaN)
    # for i in range(len(trans_mat.values)):
    #     for j in range(len(trans_mat.values[i,])):
    #         if j > i and np.isnan(trans_mat.values[i,j]):
    #             trans_mat.values[i,j] = 0
    trans_mat = trans_mat.fillna(0)
    trans_mat = trans_mat.drop(index = states)
#     norm_trans_mat = trans_mat.div(trans_mat.sum(axis = 1), axis = 0)
#     norm_trans_mat.fillna(0, inplace = True)

    new_index = []
    new_columns = []
    prog_index = []
    for i in trans_mat.index:
        out = [ct for ct in i.split("|")]
        prog_index.append([index_map[ct] for ct in i.split("|")])
        out = sorted(out, key = lambda x: singleton_state_order.index(x))
        new_index.append(",".join([index_map[ct] for ct in out]))
    for i in trans_mat.columns:
        new_columns.append(",".join([index_map[ct] for ct in i.split("|")]))
    trans_mat.index = new_index
    trans_mat.columns = new_columns

    return trans_mat.T, prog_index