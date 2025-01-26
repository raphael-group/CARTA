"""Utilities for building the transition frequency DAG."""

import ete3
import itertools
import numpy as np
import networkx as nx
import pandas as pd
import cassiopeia as cas
import json

from collections import defaultdict


def label_tree_with_leaf_states(
    tree,
    state_labels,
    state_transform_mapping=None,
    state_to_stateid=None,
    permute=False,
    attribute_name="state_labels",
    state_column="cell_state",
):
    """This function labels the leaves of a tree with given cell state labels.

    Args:
        tree: The tree to have its leaves labeled
        state_labels: Either the file or dataframe containing the labels for each cell
        state_transform_mapping: An optional dictionary mapping states to
            states, e.g. to consolidate substates into larger state groups
        state_to_stateid: An optional dictionary mapping states to state IDs,
            e.g. to map state names to indexes identifying those states
        permute: Whether to permute the cell state labels on the leaves
        attribute_name: The CassiopeiaTree attribute to store state labels
        state_column: The column in the dataframe containing the state info

    Returns:
        None. Modifies the tree in place

    """
    if type(state_labels) == str:
        state_labels = pd.read_csv(state_labels, sep="\t")
    
    if permute:
        state_labels[state_column] = np.random.permutation(state_labels[state_column])
    if state_transform_mapping is not None:
        state_labels = transform_states(state_labels, state_transform_mapping)

    nwklabel_to_state = dict(zip(state_labels["cellBC"], state_labels[state_column]))

    if state_to_stateid is not None:
        temp = {}
        for cell in nwklabel_to_state:
            temp[cell] = state_to_stateid[nwklabel_to_state[cell]]
        nwklabel_to_state = temp

    for node in tree.nodes:
        if tree.is_leaf(node):
            tree.set_attribute(node, attribute_name, [str(nwklabel_to_state[node])])
        else:
            tree.set_attribute(node, attribute_name, [])


def impute_states_from_children(tree, attribute_name="state_labels"):
    """Imputes the progenitor field for each internal node in the tree.

    The state progenitor field for each internal node is imputed as the set of
    states in the leaf children of a node. This is propogated by set addition
    of the fields of the children of each node.

    Args:
        tree: The tree on which to impute the states of internal nodes
        attribute_name: The CassiopeiaTree attribute storing the state
            progenitor field for each node

    Returns:
        None, modifies the tree in place

    """
    for node in tree.depth_first_traverse_nodes(postorder=True):
        if len(tree.get_attribute(node, attribute_name)) == 0:
            states_at_curr_node = set()
            for child in tree.children(node):
                states_at_child = tree.get_attribute(child, attribute_name)
                states_at_curr_node = states_at_curr_node.union(set(states_at_child))
            tree.set_attribute(node, attribute_name, sorted(list(states_at_curr_node)))


def get_transition_frequencies_over_trees(trees, states, attribute_name="state_labels"):
    """Counts the transitions between state progenitor fields (set of states).

    For each edge (u, v), adds 1 to the count f(s, s'), when u has the
    set of leaf descendants s and v has the set of leaf descendants s'.
    f(s, s') is calculated across all edges in all trees.

    The output format is a nested dictionary with state progenitor fields
    represented as a string joined by "|". An example output is given for
    this list of 3 states: [state1, state2, state3] is given as:

    {
        "state1|state2|state3": defaultdict(int,
             {"state1|state2": 5}),
        "state1|state2": defaultdict(int,
             {"state1": 4})
    }

    Args:
        trees: A list of trees containing the trees to be counted over
        states: A list of states
        attribute_name: The CassiopeiaTree attribute storing the state
            progenitor field for each node

    Returns:
        A nested dictionary storing the counts, with fields represented as
            strings

    """
    counts = {}
    for tree in trees:
        for node1, node2 in tree.edges:
            start_states_string = "|".join(
                sorted(list(map(str, tree.get_attribute(node1, attribute_name))))
            )
            end_states_string = "|".join(
                sorted(list(map(str, tree.get_attribute(node2, attribute_name))))
            )

            if start_states_string not in counts:
                counts[start_states_string] = defaultdict(int)
            counts[start_states_string][end_states_string] += 1

        # root_states_string = "|".join(sorted(list(map(str, tree.get_attribute(tree.root, attribute_name)))))
        # if root_states_string != "|".join(sorted(states)):
        #     if root_states_string not in counts:
        #         counts[root_states_string] = defaultdict(int)
        #     counts["|".join(sorted(states))][root_states_string] += len(tree.leaves)

    return counts


def get_transition_frequencies_weighted_by_subtree_size(
    trees, states, attribute_name="state_labels"
):
    """Counts the transitions between state progenitor fields, weighted by subtree size.

    For each edge (u, v), adds the number of leaf descendants to the count
    f(s, s'), when u has the set of leaf descendants s and v has the set of
    leaf descendants s'. f(s, s') is calculated across all edges in all trees.
    This weighting allows flow conservation on the resultant graph built from
    these transitions.

    If the state progenitor field of the root of a tree is not equal to the
    set of all states, an additional implicit edge is added from the set of
    all states to the root progenitor field, adding to that transition count
    the number of the leaves in the tree. This again maintains flow emerging
    from the fully potent progenitor.

    The output format is a nested dictionary with state progenitor fields
    represented as a string joined by "|". An example output is given for
    this list of 3 states: [state1, state2, state3] is given as:

    {
        "state1|state2|state3": defaultdict(int,
             {"state1|state2": 5}),
        "state1|state2": defaultdict(int,
             {"state1": 4})
    }

    Args:
        trees: A list of trees containing the trees to be counted over
        states: A list of states
        attribute_name: The CassiopeiaTree attribute storing the state
            progenitor field for each node

    Returns:
        A nested dictionary storing the counts, with fields represented as
            strings

    """
    counts = {}
    for tree in trees:
        num_leaves = {}
        for n in tree.depth_first_traverse_nodes(postorder=True):
            num_leaves[n] = len(tree.leaves_in_subtree(n))
        for node1, node2 in tree.edges:
            start_states_string = "|".join(
                sorted(list(map(str, tree.get_attribute(node1, attribute_name))))
            )
            end_states_string = "|".join(
                sorted(list(map(str, tree.get_attribute(node2, attribute_name))))
            )

            if start_states_string not in counts:
                counts[start_states_string] = defaultdict(int)
            counts[start_states_string][end_states_string] += num_leaves[node2]

        root_states_string = "|".join(
            sorted(list(map(str, tree.get_attribute(tree.root, attribute_name))))
        )
        if root_states_string != "|".join(sorted(states)):
            if "|".join(sorted(states)) not in counts:
                counts["|".join(sorted(states))] = defaultdict(int)
            counts["|".join(sorted(states))][root_states_string] += len(tree.leaves)

    return counts


def transform_states(df, state_transform_mapping, state_column="cell_state"):
    """A helper function that maps states to other states in a dataframe.

    Args:
        df: The dataframe to replace states in
        state_transform_mapping: A dictionary mapping states to states
        state_column: The column in the dataframe containing the state info

    Returns:
        df: The dataframe with states replaced
    """

    df = df.copy()
    for final_state in state_transform_mapping:
        for state_to_map in state_transform_mapping[final_state]:
            df[state_column].replace(state_to_map, final_state, inplace=True)

    return df


def prune_unwanted_states(tree, states, attribute_name="state_labels"):
    """Prunes lineages that end in leaf cells with unwanted state labels.

    Args:
        tree: The tree to prune unwanted states from
        states: A list of states to keep. All other states will be discarded
        attribute_name: The CassiopeiaTree attribute that stores the cell state
            labels on the tree

    Returns:
        None, modifies tree in-place

    """

    to_remove = []
    for leaf in tree.leaves:
        state = tree.get_attribute(leaf, attribute_name)[0]
        if state not in states:
            to_remove.append(leaf)
    tree.remove_leaves_and_prune_lineages(to_remove)
    tree.collapse_unifurcations()


def generate_hasse_diagram(lst):
    """Takes a list of all states in the state space and returns the Hasse Diagram on those states.

    Args:
        lst: A list containing the space of extant states

    Returns:
        The Hasse Diagram on the states, with all edge weights initialized as 0

    """

    lst = sorted(lst)

    G = nx.DiGraph()
    G.add_node("|".join(lst))

    node_queue = [lst]

    while len(node_queue) != 0:
        curr_node = node_queue.pop(0)

        for child in list(itertools.combinations(curr_node, len(curr_node) - 1)):
            child = list(child)
            G.add_node("|".join(child))
            G.add_edge("|".join(curr_node), "|".join(child), weight=0)
            if len(child) > 1:
                node_queue.append(child)
    return G


def annotate_hasse_with_transitions(
    G, states, transition_counts, remove_zero_weight_edges=False
):
    """Annotates the Hasse Diagram on the states with the frequency counts.

    Adds the frequency counts to each edge weight for each direct edge
    in the Hasse Diagram, and adds observed transitive edges to the
    Hasse Diagram weighted by their frequency count.

    Args:
        G: The unannotated Hasse Diagram
        transition_counts: A nested dictionary, where the first layer
            contains the daughter state node for each parent state node,
            and the second layer contains the frequency count for each edge
        remove_zero_weight_edges: Whether to remove zero weight edges from
            the resultant graph

    Returns:
        G: The annotated Hasse Diagram

    """
    for parent in transition_counts:
        sorted_parent = "|".join(sorted(parent.split("|")))
        for child in transition_counts[parent]:
            sorted_child = "|".join(sorted(child.split("|")))
            if sorted_parent == sorted_child:
                continue
            if G.has_edge(sorted_parent, sorted_child):
                G[sorted_parent][sorted_child]["weight"] = transition_counts[parent][
                    child
                ]
            else:
                G.add_edge(
                    sorted_parent, sorted_child, weight=transition_counts[parent][child]
                )

    if remove_zero_weight_edges:
        to_remove = [(u, v) for (u, v) in G.edges if G[u][v]["weight"] == 0]
        G.remove_edges_from(to_remove)
        to_remove = [n for n in G.nodes if G.in_degree(n) == 0]
        to_remove.remove("|".join(sorted(states)))
        G.remove_nodes_from(to_remove)

    return G

def generate_perfect_phylogeny(df_binary):

    solT_mut = nx.DiGraph()
    solT_mut.add_node('root')

    solT_cell = nx.DiGraph()
    solT_cell.add_node('root')

    df_binary = df_binary[df_binary.sum().sort_values(ascending=False).index]    

    for cell_id, row in df_binary.iterrows():
        if cell_id == 'root':
            continue

        curr_node = 'root'
        for column in df_binary.columns[row.values == 1]:
            if column in solT_mut[curr_node]:
                curr_node = column
            else:
                if column in solT_mut.nodes:
                    raise NameError(f'{column} is being repeated')
                solT_mut.add_edge(curr_node, column)
                solT_cell.add_edge(curr_node, column)
                curr_node = column

        solT_cell.add_edge(curr_node, cell_id)   

    return solT_mut, solT_cell

def tree_to_newick(T, root=None):
    if root is None:
        roots = list(filter(lambda p: p[1] == 0, T.in_degree()))
        assert 1 == len(roots)
        root = roots[0][0]
    subgs = []
    while len(T[root]) == 1:
        root = list(T[root])[0]
    for child in T[root]:
        while len(T[child]) == 1:
            child = list(T[child])[0]
        if len(T[child]) > 0:
            child_newick = tree_to_newick(T, root=child)
            if child_newick != '()':
                subgs.append(child_newick)
        else:
            subgs.append(child)
    if len(subgs) == 1:
        return str(subgs[0])
    else:
        return "(" + ','.join(map(str, subgs)) + ");"


def score_rf(T1, T2):

    T1 = ete3.Tree(T1, format=1)
    T2 = ete3.Tree(T2, format=1)

    (
        rf,
        rf_max,
        names,
        edges_t1,
        edges_t2,
        discarded_edges_t1,
        discarded_edges_t2,
    ) = T1.robinson_foulds(T2, unrooted_trees=True)

    return rf, rf_max

def get_progens_from_tree_input(file):
    tree = cas.data.CassiopeiaTree(tree = file)
    for l in tree.leaves:
        tree.set_attribute(l, "state_labels", [l])
    for n in tree.internal_nodes:
        tree.set_attribute(n, "state_labels", [])
    impute_states_from_children(tree)

    return [tree.get_attribute(n, "state_labels") for n in tree.internal_nodes]

def get_progens_from_ilp_output(file):
    progens = []
    with open(file, 'r') as inp:
        for line in inp:
            progens.append([x.lstrip("'").rstrip("'") for x in line.rstrip('}\n').lstrip('{').split(', ')])

    return progens

def jaccard_dist(s1, s2, states):
    s1 = set([state_set_to_int(s, states) for s in s1])
    s2 = set([state_set_to_int(s, states) for s in s2])

    return 1 - len(s1.intersection(s2))/len(s1.union(s2))

def state_set_to_int(state_set, states):
    tot = 0
    for i in range(len(states)):
        if states[i] in state_set:
            tot += 2**i
    return tot

def sum_of_min_symmetric_difference(s1, s2):
    tot_dist = 0
    for s in s1:
        s = set(s)
        min_dist = np.inf
        for s_ in s2:
            s_ = set(s_)
            dist = len(s.union(s_) - s.intersection(s_))
            if dist < min_dist:
                min_dist = dist
        tot_dist += min_dist
    return tot_dist

def classification_scores(ground, recon, states):
    ground = set([state_set_to_int(s, states) for s in ground])
    recon = set([state_set_to_int(s, states) for s in recon])

    true_positives = len(recon.intersection(ground))

    possible_progens = 2**len(states) - len(states) - 1

    precision = true_positives/len(recon)
    recall = true_positives/len(ground)
    false_positive_rate = len(recon - ground)/(possible_progens - len(ground))

    return precision, recall, false_positive_rate

def get_progens_qfm_json_format(fatemap_fname):
    with open(fatemap_fname) as f:
        fatemap = json.load(f)

    fatemap_nxtree = nx.DiGraph()
    for edge in fatemap['edges']:
        fatemap_nxtree.add_edge(edge['in_node'], edge['out_node'])
    
    for node in fatemap_nxtree:
        if fatemap_nxtree.in_degree(node) == 0:
            fatemap_root = node
    
    progenitors = {}
    for node in nx.dfs_postorder_nodes(fatemap_nxtree, fatemap_root):
        if len(fatemap_nxtree[node]) == 0:
            progenitors[node] = set([node])
        else:
            progenitors[node] = set()
            for child in fatemap_nxtree[node]:
                progenitors[node] = progenitors[node] | progenitors[child]

    progens = [i for i in list(progenitors.values()) if len(i) > 1]
    new_progens = []
    for i in progens:
        new_progens.append(["type_" + s for s in i])

    return new_progens

def min_unrealizations_for_set(tree, states, potency_set, attribute_name="state_labels"):
    """Given a progeny set, returns the minimum unrealizations on a tree.

    Implements a DP algorithm that finds the minimum number of unrealizations 
    on a tree. The singleton state progeny sets (sets of size 1) and the state
    progeny that is the full union of all states are implictly included in the
    progency set.

    Note: currently the root node is implicitly labeled as the state progeny
    that is the full union of all states, BUT the number of unrealizations at
    the root is always counted as 0.

    Args:
        tree: The tree over which to calculate unrealizations
        states: A list of all of the cell states, with the union of these states
            being the potency/label of the root node
        potency_set: A set of possible progenitors, which are sets of states
        attribute_name: The CassiopeiaTree attribute to store state labels


    Returns:
        A tuple containing 1) the minimum realizations 2) A dictionary of state
        progenitor indexes to their unrealization scores 3) A dictionary 
        mapping state progenitor indexes to the progenitor set
    """

    potency_set = [set(i) for i in potency_set]
    if set(states) in potency_set:
        potency_set.remove(set(states))
    potency_set.insert(0, set(states))
    for i in states:
        if set([i]) not in potency_set:
            potency_set.append(set([i]))

    ind_to_potency_set = dict(zip(range(len(potency_set)), potency_set))

    leaf_set_at_nodes = {}

    for n in tree.nodes:
        leaf_set_at_nodes[n] = set(
            [
                tree.get_attribute(l, attribute_name)[0]
                for l in tree.leaves_in_subtree(n)
            ]
        )

    potency_sets_at_nodes = defaultdict(defaultdict)
    potency_sets_at_nodes[tree.root][0] = (0, 0)

    for n in tree.depth_first_traverse_nodes(postorder=False):
        if n == tree.root or n in tree.leaves:
            continue

        leaf_set = leaf_set_at_nodes[n]

        for ind, potency_set in ind_to_potency_set.items():
            if leaf_set.issubset(potency_set):
                potency_sets_at_nodes[n][ind] = (len(potency_set) - len(leaf_set), 0)

    for n in tree.depth_first_traverse_nodes(postorder=True):
        if tree.is_leaf(n):
            continue

        for ind, scores in potency_sets_at_nodes[n].items():

            total_at_this_node = scores[0]

            for child in tree.children(n):
                if tree.is_leaf(child):
                    continue

                best_child_state_score = np.inf

                for c_ind, c_scores in potency_sets_at_nodes[child].items():

                    if ind_to_potency_set[c_ind].issubset(ind_to_potency_set[ind]):

                        if c_scores[1] < best_child_state_score:

                            best_child_state_score = c_scores[1]

                total_at_this_node += best_child_state_score

            potency_sets_at_nodes[n][ind] = (scores[0], total_at_this_node)

    return (potency_sets_at_nodes[tree.root][0][1], potency_sets_at_nodes, ind_to_potency_set)