"""Algorithms to post-process the transition frequency DAG."""

import networkx as nx
import numpy as np

from collections import defaultdict


def max_path_tree_from_DAG(G, states, alternate_edge_weights=None):
    """Builds the max path tree on the DAG."

    Finds the union of all max weight paths from the root node to each supplied
    terminal state. This union produces the tree that maximizes total weight
    for a tree that covers all states. Note edge weights are counted only once
    in the union.

    Args:
        G: The supplied DAG
        states: A list of all of the cell states, with the union of these states
            being the potency/label of the root node
        alternate_edge_weights: An alternate dictionary of edge weights that
            can be used instead of the weights on the graph

    Returns:
        The max path tree

    """
    max_paths = []

    for state in states:
        paths = list(nx.all_simple_paths(G, "|".join(sorted(states)), state))
        max_val = 0
        max_path = []
        for path in paths:
            if alternate_edge_weights is not None:
                path_length = sum(
                    [
                        alternate_edge_weights[(path[i], path[i + 1])]
                        for i in range(len(path) - 1)
                    ]
                )
            else:
                path_length = sum(
                    [G[path[i]][path[i + 1]]["weight"] for i in range(len(path) - 1)]
                )
            if path_length > max_val:
                max_val = path_length
                max_path = path

        max_paths.append(max_path)

    path_G = nx.DiGraph()
    path_G.add_node("|".join(sorted(states)))

    for path in max_paths:
        for i in range(len(path) - 1):
            if (path[i], path[i + 1]) not in path_G.edges:
                if alternate_edge_weights is not None:
                    path_G.add_edge(
                        path[i],
                        path[i + 1],
                        weight=alternate_edge_weights[(path[i], path[i + 1])],
                    )
                else:
                    path_G.add_edge(
                        path[i], path[i + 1], weight=G[path[i]][path[i + 1]]["weight"]
                    )

    return path_G


def max_flow_graph_through_all_states(
    G,
    states,
    min_flow_to_keep_edge=1,
    prune_unreachable_nodes=True,
    capacity_attribute="weight",
):
    """Builds the max flow graph through all specified states simultaneously.

    Creates a dummy node T attached to all specified terminal states and then
    calculates the max flow graph from the root to T, only keeping edges that
    have more flow than a specified flow threshold. Then prunes the dummy node.

    Args:
        G: The graph to calculate the max flow graph on
        states: A list of all of the cell states, with the union of these states
            being the potency/label of the root node
        prune_unreachable_nodes: Whether or not to prune away nodes that are
            not reachable from the root and their adjacent edges
        capacity_attribute: The attribute on the graph that contains the
            capacity of each edge

    Returns:
        The max flow graph
    """

    G.add_node("T")
    for i in states:
        G.add_edge(i, "T", weight=100000000)

    flow_value, flow_dict = nx.maximum_flow(
        G, "|".join(sorted(states)), "T", capacity=capacity_attribute
    )

    flow_G = nx.DiGraph()
    flow_G.add_node("|".join(sorted(states)))

    for par in flow_dict:
        for child in flow_dict[par]:
            flow = flow_dict[par][child]
            if flow >= min_flow_to_keep_edge and par.count("|") > 0:
                spar = "|".join(sorted(par.split("|")))
                schild = "|".join(sorted(child.split("|")))
                if (spar, schild) not in flow_G.edges:
                    flow_G.add_edge(spar, schild, weight=flow)
                else:
                    curr_weight = flow_G[spar][schild]["weight"]
                    flow_G[spar][schild]["weight"] = curr_weight + flow

    if prune_unreachable_nodes:
        for n in flow_G.nodes:
            if flow_G.in_degree(n) == 0 and n != "|".join(sorted(states)):
                flow_G.remove_node(n)

    G.remove_node("T")

    return flow_G


def max_flow_graph_through_each_state(
    G,
    states,
    min_flow_to_keep_edge=1,
    prune_unreachable_nodes=True,
    capacity_attribute="weight",
):
    """Builds the max flow graph through each specified state and overlays them.

    First, constructs max flow graphs through each specified state. Then,
    overlays these flow graphs, taking only edges from each flow graph that
    have more flow than a specified flow threshold. In overlaying these flow
    graphs, the flow from shared edges across graphs is summed in the final
    graph.

    Args:
        G: The graph to calculate the max flow graph on
        states: A list of all of the cell states, with the union of these states
            being the potency/label of the root node
        prune_unreachable_nodes: Whether or not to prune away nodes that are
            not reachable from the root and their adjacent edges
        capacity_attribute: The attribute on the graph that contains the
            capacity of each edge

    Returns:
        The graph that is the overlay of max flow graphs through each state
    """

    flow_dicts = []

    for state in states:
        flow_value, flow_dict = nx.maximum_flow(
            G, "|".join(sorted(states)), state, capacity="weight"
        )
        flow_dicts.append(flow_dict)

    flow_G = nx.DiGraph()
    flow_G.add_node("|".join(sorted(states)))

    for flow_dict in flow_dicts:
        for par in flow_dict:
            for child in flow_dict[par]:
                flow = flow_dict[par][child]
                if flow >= min_flow_to_keep_edge:
                    spar = "|".join(sorted(par.split("|")))
                    schild = "|".join(sorted(child.split("|")))
                    if (spar, schild) not in flow_G.edges:
                        flow_G.add_edge(spar, schild, weight=flow)
                    else:
                        curr_weight = flow_G[spar][schild]["weight"]
                        flow_G[spar][schild]["weight"] = curr_weight + flow

    if prune_unreachable_nodes:
        for n in flow_G.nodes:
            if flow_G.in_degree(n) == 0 and n != "|".join(sorted(states)):
                flow_G.remove_node(n)

    return flow_G


def path_constrained_max_flow_tree(orig_G, states):
    """Finds the union of the paths with max flow from the root to each state.

    First, finds the path with the most flow from the root to each state. Then,
    takes the tree that is the union of these paths. In the case of ties, the
    path with the most total edge weight is taken. In taking the union, flows
    are not summed for edges traversed by more than one path.

    Then, the remainder of the originally provided graph is returned. For every
    edge in the new graph, the flow through that edge is removed from the
    original. Then, nodes not reachable by the root in the original graph are
    pruned.

    Args:
        orig_G: The original graph to find the path constrained max flow tree on
        states: A list of all of the cell states, with the union of these states
            being the potency/label of the root node

    Returns:
        The path constrained max flow tree. The original graph is modified in place

    """
    G = orig_G.copy()
    max_paths = []

    for state in states:
        paths = list(nx.all_simple_paths(G, "|".join(states), state))
        max_flow = 0
        max_frequency = 0
        max_path = []
        for path in paths:
            path_weights = [
                G[path[i]][path[i + 1]]["weight"] for i in range(len(path) - 1)
            ]
            path_flow = min(path_weights)
            path_frequency = sum(path_weights)
            if path_flow > max_flow:
                max_flow = path_flow
                max_path = path
                max_frequency = path_frequency
            elif path_flow == max_flow:
                if path_frequency > max_frequency:
                    max_flow = path_flow
                    max_path = path
                    max_frequency = path_frequency

        max_paths.append((max_path, max_flow))

    path_G = nx.DiGraph()
    path_G.add_node("|".join(sorted(states)))

    for path, flow in max_paths:
        for i in range(len(path) - 1):
            if (path[i], path[i + 1]) not in path_G.edges:
                path_G.add_edge(path[i], path[i + 1], weight=flow)

    # Remove the flow on the produced graph from the original graph
    for u, v in path_G.edges:
        G[u][v]["weight"] -= path_G[u][v]["weight"]

    to_remove = [(u, v) for (u, v) in G.edges if G[u][v]["weight"] == 0]
    G.remove_edges_from(to_remove)
    to_remove = [n for n in G.nodes if G.in_degree(n) == 0]
    to_remove.remove("|".join(sorted(states)))
    G.remove_nodes_from(to_remove)

    return path_G


def breakdown_flow_accounted(G, cells_per_state):
    """Calculates the proportion of cells accounted by flow into terminal states

    Args:
        G: The graph on which to determine how much flow is accounted
        cells_per_state: A dictionary containing the total number of cells
            for each state in the data

    Returns:
        The proportion of all cells accounted and a dictionary containing the
            proportion accounted for each cell state

    """
    prop_cells_per_state = dict(
        zip(
            cells_per_state.keys(),
            [
                G.in_degree(i, weight="weight") / cells_per_state[i]
                for i in cells_per_state
            ],
        )
    )

    return (
        sum([G.in_degree(i, weight="weight") for i in cells_per_state])
        / sum(list(cells_per_state.values())),
        prop_cells_per_state,
    )


def move_flow_between_graphs_in_paths(orig_G, new_G, cells_per_state, threshold):
    """Iteratively adds flow in paths from one graph to another.

    Until a specified proportion of the total cells are explained, finds the
    path from the root to a dummy terminal node in the original graph that
    provides the most flow. Then, adds this flow to the new graph and
    subtracts the flow of this path from the original graph. The nodes in the
    original graph must be a subset of the nodes in the new graph. At each
    iteration, nodes not reachable in the original graph are pruned.

    Can be used to continue to find the best ways to further account for flow
    in the transition frequency graph until a certain proportion of cells are
    explained.

    Args:
        orig_G: The original graph to draw flow from
        new_G: The new graph to add flow to
        cells_per_state: A dictionary containing the total number of cells
            for each state in the data
        threshold: The proportion of total cells to be explained
        dummy_terminal_node_name: The name of the terminal dummy node in the
            original graph used to find flow paths

    Returns:
        None, modifies the graphs inplace
    """
    states = sorted(list(cells_per_state.keys()))

    orig_G.add_node("T")
    for i in states:
        orig_G.add_edge(i, "T", weight=orig_G.in_degree(i, weight="weight"))

    explained = breakdown_flow_accounted(new_G, cells_per_state)[0]

    while explained < threshold:
        paths = list(nx.all_simple_paths(orig_G, "|".join(states), "T"))
        max_val = 0
        max_path = None
        for path in paths:
            path_flow = min(
                [orig_G[path[i]][path[i + 1]]["weight"] for i in range(len(path) - 1)]
            )
            if path_flow > max_val:
                max_val = path_flow
                max_path = path

        for i in range(len(max_path) - 1):
            if (max_path[i], max_path[i + 1]) not in new_G.edges:
                new_G.add_edge(max_path[i], max_path[i + 1], weight=max_val)
            else:
                new_G[max_path[i]][max_path[i + 1]]["weight"] += max_val
            orig_G[max_path[i]][max_path[i + 1]]["weight"] -= max_val

        to_remove = [(u, v) for (u, v) in orig_G.edges if orig_G[u][v]["weight"] == 0]
        orig_G.remove_edges_from(to_remove)
        to_remove = [n for n in orig_G.nodes if orig_G.in_degree(n) == 0]
        to_remove.remove("|".join(states))
        orig_G.remove_nodes_from(to_remove)

        explained = breakdown_flow_accounted(new_G, cells_per_state)[0]

    new_G.remove_node("T")
    orig_G.remove_node("T")


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

def min_abundance_normalized_unrealizations_for_set(tree, states, potency_set, state_abundances, attribute_name="state_labels"):
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

def sankoff(tree, state_transition, attribute_name="state_labels"):
    # Convert state_transition to indices for speed
    num_states = len(state_transition.index)
    states_to_ind = dict(zip(state_transition.index, range(num_states)))
    transition_costs = state_transition.values

    # Instantiate data structure to hold costs for states at nodes
    costs_at_nodes = {}
    for n in tree.nodes:
        costs_at_nodes[n] = np.ones(num_states) * np.inf

    # Bottom-Up step
    for n in tree.depth_first_traverse_nodes(postorder=True):
        if tree.is_leaf(n):
            s = tree.get_attribute(n, attribute_name)[0]
            costs_at_nodes[n][states_to_ind[s]] = 0
        else:
            for s in range(num_states):
                total_cost = 0
                for v in tree.children(n):
                    min_cost_s_ = np.inf 
                    for s_ in range(num_states):
                        cost_s_ = transition_costs[s][s_] + costs_at_nodes[v][s_]
                        if cost_s_ < min_cost_s_:
                            min_cost_s_ = cost_s_
                    total_cost += min_cost_s_
                costs_at_nodes[n][s] = total_cost

    chosen_solution = {}

    # Top-Down
    for n in tree.depth_first_traverse_nodes(postorder=False):
        if tree.is_root(n):
            chosen_solution[n] = np.argmin(costs_at_nodes[n])
        else:
            p = tree.parent(n)
            s = chosen_solution[p]
            chosen_s_ = 0
            min_cost_s_ = np.inf
            for s_ in range(num_states):
                cost_s_ = transition_costs[s][s_] + costs_at_nodes[n][s_]
                if cost_s_ < min_cost_s_:
                    chosen_s_ = s_
                    min_cost_s_ = cost_s_
            chosen_solution[n] = chosen_s_
    
    final_chosen_solution = {}
    for n in chosen_solution:
        final_chosen_solution[n] = state_transition.index[chosen_solution[n]]

    return final_chosen_solution