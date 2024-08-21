from itertools import chain, combinations
import numpy as np
import networkx as nx
from collections import defaultdict
import cassiopeia as cas
from cassiopeia.simulator.TreeSimulator import TreeSimulatorError

import utils

def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def build_graph(states, k, root_time, seed):
    
    np.random.seed(seed)
    
    poset = [set(i) for i in powerset(states) if len(i) > 1 and len(i) < len(states)]
    progens = list(np.random.choice(poset, size = k, replace = False))
    
    id_to_progen = dict(zip(range(len(states) + 1, len(progens) + len(states) + 2), (list(p) for p in progens)))
    all_nodes = [set([i]) for i in states] + [set(states)] + progens
    G = nx.DiGraph()
    
    for ind, p in enumerate(all_nodes):
        G.add_node(ind)
    for ind, p in enumerate(all_nodes):
        for ind_, p_ in enumerate(all_nodes):
            if ind_ <= ind:
                continue
            if p.issubset(p_):
                G.add_edge(ind_, ind)
            if p_.issubset(p):
                G.add_edge(ind, ind_)
    G = nx.transitive_reduction(G)
    G.add_node('root')
    G.add_edge('root', len(states), length = root_time)
    
    paths_from_root = []
    for i in range(len(states)):
        paths_from_root += list(nx.all_simple_paths(G, len(states), i))
    paths_from_root = sorted(paths_from_root, key = len)

    total_times = defaultdict(int)
    total_times[len(states)] = root_time

    for p in paths_from_root:
        ind = 1
        while p[-ind] not in total_times:
            ind += 1
        total_remaining_time = 1 - total_times[p[-ind]]
        for i in range(ind - 2):
            curr_ind = len(p) - ind + i
            total_times[p[curr_ind + 1]] = total_times[p[-ind]] + (i + 1) * total_remaining_time/(ind - 1)

    for i in range(len(states)):
        total_times[i] = 1.0
        id_to_progen[i] = states[i]

    for p in paths_from_root:
        for i in range(len(p) - 1):
            G[p[i]][p[i + 1]]["length"] = total_times[p[i + 1]] - total_times[p[i]]
            
    id_to_progen[len(states)] = states
    
    edge_dict = defaultdict(list)
    for u, v in G.edges:
        edge_dict[u].append(v)
    for i in range(len(states)):
        edge_dict[i] = []
        
    return dict(edge_dict), total_times, id_to_progen

def build_tree_graph(states, k, root_time, seed):

    np.random.seed(seed)

    G = nx.DiGraph()

    G.add_node(len(states))

    flex = len(states) - k - 2
    node_name_gen = len(states)
    internal = []
    nodes = [len(states)]

    while len(internal) < k + 1 or len(nodes) < len(states) + k + 1:
        curr_node = np.random.choice(nodes)
        if len(list(G.successors(curr_node))) == 0:
            if len(internal) >= k + 1:
                continue
            for _ in range(2):
                node_name_gen += 1
                G.add_edge(curr_node, node_name_gen)
                nodes.append(node_name_gen)
            internal.append(curr_node)
        else:
            if flex == 0:
                continue
            flex -= 1
            node_name_gen += 1
            G.add_edge(curr_node, node_name_gen)
            nodes.append(node_name_gen)

    G.add_node('root')
    G.add_edge('root', len(states), length = root_time)

    leaves = [l for l in G.nodes if G.out_degree(l) == 0]
    node_mapping = dict(zip(leaves, [int(l) for l in np.random.permutation(range(len(states)))]))

    G = nx.relabel_nodes(G, node_mapping)

    paths_from_root = []
    for i in range(len(states)):
        paths_from_root += list(nx.all_simple_paths(G, len(states), i))
    paths_from_root = sorted(paths_from_root, key = len, reverse = True)

    total_times = defaultdict(int)
    total_times[len(states)] = root_time
    id_to_progen = {}
    id_to_progen[len(states)] = states

    for p in paths_from_root:
        ind = 1
        while p[-ind] not in total_times:
            ind += 1
        total_remaining_time = 1 - total_times[p[-ind]]
        for i in range(ind - 2):
            curr_ind = len(p) - ind + i
            total_times[p[curr_ind + 1]] = total_times[p[-ind]] + (i + 1) * total_remaining_time/(ind - 1)

    for p in paths_from_root:
        for i in range(len(p) - 1):
            G[p[i]][p[i + 1]]["length"] = total_times[p[i + 1]] - total_times[p[i]]

    for i in range(len(states)):
        total_times[i] = 1.0
        id_to_progen[i] = states[i]

    tree = cas.data.CassiopeiaTree(tree = G)

    for l in tree.leaves:
        tree.set_attribute(l, "state_labels", [states[int(l)]])
    for n in tree.internal_nodes:
        tree.set_attribute(n, "state_labels", [])
    utils.impute_states_from_children(tree)
    for n in tree.internal_nodes:
        if n == tree.root:
            continue
        id_to_progen[int(n)] = tree.get_attribute(n, "state_labels")

    edge_dict = defaultdict(list)
    for u, v in G.edges:
        edge_dict[u].append(v)
    for i in range(len(states)):
        edge_dict[i] = []

    return dict(edge_dict), total_times, id_to_progen

def mle_like_estimator_of_net_diversification_rate(num_extant, relative_extinction_rate, c):
    t = 1 - c * np.log2(num_extant) 
    r = np.log(num_extant)/t
    birth_rate = r/(1 - relative_extinction_rate)
    death_rate = relative_extinction_rate * birth_rate
    
    return birth_rate, death_rate

def build_tree(num_extant, seed):
    c = 0.01
    
    birth_rate, death_rate = mle_like_estimator_of_net_diversification_rate(num_extant, 0.1, c)
    
    initial_birth_scale = 1/birth_rate
    birth_waiting_distribution = lambda scale: np.random.exponential(scale) + c
    death_waiting_distribution = lambda: np.random.exponential(1/death_rate + c)

    size = 0
    while size == 0:
        try:
            bd_sim = cas.simulator.BirthDeathFitnessSimulator(
                birth_waiting_distribution = birth_waiting_distribution,
                initial_birth_scale = initial_birth_scale,
                death_waiting_distribution = death_waiting_distribution,
                num_extant = num_extant,
                random_seed = seed
            )
            tree = bd_sim.simulate_tree()
            size = tree.n_cell
        except TreeSimulatorError:
            size = 0
            np.random.seed(seed)
            seed = np.random.choice(10000)

    tree.relabel_nodes(dict(zip(tree.leaves, ["c" + i for i in tree.leaves])))

    total_time = 0
    for l in tree.leaves:
        total_time = tree.get_time(l)

    branch_length_dict = {}
    for e1, e2 in tree.edges:
        branch_length_dict[(e1, e2)] = tree.get_branch_length(e1, e2)/total_time
    tree.set_branch_lengths(branch_length_dict)
    
    return tree, seed

def overlay_fate_map_on_tree(edge_dict, tree, states, total_times, num_leaves_per_cell_type, seed):

    tree.set_attribute(tree.root, "progen_label", len(states))

    not_enough_to_subsample_flag = True
    while not_enough_to_subsample_flag:

        np.random.seed(seed)

        for i in tree.depth_first_traverse_nodes(postorder = False):
            if i == tree.root:
                continue
            time = tree.get_time(i)
            parent_progen = tree.get_attribute(tree.parent(i), "progen_label")
            possible_progens = [n for n in edge_dict[parent_progen] if total_times[parent_progen] < time]
            if len(possible_progens) == 0:
                chosen_state = parent_progen
            else:
                chosen_state = np.random.choice(possible_progens)
            tree.set_attribute(i, "progen_label", chosen_state)
            
        leaf_types = defaultdict(list)
        for l in tree.leaves:
            leaf_types[tree.get_attribute(l, "progen_label")].append(l)
        not_enough_to_subsample_flag = False
        for i in range(len(states)):
            if len(leaf_types[i]) < num_leaves_per_cell_type:
                not_enough_to_subsample_flag = True
        if not_enough_to_subsample_flag:
            np.random.seed(seed)
            seed = np.random.choice(10000)

    return tree, seed

def subsample_tree(tree, states, num_leaves_per_cell_type, seed):
    
    np.random.seed(seed)
    
    leaf_types = defaultdict(list)

    for l in tree.leaves:
        leaf_types[tree.get_attribute(l, "progen_label")].append(l)

    to_keep = []

    for i in range(len(states)):
        to_keep.extend(np.random.choice(leaf_types[i], size = num_leaves_per_cell_type, replace = False))

    to_remove = list(set(tree.leaves) - set(to_keep))

    tree.remove_leaves_and_prune_lineages(to_remove)
    tree.collapse_unifurcations()
    
def count_unrealizations(tree, id_to_progen):
    for l in tree.leaves:
        tree.set_attribute(l, "state_labels", [id_to_progen[tree.get_attribute(l, "progen_label")]])
    for n in tree.internal_nodes:
        tree.set_attribute(n, "state_labels", [])
    utils.impute_states_from_children(tree)

    num_unrealizations = 0
    for n in tree.internal_nodes:
        min_label = len(tree.get_attribute(n, "state_labels"))
        true_label = len(id_to_progen[tree.get_attribute(n, "progen_label")])
        num_unrealizations += true_label - min_label
    
    return num_unrealizations
