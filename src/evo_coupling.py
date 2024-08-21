import networkx as nx
import cassiopeia as cas
import pandas as pd
import itertools
import numpy as np

def avg_evolutionary_coupling(trees, states, leaf_state_labelings, use_branch_lengths = False):
    df_minBranchDist = pd.DataFrame(0, index = states, columns = states)
    df_numcells = pd.DataFrame(0, index = states, columns = states)
    
    for i in range(len(trees)):
        tree = trees[i]
        leaf_state_labeling = leaf_state_labelings[i]
        if use_branch_lengths:
            dists = dict(nx.all_pairs_dijkstra_path_length(tree._CassiopeiaTree__network.to_undirected(), weight = "length"))
            # diameter = nx.diameter(tree._CassiopeiaTree__network.to_undirected(), weight = "length")
        else:
            dists = dict(nx.all_pairs_shortest_path_length(tree._CassiopeiaTree__network.to_undirected()))
            # diameter = nx.diameter(tree._CassiopeiaTree__network.to_undirected())
        cell_pairs = list(itertools.combinations(range(len(tree.leaves)), 2))
        dist_mat = np.zeros([len(tree.leaves), len(tree.leaves)])

        # print(diameter)

        for index, col in cell_pairs:
            dist = dists[tree.leaves[index]][tree.leaves[col]]

            dist_mat[index, col] = dist#/diameter
            dist_mat[col, index] = dist#/diameter
            
        df_pairwiseTreeDist = pd.DataFrame(dist_mat, index = tree.leaves, columns = tree.leaves)
        np.fill_diagonal(df_pairwiseTreeDist.values, df_pairwiseTreeDist)
        
        for c1, c2 in itertools.combinations(states, 2):
            cellList1 = leaf_state_labeling[leaf_state_labeling["cell_state"] == c1].index
            cellList2 = leaf_state_labeling[leaf_state_labeling["cell_state"] == c2].index
            mean_dist = np.sum(df_pairwiseTreeDist.loc[cellList1, cellList2].sum())
            df_minBranchDist.loc[c1, c2] += mean_dist
            df_minBranchDist.loc[c2, c1] += mean_dist
            df_numcells.loc[c1, c2] += (len(cellList1) * len(cellList2))
            df_numcells.loc[c2, c1] += (len(cellList1) * len(cellList2))

    print(df_minBranchDist)
    print(df_numcells)

    df_minBranchDist = df_minBranchDist.div(df_numcells)

    print(df_minBranchDist)

    df_minBranchDist.fillna(-1, inplace = True)
    df_Dist = df_minBranchDist.loc[states, states]
    
    upgma_tree = cas.data.CassiopeiaTree()
    upgma_tree.set_dissimilarity_map(df_Dist)
    
    node_name_generator = cas.solver.solver_utilities.node_name_generator()
    upgma_solver = cas.solver.UPGMASolver()
    
    _dissimilarity_map = upgma_tree.get_dissimilarity_map().copy()
    N = _dissimilarity_map.shape[0]

    tree = nx.Graph()
    tree.add_nodes_from(_dissimilarity_map.index)

    while N > 2:

        i, j = upgma_solver.find_cherry(_dissimilarity_map.to_numpy())
        
        node_i, node_j = (
            _dissimilarity_map.index[i],
            _dissimilarity_map.index[j],
        )

        new_node_name = next(node_name_generator)
        tree.add_node(new_node_name)
        tree.add_edges_from(
            [(new_node_name, node_i), (new_node_name, node_j)]
        )

        _dissimilarity_map = upgma_solver.update_dissimilarity_map(
            _dissimilarity_map, (node_i, node_j), new_node_name
        )

        N = _dissimilarity_map.shape[0]

    tree = upgma_solver.root_tree(
        tree,
        upgma_tree.root_sample_name,
        _dissimilarity_map.index.values,
    )

    upgma_tree.populate_tree(tree)
    upgma_tree.collapse_unifurcations()
    
    return upgma_tree.get_newick()