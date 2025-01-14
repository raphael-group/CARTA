import argparse
import json
import pandas as pd
import numpy as np
import networkx as nx
import gurobipy as gp
import time
 
from collections import defaultdict
from gurobipy import GRB
from itertools import chain, combinations
from loguru import logger
from Bio import Phylo
from tqdm import tqdm

def powerset(iterable):
    """
    Computes the powerset of a given iterable.
    """
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def check_laminarity(p1, p2):
    """
    Checks if two progenitors are laminar.
    """
    is_laminar = False
    if len(p1.intersection(p2)) == 0:
        is_laminar = True
    if p1.issubset(p2):
        is_laminar = True
    if p2.issubset(p1):
        is_laminar = True
    return is_laminar

def construct_progenitors(cell_states, descendant_labels, radius=-1):
    """
    Constructs the set of candidate progenitors for the CDMP model,
    if the radius is -1, then all possible progenitors are considered.
    Otherwise, the set of progenitors is constructed by taking all
    progenitors within a radius of the `tight` progenitor set.
    """
    if radius == -1:
        progenitors = list(powerset(cell_states))
        progenitor_ids = list(range(len(progenitors)))
        return progenitors, progenitor_ids

    # counter = 0
    # progenitors = []
    # for _, progenitor in descendant_labels.items():
    #     progenitors.append(frozenset(progenitor))
    #     for i in range(1, radius + 1):
    #         for subset in combinations(progenitor, i):
    #             progenitors.append(frozenset(progenitor - set(subset)))

    progenitors = []
    for _, progenitor in descendant_labels.items():
        progenitors.append(frozenset(progenitor))
        if len(progenitor) > 2:
            for i in range(1, radius + 1):
                for subset in combinations(set(cell_states) - progenitor, i):
                    progenitors.append(frozenset(progenitor.union(set(subset))))

    progenitors = list(set(progenitors))
    progenitors = [tuple(p) for p in progenitors]

    progenitor_ids = list(range(len(progenitors)))
    return progenitors, progenitor_ids

def build_cdmp_model(T, labeling, descendants, enforce_tree = False, root='root', radius=-1, mode='exact'):
    """
    Builds the cell differentiation mapping problem (CDMP) model
    for the given tree T and labeling using the provided
    parameters. If mode is 'exact', then the model contains O(n)
    binary variables, otherwise, the model contains no integer 
    variables, and a rounding procedure is required.
    """

    logger.info('Constructing LP model for cell differentiation mapping problem (CDMP)...')
 
    leaves = [u for u in T.nodes if T.out_degree(u) == 0]
    cell_states = set(labeling[u] for u in leaves)
    descendant_labels = dict((u, set(labeling[v] for v in descendants[u] if T.out_degree(v) == 0)) for u in T.nodes)

    # add a dummy root node to enforce agreement at the root
    T.add_node("dummy_root")
    T.add_edge("dummy_root", root)

    descendant_labels['dummy_root'] = descendant_labels[root]

    # valid transitions are pairs of progenitors (p1, p2) such that p2 is a subset of p1
    progenitors, progenitor_ids = construct_progenitors(cell_states, descendant_labels, radius=radius)
    valid_transitions = set()
    valid_transitions_edges = defaultdict(set)
    children_transitions = defaultdict(list)
    parent_transitions = defaultdict(list)

    for p1 in progenitor_ids:
        for p2 in progenitor_ids:
            prog1 = progenitors[p1]
            prog2 = progenitors[p2]
            prog1_set = set(prog1)
            prog2_set = set(prog2)

            if not prog1_set.issubset(prog2_set):
                continue 

            for e in T.edges:
                descendants_1 = descendant_labels[e[1]]
                if not descendants_1.issubset(prog1_set):
                    continue 

                descendants_2 = descendant_labels[e[0]]
                if not descendants_2.issubset(prog2_set):
                    continue

                valid_transitions.add((e[0], e[1], p2, p1))
                valid_transitions_edges[e].add((p2, p1))
                children_transitions[(e, p2)].append(p1)
                parent_transitions[(e, p1)].append(p2)

    logger.info(f'{len(cell_states)} cell states: {cell_states}')
    logger.info(f'{len(progenitors)} candidate progenitors.')
    logger.info(f'Number of cells: {len(leaves)}')
    logger.info(f'Number of decision variables: {len(valid_transitions) + len(progenitors)}')

    model = gp.Model()

    """
    The variable x_{e,p1,p2} = 1 iff edge e is labeled by progenitor p1 and p2, and is
    x_{e,p1,p2} = 0 otherwise.
    """
    x = model.addVars(valid_transitions, vtype=GRB.CONTINUOUS, lb=0.0, name='x')

    logger.info("Adding edge agreement constraints...")

    """ Require \sum_{p'} x_{u,v,p,p'} = \sum_{p'}x_{v,w,p',p} for all (u,v),(v,w),p"""
    for idx, (u, v) in enumerate(T.edges):
        if idx % 50 == 0:
            logger.info(f'Added edge agreement constraints for {idx}/{len(T.edges)} edges...')

        for w in T[v]:
            for p in progenitor_ids:
                model.addConstr(
                    gp.quicksum(x[u, v, p1, p] for p1 in parent_transitions[((u,v),p)]) == gp.quicksum(x[v, w, p, p1] for p1 in children_transitions[((v,w),p)]),
                    name=f'conservation_{u}_{v}_{w}_{p}'
                )

    for u, v in T.edges:
        model.addConstr(
            gp.quicksum(x[u, v, p1, p2] for (p1, p2) in valid_transitions_edges[(u,v)]) == 1,
            name=f'sum_to_one_{u}_{v}'
        )

    """
    Set \sum_{c} x_{u,v,p,p'} = 0 for all e=(u,v), v is a leaf, and p' is not the label of v
    """
    for v in leaves:
        cell_state = labeling[v]
        u = list(T.predecessors(v))[0]
        for progenitor_id in progenitor_ids:
            p = progenitors[progenitor_id]
            if len(p) == 1 and p[0] == cell_state:
                continue

            e = (u,v)
            model.addConstr(
                gp.quicksum(x[u, v, p1, progenitor_id] for p1 in parent_transitions[(e, progenitor_id)]) == 0,
                name=f'leaves_{u}_{v}_{progenitor_id}'
            )

    logger.info("Adding linking constraints...")
    if mode == 'exact':
        b = model.addVars(progenitor_ids, vtype=GRB.BINARY, name='b')
    else:
        b = model.addVars(progenitor_ids, vtype=GRB.CONTINUOUS, lb=0.0, name='b')

    for idx, (u, v, p1, p2) in enumerate(valid_transitions):
        if idx % 1000 == 0:
            logger.info(f'Added {idx}/{len(valid_transitions)} linking constraints...')

        model.addConstr(
            x[u, v, p1, p2] <= b[p1],
            name=f'linking_{u}_{v}_{p1}_{p2}_1'
        )

        model.addConstr(
            x[u, v, p1, p2] <= b[p2],
            name=f'linking_{u}_{v}_{p1}_{p2}_2'
        )

    if enforce_tree: 
        logger.info("Enforcing Tree: adding tree constraints...")
        for p1, p2 in combinations(progenitor_ids, 2):
            if not check_laminarity(set(progenitors[p1]), set(progenitors[p2])):
                model.addConstr(
                    b[p1] + b[p2] <= 1,
                    name=f'tree_{p1}_{p2}_not_laminar'
            )

    for p in progenitor_ids:
        if len(progenitors[p]) == 1:
            b[p].lb = 1.0
            b[p].ub = 1.0

    cost_dict = {}
    for u, v, p1, p2 in valid_transitions:
        desc = descendant_labels[v]
        prog2 = progenitors[p2]
        cost_dict[(u, v, p1, p2)] = max(0, len(prog2) - len(desc))

    model.setObjective(
        gp.quicksum(cost_dict[(u, v, p1, p2)] * x[u, v, p1, p2]
                    for u, v, p1, p2 in valid_transitions),
        GRB.MINIMIZE
    )

    return model, cost_dict, progenitors, progenitor_ids, b, x

"""
Loads a rooted tree from a newick file and returns 
a directed tree in networkx format, the root of which 
is labeled 'root'.
"""
def from_newick_get_nx_tree(tree_path):
    phylo_tree = Phylo.read(tree_path, 'newick')
    net_tree = Phylo.to_networkx(phylo_tree)

    node_renaming_mapping = {}
    idx = 0
    for node in net_tree.nodes:
        if str(node) == 'Clade':
            node_renaming_mapping[node] = f'clade_{idx}'
            idx = idx + 1
        else:
            node_renaming_mapping[node] = node.name
    node_renaming_mapping[list(net_tree.nodes)[0]] = 'root'
    
    net_tree = nx.relabel_nodes(net_tree, node_renaming_mapping)

    directed_tree = nx.DiGraph()
    directed_tree.add_edges_from(list(nx.bfs_edges(net_tree, 'root')))
    return directed_tree