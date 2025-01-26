import argparse
import json
import pandas as pd
import numpy as np
import networkx as nx
import gurobipy as gp
import time
 
from gurobipy import GRB
from itertools import chain, combinations
from loguru import logger
from Bio import Phylo
from tqdm import tqdm

import fastCARTA as fc

def parse_args():
    parser = argparse.ArgumentParser(description='fastCARTA')
    parser.add_argument('-t', '--tree', type=str, required=True, help='tree file path (newick format)')
    parser.add_argument('-o', '--output', type=str, required=True, help='output files prefix')
    parser.add_argument('-l', '--label', type=str, help='leaf labeling file path (tsv format)', default=None)
    parser.add_argument('-ks','--klist', nargs='*', type = int, help='<Required> Set flag',default=None)
    parser.add_argument('--no_enforce_tree', action='store_false')
    parser.add_argument("--time_limit_sec", type=int, default=43200, help="time limit on ilp runs")
    parser.add_argument('-m', '--mode', help='mode for solver', choices=['exact', 'round'], default='exact')
    parser.add_argument('--radius', type=int, help='radius for progenitor set construction', default=0)
    parser.add_argument('--prune-states', type=str, nargs='+', help='cell states to prune from the tree', default=[None])
    parser.add_argument('--graph_file', type=str, default=None)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()

    if not args.no_enforce_tree:
        logger.info(f'Not Enforcing Tree')
        print(f'Not Enforcing Tree')
    else:
        logger.info(f'Enforcing Tree')
        print(f'Enforcing Tree')

    T = fc.from_newick_get_nx_tree(args.tree)
    leaves = [u for u in T.nodes if T.out_degree(u) == 0]
    
    if args.label is None:
        L = pd.DataFrame([leaves, ["type_" + i.split("_")[1] for i in leaves]]).T
        L.columns = ["cellBC", "cell_state"]
    else:
        L = pd.read_csv(args.label, sep='\t')

    labeling = dict(zip(L['cellBC'], L['cell_state'])) # map from cell name -> cell state
    descendants = dict((u, set([u]) | nx.descendants(T, u)) for u in T.nodes) # map from node -> set of descendants (including itself)

    prune_states = set(args.prune_states)
    prune_leaves = [u for u in leaves if labeling[u] in prune_states]

    logger.info(f'Dropping {len(prune_leaves)} leaves labeled by states {prune_states}...')
    for u in prune_leaves:
        del labeling[u]
        while T.out_degree(u) == 0 and T.in_degree(u) == 1:
            v = list(T.predecessors(u))[0]
            T.remove_node(u)
            u = v

    logger.info('Building CDMP model...')

    leaves = [u for u in T.nodes if T.out_degree(u) == 0]
    cell_states = set(labeling[u] for u in leaves)

    if len(cell_states) < 8:
        args.radius = -1
        print("Exhaustive radius")

    if args.klist is None:
        k_list = [len(cell_states) - 1]
    else:
        k_list = args.klist

    if args.graph_file is not None:
        if args.graph_file == "16_flag":
            base_k = 15
        else:
            graph_dicts = []
            with open(args.graph_file) as g_file:
                for line in g_file:
                    graph_dicts.append(json.loads(line))

            edge_dict = {}
            for key, value in graph_dicts[0].items():
                edge_dict[str(key)] = value

            base_k = len([i for i in edge_dict.values() if len(i) > 1])
        k_list = list(range(max(1, base_k - 3), base_k + 4))

    start_time = time.time()
    model, cost_dict, progenitors, progenitor_ids, progenitor_variables, decision_variables = fc.build_cdmp_model(T, labeling, descendants, enforce_tree=args.no_enforce_tree, radius=args.radius, mode=args.mode)
    end_time = time.time()
    build_time = end_time - start_time
    print(build_time)


    for k in k_list:
        logger.info(f'Solving for k = {k} ...')

        model.addConstr(
            gp.quicksum(progenitor_variables[p] for p in progenitor_ids) <= k + len(cell_states),
            name='total_progenitors'
        )

        logger.info('Solving model using Gurobi...')
        start_time = time.time()
        model.setParam('TimeLimit', args.time_limit_sec)
        model.optimize()
        end_time = time.time()
        solve_time = end_time - start_time

        logger.info('Writing output files...')

        decision_variables_t = {k: v.X for k, v in decision_variables.items()}
        with open(f'{args.output}_edge_labeling_{k}.txt', 'w') as f:
            f.write('src,dst,src_progenitor,dst_progenitor\n')
            for (u, v, p1, p2), val in decision_variables_t.items():
                if val > 0:
                    f.write(f'{u},{v},"{progenitors[p1]}","{progenitors[p2]}"\n')

        with open(f'{args.output}_progenitors_{k}.txt', 'w') as f:
            for i, p in enumerate(progenitors):
                if progenitor_variables[i].X > 0:
                    f.write(f'{set(p)}\n')

        with open(f'{args.output}_results_{k}.txt', 'w') as f:
            f.write('objective_score\tmodel_runtime\tmodel_buildtime\n')
            f.write(f'{model.ObjVal}\t{solve_time}\t{build_time}\n')

        c = model.getConstrByName('total_progenitors')
        model.remove(c)
        model.update()

    with open(f'{args.output}_dummy.txt', 'w') as f:
        f.write('Done!\n')

