from argparse import ArgumentParser

import sys
import os
import cassiopeia as cas
import time
import gurobipy as gp
from loguru import logger
import pandas as pd
import networkx as nx

src = os.path.join(os.path.dirname(os.path.abspath(__file__)))
if not os.path.isdir(src):
    raise ValueError("carta directory cannot be found i.e. {}, is anything been moved?".format(src))
sys.path.append(src)

import utils
import ilp
import fastCarta as fc

def get_options():
    parser = ArgumentParser()
    parser.add_argument("--prefix", type=str, help="filepath prefix for this run")
    parser.add_argument("-ks", "--klist", nargs='*', type=int, help="list of number of progenitors", default=None)
    parser.add_argument("--file_locations", type=str, help="a file that on each line specifies the filepath of a newick file and a corresponding state metadata file in that order, seperated by a tab")
    parser.add_argument("--states_file", type=str, help="a file containing the cell states present in this dataset")
    parser.add_argument("--normalize_method", default="no_normalization", choices=["no_normalization", "cell_proportion_before_pruning", "cell_proportion_after_pruning"], help="the scheme by which to normalize the objective function",)
    parser.add_argument("--time_limit_sec", default=21600, help="time limit on ilp runs")
    parser.add_argument("--enforce_tree", default=False, action='store_true', help="whether or not to enforce that the progenitors can create a tree")
    parser.add_argument('-m', '--mode', help='mode for heuristic progenitor set solver', choices=['exact', 'round'], default='exact')
    parser.add_argument('--progen_matrix', type=str, default=None)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    return args

def main(args):
    states = []
    with open(args.states_file, "r") as f:
        for line in f:
            states.append(line.rstrip())

    if args.progen_matrix is not None:
        trees = []
        with open(args.file_locations) as file_locations:
            for line in file_locations:
                nw, metadata = line.rstrip().split("\t")
                trees.append(nw)

        if len(trees) > 1:
            raise Exception("Carta with progenitor set input only supports one cell lineage tree at this time.")

        T = fc.from_newick_get_nx_tree(trees[0])
        leaves = [u for u in T.nodes if T.out_degree(u) == 0]
        L = pd.read_csv(metadata, sep='\t')

        labeling = dict(zip(L['cellBC'], L['cell_state'])) # map from cell name -> cell state
        descendants = dict((u, set([u]) | nx.descendants(T, u)) for u in T.nodes) # map from node -> set of descendants (including itself)

        cell_states = set(labeling[u] for u in leaves)
        prune_states = cell_states - set(states)
        prune_leaves = [u for u in leaves if labeling[u] in prune_states]

        logger.info(f'Dropping {len(prune_leaves)} leaves labeled by states {prune_states}...')
        for u in prune_leaves:
            del labeling[u]
            while T.out_degree(u) == 0 and T.in_degree(u) == 1:
                v = list(T.predecessors(u))[0]
                T.remove_node(u)
                u = v

        logger.info('Building CDMP model...')

        progen_df = pd.read_csv(args.progen_matrix, index_col = 0)

        if set(progen_df.columns) != set(states):
            raise Exception("States in progen matrix do not match states in states file.")

        progenitors, progenitor_ids = fc.construct_progenitors_from_matrix(progen_df)

        start_time = time.time()
        model, cost_dict, progenitors, progenitor_ids, progenitor_variables, decision_variables = fc.build_cdmp_model(T, labeling, descendants, progenitors, progenitor_ids, enforce_tree=args.enforce_tree, mode=args.mode)
        end_time = time.time()
        build_time = end_time - start_time


        for k in args.klist:
            logger.info(f'Solving for k = {k} ...')

            model.addConstr(
                gp.quicksum(progenitor_variables[p] for p in progenitor_ids) <= k + len(states),
                name='total_progenitors'
            )

            logger.info('Solving model using Gurobi...')
            start_time = time.time()
            model.setParam('TimeLimit', args.time_limit_sec)
            model.optimize()
            end_time = time.time()
            solve_time = end_time - start_time

            logger.info('Writing output files...')

            decision_variables_t = {k_: v.X for k_, v in decision_variables.items()}
            with open(f'{args.prefix}_edge_labeling_{k}.txt', 'w') as f:
                f.write('src,dst,src_progenitor,dst_progenitor\n')
                for (u, v, p1, p2), val in decision_variables_t.items():
                    if val > 0:
                        f.write(f'{u},{v},"{progenitors[p1]}","{progenitors[p2]}"\n')

            with open(f'{args.prefix}_progenitors_{k}.txt', 'w') as f:
                for i, p in enumerate(progenitors):
                    if progenitor_variables[i].X > 0:
                        f.write(f'{set(p)}\n')

            with open(f'{args.prefix}_results_{k}.txt', 'w') as f:
                f.write('objective_score\tmodel_runtime\tmodel_buildtime\n')
                f.write(f'{model.ObjVal}\t{solve_time}\t{build_time}\n')

            c = model.getConstrByName('total_progenitors')
            model.remove(c)
            model.update()

    else:
        labeled_trees = []
        with open(args.file_locations) as file_locations:
            for line in file_locations:
                nw, metadata = line.rstrip().split("\t")
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

        for k in args.klist:
            if args.enforce_tree:
                solved_model, observed_potencies = ilp.solve_large_k_problem_tree(
                    labeled_trees, states, k - 1, weights, int(args.time_limit_sec)
                )
                out = ilp.post_process_solution_tree(solved_model, observed_potencies, states, labeled_trees)
                # solved_model = ilp.solve_large_k_problem_tree(
                #     labeled_trees, states, args.k, weights, int(args.time_limit_sec)
                # )
                # out = ilp.post_process_solution_tree(solved_model, states)

            else:
                solved_model = ilp.solve_large_k_problem(
                    labeled_trees, states, k - 1, weights, int(args.time_limit_sec)
                )
                out = ilp.post_process_solution_from_node_state_labels(solved_model)

            with open(f"{args.prefix}_results_{k}.txt", "w") as f:
                f.write("objective_score\tmodel_runtime\tk\n")
                f.write(f"{out[0]}\t{solved_model.runtime}\t{k}\n")

            with open(f"{args.prefix}_progenitors_{k}.txt", "w") as f:
                for progen in out[1]:
                    f.write(f"{progen}\n")

            with open(f"{args.prefix}_nodeLabels_{k}.txt", "w") as f:
                f.write("tree_index\tnode_name\tlabel\n")
                for node, label in out[2].items():
                    tree_index = node.split("-")[-1]
                    node_name = "-".join(node.split("-")[:-1])
                    f.write(f"{tree_index}\t{node_name}\t{label}\n")


def main_cli():
    """Entry point for command-line script"""
    arguments = get_options()
    main(arguments)

if __name__ == "__main__":
    arguments = get_options()
    main(arguments)
