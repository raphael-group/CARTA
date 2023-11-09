import re
import numpy as np

from collections import defaultdict
from gurobipy import *
from itertools import combinations, product

import utils


def define_large_k_model(
    k, states, V, E, leaf_set_at_nodes, state_weights=None,
):
    k_list = list(range(k + 1))

    ### Begin defining model
    model = Model("large_k")

    ### Add variables

    # Add labeling variables
    inclusion_vars = model.addVars(
        k_list, states, vtype=GRB.BINARY, name="progen_inclusion"
    )
    label_vars = model.addVars(V, k_list, lb=0, ub=1, name="node_progen_label")
    selection_vars = model.addVars(V, states, vtype=GRB.BINARY, name="node_state_label")

    ### Add constraints

    # Enforce the root node always appears
    for s in states:
        model.addConstr(inclusion_vars[0, s] == 1, f"root_progen_inclusion")

    # Constraints to enforce that each node only has 1 progenitor labeling
    for v in V:
        model.addConstr(label_vars.sum(v, "*") == 1, f"one_progen_per_node_{v}")

    # Constraints to enforce that the size of the potency set is k, by enforcing that y_{v, p} == 0
    # iff x_{p, s} =/= z_{v, s}
    for v in V:
        for p in k_list:
            for s in states:
                model.addConstr(
                    label_vars[v, p] <= 1 - inclusion_vars[p, s] + selection_vars[v, s],
                    f"cardinality1_{v}_{s}_{p}",
                )
                model.addConstr(
                    label_vars[v, p] <= 1 + inclusion_vars[p, s] - selection_vars[v, s],
                    f"cardinality2_{v}_{s}_{p}",
                )

    # Constraints to enforce feasibility on directed edges u, v
    for u, v in E:
        for s in states:
            model.addConstr(
                selection_vars[u, s] >= selection_vars[v, s], f"feasibility_{u}_{v}_{s}"
            )

    # Constraints to enforce consistency
    for v in V:
        for s in sorted(leaf_set_at_nodes[v]):
            model.addConstr(selection_vars[v, s] == 1, f"consistency_{v}_{s}")

    # Symmetry breaking constraints to enforce rows are unique
    states_to_coeff = {}
    for i in range(len(states)):
        states_to_coeff[states[i]] = 2**i
    for p in k_list[:-1]:
        model.addConstr(
            sum([states_to_coeff[s] * inclusion_vars[p, s] for s in states])
            >= sum([states_to_coeff[s] * inclusion_vars[p + 1, s] for s in states]),
            f"symmetry_{p}_{p + 1}",
        )

    ### Define objective function, choosing the normalized version or not

    if state_weights is not None:
        objective_expression = quicksum(
            state_weights[s] * selection_vars[v, s]
            for v, s in sorted(product(V, states))
            if s not in leaf_set_at_nodes[v]
        )
    else:
        objective_expression = quicksum(
            selection_vars[v, s]
            for v, s in sorted(product(V, states))
            if s not in leaf_set_at_nodes[v]
        )

    model.setObjective(objective_expression, GRB.MINIMIZE)

    return model

def define_large_k_model_tree(
    k, states, observed_potencies, potency_counts, state_weights=None,
):
    N = len(observed_potencies)
    S = len(states)

    ### Begin defining model
    model = Model("large_k_tree")

    ### Add variables

    # Add labeling variables
    inclusion_vars = model.addVars(
        k + 1, S, vtype=GRB.BINARY, name="progen_inclusion"
    )
    analogue_vars = model.addVars(N, k + 1, vtype=GRB.BINARY, name="closest_progen_to_observed_potency")
    selection_vars = model.addVars(N, S, lb=0, ub=1, name="observed_potency_state_label")

    # Add tree constraint variables
    subset_vars = {}
    sum_to_one_vars = {}

    for c in range(k + 1):
        for d in range(k + 1):
            if c == d:
                continue
            subset_vars[c, d] = model.addVar(lb=0, ub=1, name=f"pairwise_column_inclusion_{c}_{d}")
            if c < d:
                sum_to_one_vars[c, d] = model.addVar(lb=0, ub=1, name=f"pairwise_column_disjoint_{c}_{d}")

    ### Add constraints

    # Enforce the root node always appears
    for s in range(S):
        model.addConstr(inclusion_vars[0, s] == 1, f"root_progen_inclusion")

    # Constraints to enforce that each observed potency has one closest analogue
    for n in range(N):
        model.addConstr(analogue_vars.sum(n, "*") == 1, f"one_analogue_per_potency_{n}")

    # Constraints to indicate the variables to be flipped for an observed potency
    # and its closest chosen progenitor analogue
    for n in range(N):
        for s in range(S):
            for p in range(k + 1):
                model.addConstr(
                    selection_vars[n, s] >= inclusion_vars[p, s] + analogue_vars[n, p] - 1,
                    f"flip_indicator_{p}_{n}_{s}",
                )
                model.addConstr(
                    inclusion_vars[p, s] >= selection_vars[n, s] + analogue_vars[n, p] - 1,
                    f"inclusion_constraint_{p}_{n}_{s}",
                )
    # Constraints to encode each observed potency
    for n in range(N):
        for s in range(S):
            if observed_potencies[n][s] == 1:
                model.addConstr(selection_vars[n, s] == 1, f"observed_potency_{n}_{s}")

    # Constraints to enforce tree-ness
    for c in range(k + 1):
        for d in range(k + 1):
            if c == d:
                continue
            
            for s in range(S):
                model.addConstr(inclusion_vars[c, s] <= inclusion_vars[d, s] + (1 - subset_vars[c, d]))
            
            if c < d:
                for s in range(S):
                    model.addConstr(inclusion_vars[c, s] + inclusion_vars[d, s] <= 2 - sum_to_one_vars[c, d])
                model.addConstr(1 <= subset_vars[c, d] + subset_vars[d, c] + sum_to_one_vars[c, d])

    # Symmetry breaking constraints to enforce rows are unique
    for p in range(k):
        model.addConstr(
            sum([2**s * inclusion_vars[p, s] for s in range(S)])
            >= sum([2**s * inclusion_vars[p + 1, s] for s in range(S)]),
            f"symmetry_{p}_{p + 1}",
        )

    ### Define objective function, choosing the normalized version or not

    if state_weights is not None:
        objective_expression = quicksum(
            selection_vars[n, s] * potency_counts[n] * state_weights[s]
            for n, s in sorted(product(range(N), range(S)))
            if observed_potencies[n][s] == 0
        )
    else:
        objective_expression = quicksum(
            selection_vars[n, s] * potency_counts[n]
            for n, s in sorted(product(range(N), range(S)))
            if observed_potencies[n][s] == 0
        )

    model.setObjective(objective_expression, GRB.MINIMIZE)

    return model

def post_process_solution_from_node_state_labels(model):
    vals = {}

    for var in model.getVars():
        vals[var.varName] = var.x

    node_progens = dict([(key, vals[key]) for key in vals if "node_state_label" in key])
    states_at_nodes = defaultdict(set)
    for key, value in node_progens.items():
        if value == 1:
            node, state = tuple(re.search("\[(.*?)\]", key)[0][1:-1].split(","))
            states_at_nodes[node].add(state)
    progens = []
    for i in states_at_nodes.values():
        if i not in progens:
            progens.append(i)

    return model.objVal, progens, dict(states_at_nodes)

def post_process_solution_tree(model, observed_potencies, states, trees, attribute_name = "state_labels"):
    vals = {}

    for var in model.getVars():
        vals[var.varName] = var.x

    progens = dict([(key, vals[key]) for key in vals if 'progen_inclusion' in key])

    final_progens = defaultdict(set)
    for key in progens:
        if progens[key] == 1:
            progen_index_by_state = re.search("\[(.*?)\]", key)[0][1: -1].split(",")
            final_progens[int(progen_index_by_state[0])].add(states[int(progen_index_by_state[1])])
    final_progens = list(final_progens.values())

    node_progens = dict([(key, vals[key]) for key in vals if "closest_progen_to_observed_potency" in key])

    potency_to_closest_progenitor = {}
    for key, value in node_progens.items():
        if value == 1:
            observed_potency_ind, progen_ind = tuple(re.search("\[(.*?)\]", key)[0][1:-1].split(","))
            potency_set = observed_potencies[int(observed_potency_ind)]
            tot = 0
            for i in range(len(potency_set)):
                if potency_set[i]:
                    tot += 2**i
            potency_to_closest_progenitor[tot] = final_progens[int(progen_ind)]


    states_at_nodes = {}
    for ind, tree in enumerate(trees):
        for n in tree.internal_nodes:
            leaf_set = set(
                [
                    tree.get_attribute(l, attribute_name)[0]
                    for l in tree.leaves_in_subtree(n)
                ]
            ) 
            if len(leaf_set) != 1:
                progen_state = potency_to_closest_progenitor[state_set_to_int(leaf_set, states)]
                states_at_nodes[n + "-" + str(ind)] = progen_state


    return model.objVal, final_progens, states_at_nodes

def solve_large_k_problem(
    trees,
    states,
    k,
    state_weights=None,
    time_limit_min=120,
    attribute_name="state_labels",
):

    all_V = []
    all_E = []
    all_leaf_set_at_nodes = {}

    for ind, tree in enumerate(trees):
        single_state_internal_nodes = set()
        leaf_set_at_nodes = {}
        for n in tree.internal_nodes:
            leaf_set = set(
                [
                    tree.get_attribute(l, attribute_name)[0]
                    for l in tree.leaves_in_subtree(n)
                ]
            )
            # Ensure that there are no single-state internal nodes, this will force a singleton
            # state into the progenitor set solution
            if len(leaf_set) == 1:
                single_state_internal_nodes.add(n)
            else:
                leaf_set_at_nodes[n] = leaf_set

        V = list(set(tree.nodes) - set(tree.leaves) - single_state_internal_nodes)
        E = [(e[0] + f"-{ind}", e[1] + f"-{ind}") for e in tree.edges if e[1] in V]

        all_V += [v + f"-{ind}" for v in V]
        all_E += E
        for key, value in leaf_set_at_nodes.items():
            all_leaf_set_at_nodes[key + f"-{ind}"] = value


    model = define_large_k_model(
        k,
        sorted(states),
        sorted(all_V),
        sorted(all_E),
        {i: all_leaf_set_at_nodes[i] for i in sorted(list(all_leaf_set_at_nodes.keys()))},
        state_weights,
    )

    model.setParam('TimeLimit', time_limit_min)
    model.optimize()

    return model

def state_set_to_int(state_set, states):
    tot = 0
    for i in range(len(states)):
        if states[i] in state_set:
            tot += 2**i
    return tot

def int_to_state_set(int_, states):
    bin_rep = bin(int_).replace("0b", "")
    bin_rep = bin_rep.rjust(len(states), '0')    
    return [int(i) for i in bin_rep][::-1]

def solve_large_k_problem_tree(
    trees,
    states,
    k,
    state_weights=None,
    time_limit_min=120,
    attribute_name="state_labels",
):
    
    all_observed_potencies = defaultdict(int)

    for ind, tree in enumerate(trees):
        for n in tree.internal_nodes:
            leaf_set = set(
                [
                    tree.get_attribute(l, attribute_name)[0]
                    for l in tree.leaves_in_subtree(n)
                ]
            )
            # Ensure that there are no single-state internal nodes, this will force a singleton
            # state into the progenitor set solution
            if len(leaf_set) != 1:
                all_observed_potencies[state_set_to_int(leaf_set, states)] += 1      
                
    observed_potencies = [int_to_state_set(int_, states) for int_ in all_observed_potencies]
    potency_counts = list(all_observed_potencies.values())

    model = define_large_k_model_tree(
        k,
        sorted(states),
        observed_potencies,
        potency_counts,
        state_weights,
    )

    model.setParam('TimeLimit', time_limit_min)
    model.optimize()

    return model, observed_potencies

