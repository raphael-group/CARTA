def define_and_solve_large_k_model_first_pass(k, states, V, E, leaf_set_at_nodes):
    k_list = list(range(k))
    
    ### Begin defining model
    model = Model("large_k")

    ### Add variables
    inclusion_vars = model.addVars(k_list, states, vtype=GRB.BINARY, name = "progen_inclusion")
    label_vars = model.addVars(V, k_list, vtype=GRB.BINARY, name = "node_progen_label")
    dummy_vars = model.addVars(V, states, k_list, lb = 0, ub = 1, name = "progen_state_and")
    selection_vars = model.addVars(V, states, lb = 0, ub = 1, name = "node_state_label")
    
    ### Add constraints

    # Constraints to enforce that each node only has 1 progenitor labeling
    for v in V:
        model.addConstr(label_vars.sum(v,'*') == 1, f"one_progen_per_node_{v}")

    # Constraints to express Q_{v, s, p} = X_{p, s} * Y_{v, p} for the dummy variables indicating 
    # if a progenitor has a state
    for v in V:
        for s in states:
            for p in k_list:
                model.addConstr(dummy_vars[v, s, p] <= inclusion_vars[p, s], f"mult_expression_1_{v}_{s}_{p}")
                model.addConstr(dummy_vars[v, s, p] <= label_vars[v, p], f"mult_expression_2_{v}_{s}_{p}")
                model.addConstr(dummy_vars[v, s, p] >= inclusion_vars[p, s] + label_vars[v, p] - 1, f"mult_expression_3_{v}_{s}_{p}")

    # Constraints to enforce that a node is labeled with a state if its labeled progenitor 
    # has that state
    for v in V:
        for s in states:
            model.addConstr(selection_vars[v, s] == dummy_vars.sum(v, s, '*'), f"node_state_label_{v}_{s}")

    # Constraints to enforce feasibility on directed edges u, v
    for u, v in E:
        for s in states:
            model.addConstr(selection_vars[u, s] >= selection_vars[v, s], f"feasibility_{u}_{v}_{s}")

    # Constraints to enforce consistency
    for v in V:        
        for s in leaf_set_at_nodes[v]:
            model.addConstr(selection_vars[v, s] == 1, f"consistency_{v}_{s}")
            
    ### Define objective function
    objective_expression = quicksum(selection_vars[v, s] for v, s in product(V, states) if s not in leaf_set_at_nodes[v])
    model.setObjective(objective_expression, GRB.MINIMIZE)
    
    model.optimize()
    
    return model

def post_process_solution_from_node_progen_labels(model):

    vals = {}

    for var in model.getVars():
        vals[var.varName] = var.x

    progens = dict([(key, vals[key]) for key in vals if 'progen_inclusion' in key])
    
    final_progens = defaultdict(set)
    for key in progens:
        if progens[key] == 1:
            progen_index_by_state = re.search("\[(.*?)\]", key)[0][1: -1].split(",")
            final_progens[int(progen_index_by_state[0])].add(progen_index_by_state[1])
    
    node_progens = dict([(key, vals[key]) for key in vals if 'node_progen_label' in key])
    final_node_progen_labels = {}
    for key in node_progens:
        if node_progens[key] == 1:
            node_by_progen_label_index = re.search("\[(.*?)\]", key)[0][1: -1].split(",")
            final_node_progen_labels[node_by_progen_label_index[0]] = final_progens[int(node_by_progen_label_index[1])]
    
    return model.objVal, dict(final_progens), final_node_progen_labels

def solve_large_k_one_tree(tree, states, k, attribute_name = "state_labels", normalize_cell_proportion = False, cell_proportions = None):
    
    if normalize_cell_proportion and cell_proportions is None:
        raise Exception("In order to normalize for cell proportions, please specify the cell proportions")

    single_state_internal_nodes = set()
    leaf_set_at_nodes = {}
    for n in tree.internal_nodes:
        leaf_set = set([tree.get_attribute(l, attribute_name)[0] for l in tree.leaves_in_subtree(n)])
        # Ensure that there are no single-state internal nodes, this will force a singleton
        # state into the progenitor set solution
        if len(leaf_set) == 1:
            single_state_internal_nodes.add(n)
        else:
            leaf_set_at_nodes[n] = leaf_set
        
    V = list(set(tree.nodes) - set(tree.leaves) - single_state_internal_nodes)
    E = [e for e in tree.edges if e[1] in V]
    
    model = define_and_solve_large_k_model(k, states, V, E, leaf_set_at_nodes, normalize_cell_proportion, cell_proportions)
    
    return model

def define_large_k_model_constrain_node_state_labels(
    k, states, V, E, leaf_set_at_nodes, state_weights=None, symmetry_breaking=True
):
    k_list = list(range(k + 1))

    ### Begin defining model
    model = Model("large_k")

    ### Add variables
    inclusion_vars = model.addVars(
        k_list, states, vtype=GRB.BINARY, name="progen_inclusion"
    )
    label_vars = model.addVars(V, k_list, vtype=GRB.BINARY, name="node_progen_label")
    selection_vars = model.addVars(V, states, lb=0, ub=1, name="node_state_label")

    ### Add constraints

    # Enforce the root node always appears
    for s in states:
        model.addConstr(inclusion_vars[0, s] == 1, f"root_progen_inclusion")

    # Constraints to enforce that each node only has 1 progenitor labeling
    for v in V:
        model.addConstr(label_vars.sum(v, "*") == 1, f"one_progen_per_node_{v}")

    # Constraints to enforce that a node is labeled with a state if its labeled progenitor
    # has that state
    for v in V:
        for s in states:
            for p in k_list:
                model.addConstr(
                    selection_vars[v, s] >= inclusion_vars[p, s] + label_vars[v, p] - 1,
                    f"node_state_label_{v}_{s}_{p}",
                )
                model.addConstr(
                    inclusion_vars[p, s] >= selection_vars[v, s] + label_vars[v, p] - 1,
                    f"progen_inclusion_{v}_{s}_{p}",
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

    if symmetry_breaking:
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

def define_large_k_model_constrain_node_progen_labels(
    k, states, V, E, leaf_set_at_nodes, state_weights=None, symmetry_breaking=True
):
    k_list = list(range(k + 1))

    ### Begin defining model
    model = Model("large_k")

    ### Add variables
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

    if symmetry_breaking:
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

def define_large_k_model_tree_restricted(
    k, states, V, E, leaf_set_at_nodes, state_weights=None, symmetry_breaking=True
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

    # Add tree constraint variables
    subset_vars = {}
    sum_to_one_vars = {}

    for c, d in itertools.combinations(states, 2):
        subset_vars[c, d] = model.addVar(lb=0, name="pairwise_column_inclusion_{c}_{d}")
        subset_vars[d, c] = model.addVar(lb=0, name="pairwise_column_inclusion_{d}_{c}")
        sum_to_one_vars[c, d] = model.addVar(ub=1, name="pairwise_column_disjoint_{c}_{d}")

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

    # Constraints to enforce tree-ness
    for (c, d), var in subset_vars.items():
        for k in k_list:
            model.addConstr(
                    inclusion_vars[k, c] <= inclusion_vars[k, d] + (1 - var),
                    f"treeness1_{c}_{d}_{k}",
                )

    for (c, d), var in sum_to_one_vars.items():
        for k in k_list:
            model.addConstr(
                inclusion_vars[k, c] + inclusion_vars[k, d] <= 2 - var,
                f"treeness2_{c}_{d}_{k}",
            )
        model.addConstr(
            1 <= subset_vars[c, d] + subset_vars[d, c] + var,
            f"treeness3_{c}_{d}",
        )

    if symmetry_breaking:
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

def solve_large_k_problem(
    trees,
    states,
    k,
    ilp_formulation="constrain_state_labels",
    state_weights=None,
    symmetry_breaking=True,
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

    if ilp_formulation == "constrain_state_labels":
        model = define_large_k_model_constrain_node_state_labels(
            k,
            sorted(states),
            sorted(all_V),
            sorted(all_E),
            {i: all_leaf_set_at_nodes[i] for i in sorted(list(all_leaf_set_at_nodes.keys()))},
            state_weights,
            symmetry_breaking,
        )
    elif ilp_formulation == "constrain_node_labels":
        model = define_large_k_model_constrain_node_progen_labels(
            k,
            sorted(states),
            sorted(all_V),
            sorted(all_E),
            {i: all_leaf_set_at_nodes[i] for i in sorted(list(all_leaf_set_at_nodes.keys()))},
            state_weights,
            symmetry_breaking,
        )
    elif ilp_formulation == "tree_restrict":
        model = define_large_k_model_tree_restricted(
            k,
            sorted(states),
            sorted(all_V),
            sorted(all_E),
            {i: all_leaf_set_at_nodes[i] for i in sorted(list(all_leaf_set_at_nodes.keys()))},
            state_weights,
            symmetry_breaking,
        )

    model.optimize()

    return model