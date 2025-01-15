# make_our_sims_from_our_maps.smk

# cell_types_list = ["graph/type_6", "graph/type_10", "graph/type_14", "tree/type_8", "tree/type_12"]
# cell_types_list = ["poly_tree/type_8", "poly_tree/type_12", "poly_tree/type_16"]
cell_types_list = ["graph/type_6", "graph/type_10"]
num_cells_list = [50, 100, 200]
fm_indices = [str(i).zfill(4) for i in list(range(2, 12)) + list(range(17, 27))]
tree_index_list = list(range(5))

rule all:
    input:
        expand("/n/fs/ragr-data/users/palash/carta-rebuttal/sims/trees_prob/{cell_type}/cells_{num_cell}/{fm_index}_tree_{tree_index}.txt", num_cell=num_cells_list, cell_type=cell_types_list, fm_index=fm_indices, tree_index=tree_index_list),

rule run:
    input:
        graph_file="/n/fs/ragr-data/users/palash/carta-rebuttal/sims/differentiation_maps_probs/{cell_type}/graph_fate_map{fm_index}.txt",
    output:
        #prog_file="results/{ntypes}_{k}_{fatemaprep}_{ncells}_{subsample_prop}_{treerep}_{seed}",
        tree_file="/n/fs/ragr-data/users/palash/carta-rebuttal/sims/trees_prob/{cell_type}/cells_{num_cell}/{fm_index}_tree_{tree_index}.txt",
    params:
        prefix="/n/fs/ragr-data/users/palash/carta-rebuttal/sims/trees_prob/{cell_type}/cells_{num_cell}/{fm_index}",
        tree_ind="{tree_index}",
        num_cell="{num_cell}"
    log:
        std="/n/fs/ragr-data/users/palash/carta-rebuttal/sims/trees_prob/{cell_type}/cells_{num_cell}/logs/{fm_index}_{tree_index}.log",
        err="/n/fs/ragr-data/users/palash/carta-rebuttal/sims/trees_prob/{cell_type}/cells_{num_cell}/logs/{fm_index}_{tree_index}.err",
    shell:
        "python /n/fs/ragr-data/users/palash/carta-rebuttal/scripts/simulate_tree_from_qfm_fate_map.py --prefix {params.prefix} "
        "--tree_ind {params.tree_ind} --graph_file_location {input.graph_file} --num_sampled_per_cell_type {params.num_cell} "
        "--min_sampled_per_cell_type {params.num_cell} --graph_format carta > {log.std} 2> {log.err}"
