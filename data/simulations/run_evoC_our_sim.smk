# run_evoC_our_sim.smk

# cell_types_list = ["graph/type_6", "graph/type_10", "graph/type_14", "tree/type_8", "tree/type_12", "tree/type_16"]
# cell_types_list = ["poly_tree/type_8", "poly_tree/type_12", "poly_tree/type_16"]
cell_types_list = ["graph/type_6", "graph/type_10"]
num_cells_list = [50, 100, 200]
fm_indices = [str(i).zfill(4) for i in list(range(2, 12)) + list(range(17, 27))]
tree_index_list = list(range(5))

rule all:
    input:
        expand("/n/fs/ragr-data/users/palash/carta-rebuttal/results/evoC_prob/our_sims/{cell_type}/cells_{num_cell}/{fm_index}_{tree_index}_EvoCGraph.nwk", num_cell=num_cells_list, cell_type=cell_types_list, fm_index=fm_indices, tree_index=tree_index_list),


rule run:
    input:
        tree_file="/n/fs/ragr-data/users/palash/carta-rebuttal/sims/trees_prob/{cell_type}/cells_{num_cell}/{fm_index}_time_estimate_tree_{tree_index}.nwk",
    output:
        tree="/n/fs/ragr-data/users/palash/carta-rebuttal/results/evoC_prob/our_sims/{cell_type}/cells_{num_cell}/{fm_index}_{tree_index}_EvoCGraph.nwk",
    params:
        prefix="/n/fs/ragr-data/users/palash/carta-rebuttal/results/evoC_prob/our_sims/{cell_type}/cells_{num_cell}/{fm_index}_{tree_index}",
        meta_file="/n/fs/ragr-data/users/palash/carta-rebuttal/sims/trees_prob/{cell_type}/cells_{num_cell}/{fm_index}_meta_{tree_index}.txt",
    log:
        std="/n/fs/ragr-data/users/palash/carta-rebuttal/results/evoC_prob/our_sims/{cell_type}/cells_{num_cell}/logs/{fm_index}_{tree_index}.log",
        err="/n/fs/ragr-data/users/palash/carta-rebuttal/results/evoC_prob/our_sims/{cell_type}/cells_{num_cell}/logs/{fm_index}_{tree_index}.err",
    benchmark: "/n/fs/ragr-data/users/palash/carta-rebuttal/results/evoC_prob/our_sims/{cell_type}/cells_{num_cell}/logs/{fm_index}_{tree_index}.benchmark",
    shell:
        "python /n/fs/ragr-data/users/palash/carta-rebuttal/scripts/run_evo_coupling.py --prefix {params.prefix} --tree_file {input.tree_file}"
        " --label_file {params.meta_file} --use_branch_lengths > {log.std} 2> {log.err}"