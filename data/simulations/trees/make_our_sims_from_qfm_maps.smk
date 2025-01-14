# make_our_sims_from_qfm_maps.smk

fm_indices = [str(i).zfill(4) for i in list(range(2, 12)) + list(range(17, 27))]
tree_index_list = list(range(5))

rule all:
    input:
        expand("/n/fs/ragr-data/users/palash/carta-rebuttal/sims/trees_prob/tree/type_16/cells_100/{fm_index}_tree_{tree_index}.txt", fm_index=fm_indices, tree_index=tree_index_list),

rule run:
    input:
        graph_file="/n/fs/ragr-data/users/palash/carta_rebuttal_ryz/supplementary_data_1_fate_map_panel/fate_map{fm_index}.json",
    output:
        #prog_file="results/{ntypes}_{k}_{fatemaprep}_{ncells}_{subsample_prop}_{treerep}_{seed}",
        tree_file="/n/fs/ragr-data/users/palash/carta-rebuttal/sims/trees_prob/tree/type_16/cells_100/{fm_index}_tree_{tree_index}.txt",
    params:
        prefix="/n/fs/ragr-data/users/palash/carta-rebuttal/sims/trees_prob/tree/type_16/cells_100/{fm_index}",
        tree_ind="{tree_index}"
    log:
        std="/n/fs/ragr-data/users/palash/carta-rebuttal/sims/trees_prob/tree/type_16/cells_100/logs/{fm_index}_{tree_index}.log",
        err="/n/fs/ragr-data/users/palash/carta-rebuttal/sims/trees_prob/tree/type_16/cells_100/logs/{fm_index}_{tree_index}.err",
    shell:
        "python /n/fs/ragr-data/users/palash/carta-rebuttal/scripts/simulate_tree_from_qfm_fate_map.py --prefix {params.prefix} "
        "--tree_ind {params.tree_ind} --graph_file_location {input.graph_file} --num_sampled_per_cell_type 100 "
        "--min_sampled_per_cell_type 100 > {log.std} 2> {log.err}"
