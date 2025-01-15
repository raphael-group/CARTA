# run_fastCARTA_many_graph_poly.smk

cell_types_list = ["poly_tree/type_8", "poly_tree/type_12", "poly_tree/type_16"]
num_cells_list = [50, 100, 200]
fm_indices = [str(i).zfill(4) for i in list(range(2, 12)) + list(range(17, 27))]
tree_index_list = list(range(5))

rule all:
    input:
        expand("/n/fs/ragr-data/users/palash/carta-rebuttal/results/carta_dag_prob/our_sims/{cell_type}/cells_{num_cell}/{fm_index}_{tree_index}_dummy.txt", num_cell=num_cells_list, cell_type=cell_types_list, fm_index=fm_indices, tree_index=tree_index_list),

rule run:
    input:
        tree_file="/n/fs/ragr-data/users/palash/carta-rebuttal/sims/trees_prob/{cell_type}/cells_{num_cell}/{fm_index}_tree_{tree_index}.txt",
    output:
        out="/n/fs/ragr-data/users/palash/carta-rebuttal/results/carta_dag_prob/our_sims/{cell_type}/cells_{num_cell}/{fm_index}_{tree_index}_dummy.txt",
    params:
        out_folder="/n/fs/ragr-data/users/palash/carta-rebuttal/results/carta_dag_prob/our_sims/{cell_type}/cells_{num_cell}/{fm_index}_{tree_index}",
        sim_folder="/n/fs/ragr-data/users/palash/carta-rebuttal/sims/trees_prob/{cell_type}/cells_{num_cell}/{fm_index}_meta_{tree_index}.txt",
        tree_ind="{tree_index}",
        fm_ind="{fm_index}",
        graph_file="/n/fs/ragr-data/users/palash/carta-rebuttal/sims/differentiation_maps_probs/{cell_type}/graph_fate_map{fm_index}.txt",
    log:
        std="/n/fs/ragr-data/users/palash/carta-rebuttal/results/carta_dag_prob/our_sims/{cell_type}/cells_{num_cell}/logs/log_{fm_index}_{tree_index}.std",
        err="/n/fs/ragr-data/users/palash/carta-rebuttal/results/carta_dag_prob/our_sims/{cell_type}/cells_{num_cell}/logs/log_{fm_index}_{tree_index}.err",
    shell:
        "python run_fastCARTA_graph.py -t {input.tree_file} --output {params.out_folder} --label {params.sim_folder} --radius 0 --graph_file {params.graph_file} --no_enforce_tree > {log.std} 2> {log.err}"
