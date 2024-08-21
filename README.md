# CARTA (cell differentiation map inference)

<!-- ![Overview of CARTA](intro.png) -->
<!--Overview of the CARTA algorithm.-->
CARTA employs a MILP to solve a constrained maximum parsimony problem to infer (i) a cell differentiatoin map and (ii) an ancestral cell type labeling for a set of cell lineage trees.

CARTA takes as input: 
* Cell lineage trees
* Terminal cell cell type annotations for each cell in each cell lineage tree
* An integer constraint specifying the number of progenitors in the inferred cell differentiation map

## Contents

  1. [Pre-requisites](#pre-requisites)
  2. [Usage instcructions](#usage)
     * [I/O formats](#io)
     * [CARTA](#carta)

<a name="pre-requisites"></a>
## Pre-requisites (see .yaml file for versions)
+ python3 (note that 3.9 is necessary in order to correctly build cassiopeia)
+ [numpy](https://numpy.org/doc/)
+ [pandas](https://pandas.pydata.org/pandas-docs/stable/index.html)
+ [gurobipy](https://www.gurobi.com/documentation/9.0/quickstart_mac/py_python_interface.html)
+ [networkx](https://networkx.org/)
+ [cassiopeia](https://github.com/YosefLab/Cassiopeia)
+ [ete3](https://pypi.org/project/ete3/)
+ (for generating simulation and real data instances) [snakemake (>=5.2.0)](https://snakemake.readthedocs.io)

<a name="usage"></a>
## Usage instructions

<a name="io"></a>
### I/O formats
The input for CARTA is 
* A tab-delimited file, which has on each line the locations of the newick and state annotation files of the set of cell lineage trees over which to infer the cell differentiation map.
    * Example: `data/gastruloid/TLS_locations.txt`
* An integer k specifying the number of progenitors in the inferred cell differentiation map. Zero-indexed; i.e. k = 0 will specify only the root progenitor.
* A file containing all terminal cell types must be provided, with each cell type on its own line.
    * example: `data/gastruloid/TLS_states.txt`

<a name="carta"></a>
### CARTA

    usage: run_ilp.py [--prefix PREFIX] [-k K] [--file_locations FILE_LOCATIONS] [--states_file STATES_FILE] [--normalize_method NORMALIZE_METHOD] [--time_limit_min TIME_LIMIT_MIN] [--enforce_tree]

    required arguments:
      --prefix PREFIX   filepath for folder at which to store output files       
      -k K  number of progenitors in the output
      --file_locations FILE_LOCATIONS   txt file with newick and state annotation file locations
      --states_file STATES_FILE  file containing the terminal states
    optional arguments:
      --normalize_method NORMALIZE_METHOD   The weights for each terminal cell state corresponding to w_s(t). Default is w_s(t) = 1 for each terminal cell type
      --time_limit_min TIME_LIMIT_MIN   The time limit in minutes. Default is 120 minutes.
      --enforce_tree    Whether to enforce that the output cell differentiation map is a tree. Default is False

An example of usage is as follows. This command can be run from the directory that contains this README file.

    python src/run_ilp.py --prefix test -k 5 --file_locations data/gastruloid/TLS_locations.txt --states_file data/gastruloid/TLS_states.txt

A cell differentiation map built from the progenitors output by CARTA and with edge weights counting the number of cells in the dataset that traverse each edge can be generated using the following command.

    python src/build_DAG_from_labeled_trees.py --prefix test --file_locations data/gastruloid/TLS_locations.txt --states_file data/gastruloid/TLS_states.txt --node_labels_file test_nodeLabels.txt

<a name="data"></a>
### Data
Currently, the newick files encoding the TLS cell lineage trees are stored in `data/gastruloid/input_trees` and the metadata files containing the cell type annotations are stored in `data/gastruloid/formatted_and_reduced_labels`.
