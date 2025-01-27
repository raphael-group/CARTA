#!/usr/bin/env Rscript
#install.packages("ape", dependencies = T, repos = "http://cran.us.r-project.org")
#install.packages("qfm", dependencies = T, repos = "http://cran.us.r-project.org")
#install.packages("optparse", dependencies = T, repos = "http://cran.us.r-project.org")
#install.packages("devtools", dependencies = T, repos = "http://cran.us.r-project.org")
#install.packages("BiocManager", dependencies = T, repos = "http://cran.us.r-project.org")
#if (!require("devtools"))
#      install.packages("devtools")

#if (!require("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")

# install dependencies from Bioconductor
#BiocManager::install("ComplexHeatmap")

#devtools::install_github("Kalhor-Lab/QFM")

library(ape)
library(qfm)
library(optparse)
require(readr)
require("stringr")

option_list = list(
  make_option(c("-p", "--prefix"), type="character", 
              help="The prefix where to store the output."),
  make_option(c("-t", "--tree_file"), type="character", 
              help="The filepath to a newick file storing simulated tree."),
  make_option(c("-l", "--label_file"), type="character", 
            help="The filepath to a newick file storing labels.", default = "NULL"),
  make_option(c("-k", "--tree_time"), type="numeric", 
          help="The filepath to a newick file storing labels.", default = 15),
  make_option(c("-r", "--root_time"), type="numeric", 
        help="The filepath to a newick file storing labels.", default = 0.6)
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

tree_list = list()
tree = read.tree(opt$tree_file)
tree_list[[length(tree_list) + 1]] = tree

meta_list = list()
if (opt$label_file == "NULL") {
    states = unlist(sapply(regmatches(tree$tip.label, gregexpr("type_-[0-9]+", tree$tip.label)), function(x){
            ifelse(identical(character(0), x), NA, x)
        })
    )
    names(states) = tree$tip.label
} else {
    meta = read.table(opt$label_file, sep = "\t", header = TRUE)
    states = unlist(as.character(meta$cell_state))
    names(states) = meta$cellBC
}
meta_list[[length(meta_list) + 1]] = states

res = ice_fase_multi(tree_list, meta_list, total_time = opt$tree_time - opt$root_time, root_time = opt$root_time)
# res = ice_fase_multi(tree_list, meta_list, total_time = 1, root_time = 0.01)

write.tree(res$gr, file = paste0(opt$prefix, "_IceFaseGraph.nwk"), append = FALSE)