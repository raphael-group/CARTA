#!/usr/bin/env Rscript

library(ape)
library(qfm)
library(optparse)

option_list = list(
  make_option(c("-p", "--prefix"), type="character", 
              help="The prefix where to store the output."),
  make_option(c("-l", "--file_locations"), type="character", 
              help="The filepath to a file storing the locations of the trees and metadata in this run.")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

tree_list = list()
meta_list = list()

con = file(opt$file_locations, "r")
while ( TRUE ) {
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  files = strsplit(line, split = "\t")[[1]]
  tree_list[[length(tree_list) + 1]] = ape::read.tree(files[1])
  meta = read.table(files[2], sep = "\t", header = TRUE)
  meta_vec = unlist(as.character(meta[["cell_state"]]))
  names(meta_vec) = meta$cellBC
  meta_list[[length(meta_list) + 1]] = meta_vec
}
close(con)

res = ice_fase_multi(tree_list, meta_list, total_time = 15 - 0.6, root_time = 0.6)
write.tree(res$gr, file = paste0(opt$prefix, "_IceFaseGraph.nwk"), append = FALSE)
