options(repos = c(CRAN = "https://cloud.r-project.org/"))
install.packages(c("optparse", "ape"))
library(ape)
library(optparse)
rm(list=ls())
option_list <- list(
  make_option(c("-i", "--input"), 
              help="Path to phylogenetic tree in Newick format.", 
              metavar="TREE"),
  make_option(c("-o", "--output"), 
              help="Path to write the VCV matrix", 
              metavar="MATRIX")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
cat("Successfully loaded phylogenetic tree:", opt$input, "\n")
phylogenetic_tree <- read.tree(opt$input, tree.names=TRUE)
vcv_matrix <- vcv.phylo(phylogenetic_tree, cor=T)
vcv_df <- as.data.frame(vcv_matrix)
taxon_ids <- rownames(vcv_df)
vcv_df <- cbind(taxa = taxon_ids, vcv_df)
rownames(vcv_df) <- NULL
write.csv(x = vcv_df, file=opt$output, row.names = F)
cat("Successfully created VCV matrix", "\n")