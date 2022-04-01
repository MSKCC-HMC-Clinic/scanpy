#!/usr/bin/env Rscript

# https://stackoverflow.com/questions/750786/whats-the-best-way-to-use-r-scripts-on-the-command-line-terminal
dir.create(Sys.getenv('R_LIBS_USER'), recursive = TRUE)  # create personal library
.libPaths(Sys.getenv('R_LIBS_USER'))  # add to the path
install.packages(c('tools'), repos = 'http://cran.us.r-project.org')
install.packages(c('rtools'), repos = 'http://cran.us.r-project.org')
install.packages(c('devtools'), repos = 'http://cran.us.r-project.org')
install.packages(c('BiocManager'), repos = 'http://cran.us.r-project.org')
library(BiocManager)
install('gage')
library(devtools)

install_github('ctlab/fgsea')

library(fgsea)
library(gage)
library(data.table)
library(ggplot2)
library(tools)

args <- commandArgs(trailingOnly = TRUE)

# test if there is 4 arguments: if not, return an error
# args[1]: --vanilla
# args[2]: temp_dir_path
# args[3]: input_gene_ranking_file
# args[4]: hallmark_gene_sets_file
if (length(args) != 4) {
  stop("Missing arguments", call. = FALSE)
}

temp_dir_path <- args[2]
gene_ranking_filename <- args[3]
hallmark_gene_filename <- args[4]

# process gene ranking file
gene_ranking_file_type <- file_ext(gene_ranking_filename)

if (gene_ranking_file_type == 'rnk') {
    # .rnk files are tab deliminated
    ranked_genes0 <- read.table(
    gene_ranking_filename,
    sep = "\t")
} else if (gene_ranking_file_type == "csv") {
    ranked_genes0 <- read.csv(gene_ranking_filename)
} else {
    stop("Missing argument: no preranked gene set provided", call. = FALSE)
}

ranked_genes <- ranked_genes0[, 2]
names(ranked_genes) <- ranked_genes0[, 1]

# process hallmark gene file
hallmark_gene_file_type <- file_ext(hallmark_gene_filename)

if (hallmark_gene_file_type == "gmt") {
    gset <- gage::readList(hallmark_gene_filename)
} else {
    stop("Hallmark gene set must be of type .gmt", call. = FALSE)
}

set.seed(42)
fgseaRes <- fgsea(pathways = gset,
                  stats    = ranked_genes,
                  minSize  = 15,
                  maxSize  = 500,
                  eps = 0.0)

write_path <- paste0(temp_dir_path, "/fgseaRes.csv")
fwrite(fgseaRes, write_path, sep = ",", sep2 = c("", " ", ""))
