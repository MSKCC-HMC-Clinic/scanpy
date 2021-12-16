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
# install_github('datapplab/gage')

library(fgsea)
library(gage)
library(data.table)
library(ggplot2)
library(tools)

args = commandArgs(trailingOnly=TRUE)

# test if there is 4 arguments: if not, return an error
if (length(args)!=4) {
  stop('Missing arguments', call.=FALSE)
}

# args[1]: --vanilla

# args[2] is always .rnk or .csv preranked gene list
filename = args[2]

file_type = file_ext(filename)

if (file_type == 'rnk') {
    # .rnk files are tab deliminated
    ranked_genes0 = read.table(
    filename,
    sep="\t")
} else if (file_type == 'csv') {
    ranked_genes0 = read.csv(filename)
} else {
    stop('Missing argument: no preranked gene set provided', call.=FALSE)
}

ranked_genes = ranked_genes0[,2]
names(ranked_genes) = ranked_genes0[,1]

hallmark_gene_type = args[3]
if (hallmark_gene_type == '--file') {
  # then args[4] is .gmt file
  # use gage to read in .gmt files
  # have to use gage read list...so we need to pass in the filename
  hallmark_gene_set = args[4]
  gset_file_type = file_ext(hallmark_gene_set)
  if (gset_file_type == 'gmt') {
      gset = gage::readList(hallmark_gene_set)
  } else {
      stop('Error: hallmark gene set must be of type .gmt or list', call.=FALSE)

  }


} else if (hallmark_gene_type == '--list') {
  # then args[4] is list
  # input is just a list and can be directly read in and provided?

  string_list = args[4]
  gset = as.list(strsplit(string_list, ",")[[1]])

} else {
  stop('Missing argument: no hallmark gene set provided', call.=FALSE)
}


# # example data: idea: could be when args == 0...
# data(examplePathways)
# data(exampleRanks)

set.seed(42)
fgseaRes <- fgsea(pathways = gset, 
                  stats    = ranked_genes,
                  minSize  = 15,
                  maxSize  = 500,
                  eps = 0.0)


# TODO: figure out where to write intermediate file
# in some cached folder, following the scanpy cachedir settings
fwrite(fgseaRes, file ="fgseaRes.csv", sep=',', sep2=c('', ' ', ''))