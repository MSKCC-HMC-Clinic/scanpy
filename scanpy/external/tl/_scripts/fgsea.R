#!/usr/bin/env Rscript


# https://stackoverflow.com/questions/750786/whats-the-best-way-to-use-r-scripts-on-the-command-line-terminal
dir.create(Sys.getenv('R_LIBS_USER'), recursive = TRUE)  # create personal library
.libPaths(Sys.getenv('R_LIBS_USER'))  # add to the path

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

args = commandArgs(trailingOnly=TRUE)
print(args)

# test if there is 4 arguments: if not, return an error
if (length(args)!=4) {
  stop('Missing arguments', call.=FALSE)
}

# args[1]: --vanilla

# args[2] is always .rnk or .csv preranked gene list
filename = args[2]
ranks = read.csv(file=filename, header=FALSE)

hallmark_gene_type = args[3]
if (hallmark_gene_type == '--file') {
# if args[3] == '-file'
# then args[4] is .gmt file
# use gage to read in .gmt files
# have to use gage read list...so we need to pass in the filename
  print('file')
} else if (hallmark_gene_type == '--list') {
  print('list')
# else if args[3] == '-list'
# then args[4] is list
# input is just a list and can be directly read in and provided?

} else {
  stop('Missing argument: no hallmark gene set provided', call.=FALSE)
}

# example data: idea: could be when args == 0...
data(examplePathways)
data(exampleRanks)

set.seed(42)
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = exampleRanks,
                  minSize  = 15,
                  maxSize  = 500,
                  eps = 0.0)

# TODO: figure out where to write intermediate file
# in a _data folder under _scripts??
fwrite(fgseaRes, file ='fgseaRes.csv', sep=',', sep2=c('', ' ', ''))