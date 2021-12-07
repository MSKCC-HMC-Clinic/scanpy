#!/usr/bin/env Rscript


# https://stackoverflow.com/questions/750786/whats-the-best-way-to-use-r-scripts-on-the-command-line-terminal
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
.libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path


install.packages(c("devtools"), repos = "http://cran.us.r-project.org")
install_github("ctlab/fgsea")
install_github("datapplab/gage")

library(devtools)
library(fgsea)
library(gage)
library(data.table)
library(ggplot2)
library(commandArgs)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least 3 arguments: if not, return an error
if (length(args)!=3) {
  stop("Missing arguments", call.=FALSE)
}

# args[1] is always .rnk or .csv preranked gene list

# if args[2] == '-file'
# input is .gmt file
# have to use gage read list...so we need to pass in the filename

# else if args[2] == '-list'
# input is just a list and can be directly read in and provided?

# use gage to read in .gmt files

# example data: idea: could be when args == 0...
# data(examplePathways)
# data(exampleRanks)

set.seed(42)
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = exampleRanks,
                  minSize  = 15,
                  maxSize  = 500,
                  eps = 0.0)

# TODO: figure out where to write intermediate file
# in a _data folder under _scripts??
fwrite(fgseaRes, file ="fgseaRes.csv", sep=",", sep2=c("", " ", ""))