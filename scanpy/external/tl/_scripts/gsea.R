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

# use gage to read in .gmt files
data(examplePathways)
data(exampleRanks)

set.seed(42)
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = exampleRanks,
                  minSize  = 15,
                  maxSize  = 500,
                  eps = 0.0)
fwrite(fgseaRes, file ="fgseaRes.csv", sep=",", sep2=c("", " ", ""))