#!/usr/bin/env Rscript


# https://stackoverflow.com/questions/750786/whats-the-best-way-to-use-r-scripts-on-the-command-line-terminal
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
.libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path


install.packages(c("devtools"), repos = "http://cran.us.r-project.org")

library(devtools)
install_github("ctlab/fgsea")
library(fgsea)
library(data.table)
library(ggplot2)
data(examplePathways)
data(exampleRanks)
set.seed(42)
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = exampleRanks,
                  minSize  = 15,
                  maxSize  = 500)
fwrite(fgseaRes, file ="fgseaRes.csv", sep=",", sep2=c("", " ", ""))