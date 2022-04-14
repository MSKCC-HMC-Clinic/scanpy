#!/usr/bin/env Rscript

# https://stackoverflow.com/questions/750786/whats-the-best-way-to-use-r-scripts-on-the-command-line-terminal
dir.create(Sys.getenv('R_LIBS_USER'), recursive = TRUE)  # create personal library
.libPaths(Sys.getenv('R_LIBS_USER'))  # add to the path

print("inside R")

install.packages(c('devtools'), repos = 'http://cran.us.r-project.org')
install.packages(c('BiocManager'), repos = 'http://cran.us.r-project.org')
library(BiocManager)
# install('gage')
library(devtools)

# install.packages('scuttle')
BiocManager::install("scuttle")
BiocManager::install("scran")

# install.packages(c('scuttle'), repos = 'http://cran.us.r-project.org')
# install.packages(c('scran'), repos = 'http://cran.us.r-project.org')
install.packages(c('anndata'), repos = 'http://cran.us.r-project.org')

# install_github('ctlab/fgsea')

# library(fgsea)
# library(gage)
# library(data.table)
# library(ggplot2)
# library(tools)
library(scuttle)
library(scran)
# library(S4Vectors)
library(anndata)


args <- commandArgs(trailingOnly = TRUE)

# test if there is 4 arguments: if not, return an error
# args[1]: --vanilla
# args[2]: temp_dir_path
# args[3]: adata_path to read
# args[4]: norm_adata_out_dir path to write
if (length(args) != 4) {
  stop("Missing arguments", call. = FALSE)
}

temp_dir_path <- args[2]
adata_read_path <- args[3]
adata_write_path <- args[4]

print("temp dir")
print(temp_dir_path)

print("adata path")
print(adata_read_path)
print("write path")
print(adata_write_path)

adata <- read_h5ad(adata_read_path)
print(adata$layers)

# set up a single cell experiment object in R, using the raw data stored in 'X' in adata anndata object
counts(adata) <- assay(adata, "X")

print("Finished Setup")

# As part of scran the cells need to be clustered
clusters <- quickCluster(adata);
print("Finished Clustering")

# Get size factors from SCRAN
adata <- computeSumFactors(adata, clusters=clusters, positive=TRUE);
print("Finished Computing Size Factors")

# Normalize using SCRAN's size factors; each cell is divided by it's cluster specific factor
# This is distinct from median library size normalization where one factor is used for all cells
adata <- logNormCounts(adata); #divide counts by scran size factors
print("Finished Normalization")

print(adata$layers)

write_h5ad(adata, adata_write_path)
