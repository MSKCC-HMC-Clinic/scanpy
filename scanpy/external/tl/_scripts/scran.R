#!/usr/bin/env Rscript

# https://stackoverflow.com/questions/750786/whats-the-best-way-to-use-r-scripts-on-the-command-line-terminal
dir.create(Sys.getenv('R_LIBS_USER'), recursive = TRUE)  # create personal library
.libPaths(Sys.getenv('R_LIBS_USER'))  # add to the path

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
BiocManager::install("scuttle")
BiocManager::install("scran")

# install.packages(c('anndata'), repos = 'http://cran.us.r-project.org')

library(scuttle)
library(scran)
library(Matrix)
# library(anndata)
# library(SummarizedExperiment)

args <- commandArgs(trailingOnly = TRUE)

# test if there is 4 arguments: if not, return an error
# args[1]: --vanilla
# args[2]: temp_dir_path
# args[3]: mtx_input_path to read
# args[4]: mtx_output_path path to write
if (length(args) != 4) {
  stop("Missing arguments", call. = FALSE)
}

temp_dir_path <- args[2]
mtx_input_path <- args[3]
mtx_output_path <- args[4]

print("temp dir")
print(temp_dir_path)

print("mtx path")
print(mtx_input_path)
print("write path")
print(mtx_output_path)

print("reading mtx")
a <- Matrix::readMM(mtx_input_path)

# transpose, because single cell experiements (sce) are gene x cell instead of cell x gene
t <- t(a)

print("dimensions")
print(dim(t))

print("creating single cell experiment")
sce <- SingleCellExperiment(list(counts = t))

print("sce:\n")
print(sce)
print("test")

# set up a single cell experiment object in R, using the raw data stored in 'X' in adata anndata object
# counts(adata) <- SummarizedExperiment::assay(adata, "X") # ERROR HERE
# print(adata)

# print("Finished Setup")

print("starting clustering")
# As part of scran the cells need to be clustered
clusters <- quickCluster(sce);

print("Finished Clustering")

# Get size factors from SCRAN
sce <- computeSumFactors(sce, clusters=clusters, positive=TRUE);
print("Finished Computing Size Factors")

# Normalize using SCRAN's size factors; each cell is divided by it's cluster specific factor
# This is distinct from median library size normalization where one factor is used for all cells
sce <- logNormCounts(sce); #divide counts by scran size factors
print("Finished Normalization")

b <- logcounts(sce)
Matrix::writeMM(b, mtx_output_path)

