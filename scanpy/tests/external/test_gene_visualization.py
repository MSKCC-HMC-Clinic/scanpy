import csv
import pytest

import os
import scanpy as sc
import numpy as np
import pandas as pd
import scanpy.external as sce
from matplotlib.testing.compare import compare_images
from scipy.sparse import csr_matrix
import anndata as ad
from scanpy import settings

# pytest.importorskip("gene_visualization")


def test_gene_visualization():
    """
    Test that gene visualization method works on sample data
    """
    # remove old test file, if applicable
    test_file_paths = ['_images/test.png', '_images/test_compare.png']
    for path in test_file_paths:
        if os.path.exists(path):
            os.remove(path)

    # create random, small sample adata
    # single cell gene expression generally follows negative binomial distribution
    counts = csr_matrix(np.random.poisson(1, size=(100, 2000)), dtype=np.float32)
    adata = ad.AnnData(counts)
    adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
    adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]

    # dimension reduction
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    adata.X = adata.X.toarray()  # convert compressed csr_matrix to ndarray

    genes_of_interest = ['Gene_0', 'Gene_1', 'Gene_2']
    sce.tl.visualize_gene_signature(adata=adata, gene_signature=genes_of_interest, title='UMAP colored by average z-scored expression', basis='umap', return_fig=False, out_pth='_images/test.png')

    # test reproducibility, file location
    sce.tl.visualize_gene_signature(adata=adata, gene_signature=genes_of_interest, title='UMAP colored by average z-scored expression', basis='umap', return_fig=False, out_pth='_images/test_compare.png')
    assert compare_images('_images/test.png', '_images/test_compare.png', tol=5) is None
