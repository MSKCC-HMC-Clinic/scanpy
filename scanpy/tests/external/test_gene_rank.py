import pytest
import numpy as np
import scanpy as sc
import pandas as pd
import scanpy.external as sce
import os
from matplotlib.testing.compare import compare_images
from scipy.sparse import csr_matrix
import anndata as ad
import random

# pytest.importorskip("rank_gene")


def create_sample_adata():
    """
    create random, small sample adata with following properties:

    AnnData object with n_obs x n_vars = 100 x 2000
    obs: 'celltype'
    uns: 'pca', 'neighbors', 'umap'
    obsm: 'X_pca', 'X_umap', 'X_factors'
    varm: 'PCs'
    obsp: 'distances', 'connectivities'

    Where 'celltype' is {A, B, C} and factors are {0, 1, 2, 3}

    """
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

    celltypes = ['A', 'B', 'C']
    
    conditions = [random.choice(celltypes) for i in range(len(adata.obs_names))]
    
    # obs: 'celltype' (Pandas Series: cell id by string val)
    adata.obs['celltype'] = pd.Series(conditions, index=adata.obs_names)
    
    print(adata)
    print(adata.obs['celltype'])

    # obsm: 'X_factors' (cell id by factors), 
    factors = pd.DataFrame(np.random.rand(len(adata.obs_names), 4), columns=['0','1','2','3'], index=adata.obs_names)
    print("factor df?", factors, factors.shape)
    adata.obsm['X_factors'] = factors
    
    # print(adata)
    # ranked_df = sce.tl.rank_gene(adata, "celltype", 20)
    # print(ranked_df)
    return adata


def test_gene_ranking():
    """
    Test basic functionality, given for_hmc.5ad dataset
    """
  
    # remove old test file, if applicable
    test_file_path = '_data/test_gene_rank_output.csv'
    if os.path.exists(test_file_path):
        os.remove(test_file_path)
    
    adata = sc.read_h5ad('_data/test_gene_rank_input.h5ad')
    ranked_df = sce.tl.rank_gene(adata, "Celltype_myeloid", 20)

    # test for correctness (compared to google colab code)
    ranked_df.to_csv('_data/test_gene_rank_output.csv')
    test_ranked_df = pd.read_csv('./_data/test_gene_rank_output.csv')
    master_ranked_df = pd.read_csv('./_data/master_gene_rank_output.csv')
    assert test_ranked_df.equals(master_ranked_df)

    # test for reproducibility
    assert sce.tl.rank_gene(adata, "Celltype_myeloid", 20).equals(ranked_df)

    # test contents and dimensions of final df:
    #       ensure row titles are as expected
    #       ensure dimensions are as expected
    # assert ranked_df.index.to_numpy().equals(np.unique(adata.obs["Celltype_myeloid"]))
    assert np.array_equal(ranked_df.index.to_numpy(), np.unique(adata.obs["Celltype_myeloid"]))
    return ranked_df


def test_heat_map():
    """
    Test reproducibility of rank_gene_heatmap function
    """

    # remove old test file, if applicable
    test_file_path = './_images/test_heatmap.png'
    if os.path.exists(test_file_path):
        os.remove(test_file_path)

    test_df_path = './_data/test_heatmap.csv'
    ranked_df = pd.read_csv(test_df_path, index_col=0)

    sce.tl.rank_gene_heatmap(ranked_df, out_dir=test_file_path)
    assert compare_images('_images/test_heatmap.png', '_images/master_heatmap.png', tol=5) is None


def test_gene_ranking_normalization():
    """
    Test functionality of normalization types
    """

    normalization_types = ["by_total_transition_probability", "by_number_cells", "L1", "L2", "feature_variance"]

    # remove old test file, if applicable
    test_file_path = '_data/test_gene_rank_output.csv'
    if os.path.exists(test_file_path):
        os.remove(test_file_path)
    
    adata = sc.read_h5ad('_data/test_gene_rank_input.h5ad')
    
    for normalization_type in normalization_types:
        print("type:", normalization_type)
        ranked_df = sce.tl.rank_gene(adata=adata, cell_group_by="Celltype_myeloid", n_components=20, normalization_type=normalization_type)
        # test contents and dimensions of final df:
        # ensure row titles are as expected
        # ensure dimensions are as expected
        assert np.array_equal(ranked_df.index.to_numpy(), np.unique(adata.obs["Celltype_myeloid"]))

    return ranked_df
