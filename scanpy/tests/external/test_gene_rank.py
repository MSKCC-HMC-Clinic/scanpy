import pytest
import numpy as np
import scanpy as sc
import pandas as pd
import scanpy.external as sce
import os

# pytest.importorskip("rank_gene")


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
    test_df_path = './_data/test_heatmap.csv'
    ranked_df = pd.read_csv(test_df_path)
    print(ranked_df)
    return sce.tl.rank_gene_heatmap(ranked_df)

def test_gene_ranking_normalization():

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

    # # test for correctness (compared to google colab code)
    # ranked_df.to_csv('_data/test_gene_rank_output.csv')
    # test_ranked_df = pd.read_csv('./_data/test_gene_rank_output.csv')
    # master_ranked_df = pd.read_csv('./_data/master_gene_rank_output.csv')
    # assert test_ranked_df.equals(master_ranked_df)

    # # test for reproducibility
    # assert sce.tl.rank_gene(adata, "Celltype_myeloid", 20).equals(ranked_df)

    # # test contents and dimensions of final df:
    # #       ensure row titles are as expected
    # #       ensure dimensions are as expected
    # # assert ranked_df.index.to_numpy().equals(np.unique(adata.obs["Celltype_myeloid"]))
    return ranked_df
