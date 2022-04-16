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
    Add sample factor data to sc.datasets.pbmc68k_reduced for testing gene ranking

    pbmc68k_reduced includes:
        obs: 'bulk_labels', 'n_genes', 'percent_mito', 'n_counts', 'S_score', 'G2M_score', 'phase', 'louvain'
        var: 'n_counts', 'means', 'dispersions', 'dispersions_norm', 'highly_variable'
        uns: 'bulk_labels_colors', 'louvain', 'louvain_colors', 'neighbors', 'pca', 'rank_genes_groups'
        obsm: 'X_pca', 'X_umap'
        varm: 'PCs'
        obsp: 'distances', 'connectivities'

    This function adds:
        obs: 'X_factors'
    """

    adata = sc.datasets.pbmc68k_reduced()

    # add obsm: 'X_factors' (cell id by factors) with factors {0, 1, 2, 3}
    factors = pd.DataFrame(np.random.rand(len(adata.obs_names), 4), columns=['0','1','2','3'], index=adata.obs_names)
    adata.obsm['X_factors'] = factors

    return adata


def test_gene_ranking():
    """
    Test basic functionality and reproducibility
    """

    # create adata
    adata = create_sample_adata()

    # we use the 'phase' observations, which includes {'G1', 'G2M', 'S'}
    ranked_df = sce.tl.rank_gene(adata, "phase", 20)

    # test for reproducibility
    assert sce.tl.rank_gene(adata, "phase", 20).equals(ranked_df)

    # test contents and dimensions of final df:
    #       ensure row titles are as expected
    #       ensure dimensions are as expected
    assert np.array_equal(ranked_df.index.to_numpy(), np.unique(adata.obs["phase"]))
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


def test_gene_ranking_all():
    """
    Test basic functionality and reproducibility, including all normalization types
    """

    # None - non-normalized output score
    normalization_types = [None, "by_total_transition_probability", "by_number_cells", "L1", "L2", "feature_variance"]

    # create adata

    adata = create_sample_adata()

    for normalization_type in normalization_types:
        print("\n\ntype:", normalization_type)

        # we use the 'phase' observations, which includes {'G1', 'G2M', 'S'}
        ranked_df = sce.tl.rank_gene(adata=adata, cell_group_by="phase", n_components=20, normalization_type=normalization_type)
        # test contents and dimensions of final df:
        # ensure row titles are as expected
        # ensure dimensions are as expected
        assert np.array_equal(ranked_df.index.to_numpy(), np.unique(adata.obs["phase"]))

    return ranked_df
