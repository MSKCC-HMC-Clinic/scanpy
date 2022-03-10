import csv
import pytest

import scanpy as sc
import pandas as pd
import scanpy.external as sce
from matplotlib.testing.compare import compare_images
import os
import shutil
from scanpy import settings

# pytest.importorskip("gene_visualization")

def test_gene_visualization():
    """
    Test that gene visualization method works
    """
    adata_path = './_data/data_for_gene_signatures.h5ad'
    adata = sc.read_h5ad(adata_path)
    genes_of_interest = ['F13A1', 'VCAN', 'CCR2', 'ITGAL', 'CX3CR1', 'CD300E', 'ACE', 'S100A8', 'S100A9', 'CSF3R', 
                        'IL1B', 'CLEC10A', 'CD74', 'FLT3', 'CSF1R']
    sce.tl.visualize_gene_signature(adata=adata, gene_signature=genes_of_interest, title='UMAP colored by average z-scored expression', basis='umap', return_fig=True, out_pth='test.png')

    # # remove old test file, if applicable
    # test_file_path = '_images/test_gsea_barplot.png'
    # if os.path.exists(test_file_path):
    #     os.remove(test_file_path)
    # gseapy_df = pd.read_csv('./_data/master_gsea_df.csv')
    # sce.tl.gsea_barplot(gseapy_df, score='nes', color=['b','r'], out_dir='_images/test_gsea_barplot.png')
    # assert compare_images('_images/test_gsea_barplot.png', '_images/master_gsea_barplot.png', tol=5) is None
