import pytest

import scanpy as sc
import pandas as pd
import scanpy.external as sce
from matplotlib.testing.compare import compare_images
import os
import shutil

# pytest.importorskip("gsea")


def test_gsea_barplot():
    """
    Test that gsea_barplot works
    """
    # remove old test file, if applicable
    test_file_path = '_images/test_gsea_barplot.png'
    if os.path.exists(test_file_path):
        os.remove(test_file_path)
    gseapy_df = pd.read_csv('./_data/master_gseapy_df.csv')
    sce.tl.gsea_barplot(gseapy_df, score='nes', color=['b','r'], out_dir='_images/test_gsea_barplot.png')
    assert compare_images('_images/test_gsea_barplot.png', '_images/master_gsea_barplot.png', tol=5) is None


def test_gseapy():
    """
    Test that type = 'gseapy' option for gsea() works
    """
    # remove old test file, if applicable
    test_file_path = '_data/test_gseapy_df.csv'
    if os.path.exists(test_file_path):
        os.remove(test_file_path)
    # test .rnk and .gmt input
    # KEGG_2016 is primary test example for gseapy: data from https://maayanlab.cloud/Enrichr/#libraries
    sce.tl.gsea(input_gene_ranking_file='./_data/gsea_data.gsea_data.rnk', hallmark_gene_sets_file='./_data/KEGG_2016.gmt', type='gseapy', out_dir='_data/test_gseapy_df.csv')
    master_gseapy_sample_df = pd.read_csv('_data/master_gseapy_df.csv')
    test_gseapy_sample_df = pd.read_csv('_data/test_gseapy_df.csv')
    assert test_gseapy_sample_df.equals(master_gseapy_sample_df)


def test_fgsea(
    rscript_path=None
):
    """
    Test that type = 'fgsea' module option for gsea() works
    """
    # remove old test file, if applicable
    test_file_path = '_data/test_fgsea_df.csv'
    if os.path.exists(test_file_path):
        os.remove(test_file_path)
    # test .rnk and .gmt input
    sce.tl.gsea(input_gene_ranking_file='./_data/gsea_data.gsea_data.rnk', hallmark_gene_sets_file='./_data/h.all.v7.4.symbols.gmt', type='fgsea', out_dir='_data/test_fgsea_df.csv', rscript_path=rscript_path, cache=True)
    master_fgsea_sample_df = pd.read_csv('_data/master_fgsea_df.csv')
    test_fgsea_sample_df = pd.read_csv('_data/test_fgsea_df.csv')
    assert test_fgsea_sample_df.equals(master_fgsea_sample_df)
