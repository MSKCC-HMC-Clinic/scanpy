import pytest

import scanpy as sc
import pandas as pd
import scanpy.external as sce
from matplotlib.testing.compare import compare_images
import os
import shutil
from scanpy import settings

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
    Test reproducibility and whether type = 'gseapy' option for gsea()
    """
    # clear cache
    sc.external.tl.clear_cache()
    cachedir = settings.cachedir

    # define paths for first and second files for comparison
    first_test_file_path = os.path.join(cachedir,'first_gseapy_df.csv')
    second_test_file_path = os.path.join(cachedir,'second_gseapy_df.csv')

    # test .rnk and .gmt input
    # KEGG_2016 is primary test example for gseapy: data from https://maayanlab.cloud/Enrichr/#libraries
    sce.tl.gsea(input_gene_ranking_file='./_data/gsea_data.gsea_data.rnk', hallmark_gene_sets_file='./_data/KEGG_2016.gmt', type='gseapy', out_dir=first_test_file_path, cache=False)
    sce.tl.gsea(input_gene_ranking_file='./_data/gsea_data.gsea_data.rnk', hallmark_gene_sets_file='./_data/KEGG_2016.gmt', type='gseapy', out_dir=second_test_file_path, cache=False)

    first_gseapy_sample_df = pd.read_csv(first_test_file_path)
    second_gseapy_sample_df = pd.read_csv(second_test_file_path)
    assert first_gseapy_sample_df.equals(second_gseapy_sample_df)


def test_fgsea(
    rscript_path=None
):
    """
    Test reproducibility that type = 'fgsea' module option for gsea() works
    """

    sc.external.tl.clear_cache()
    cachedir = settings.cachedir

    # define paths for first and second files for comparison
    first_test_file_path = os.path.join(cachedir,'first_fgsea_df.csv')
    second_test_file_path = os.path.join(cachedir,'second_fgsea_df.csv')

    # test .rnk and .gmt input
    sce.tl.gsea(input_gene_ranking_file='./_data/gsea_data.gsea_data.rnk', hallmark_gene_sets_file='./_data/h.all.v7.4.symbols.gmt', type='fgsea', out_dir=first_test_file_path, rscript_path=rscript_path, cache=True)
    sce.tl.gsea(input_gene_ranking_file='./_data/gsea_data.gsea_data.rnk', hallmark_gene_sets_file='./_data/h.all.v7.4.symbols.gmt', type='fgsea', out_dir=second_test_file_path, rscript_path=rscript_path, cache=True)
    first_fgsea_sample_df = pd.read_csv(first_test_file_path)
    second_fgsea_sample_df = pd.read_csv(second_test_file_path)
    assert first_fgsea_sample_df.equals(second_fgsea_sample_df)
