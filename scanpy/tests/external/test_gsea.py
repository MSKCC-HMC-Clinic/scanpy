import csv
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
    Test reproducibility for type = 'gseapy'
    test input_gene_ranking_file inputs for type .csv and .rnk
    test hallmark_gene_sets_file inputs for type .gmt
    """
    # clear and create cache
    sc.external.tl.clear_cache()
    sc.external.tl.create_cache()

    # test .rnk and .csv for input_gene_ranking_file
    input_gene_ranking_rnk = './_data/gsea_data.gsea_data.rnk'
    input_gene_ranking_csv = './_data/gsea_data.csv'

    # test .gmt for hallmark_gene_sets_file
    # KEGG_2016 is primary test example for gseapy: data from https://maayanlab.cloud/Enrichr/#libraries
    hallmark_gene_sets_gmt = './_data/KEGG_2016.gmt'

    # test reproducibility for gseapy
    test_reproducibility('gseapy', input_gene_ranking_rnk, hallmark_gene_sets_gmt)
    test_reproducibility('gseapy', input_gene_ranking_csv, hallmark_gene_sets_gmt)


def test_fgsea(
    rscript_path=None
):
    """
    Test reproducibility that type = 'fgsea' module option for gsea() works
    """
    # clear and create cache
    sc.external.tl.clear_cache()
    sc.external.tl.create_cache()

    # test .rnk and .csv for input_gene_ranking_file
    input_gene_ranking_rnk = './_data/gsea_data.gsea_data.rnk'
    input_gene_ranking_csv = './_data/gsea_data.csv'

    # test .gmt for hallmark_gene_sets_file
    # KEGG_2016 is primary test example for gseapy: data from https://maayanlab.cloud/Enrichr/#libraries
    hallmark_gene_sets_gmt = './_data/KEGG_2016.gmt'

    # test reproducibility for fgsea
    test_reproducibility('fgsea', input_gene_ranking_rnk, hallmark_gene_sets_gmt, rscript_path)
    test_reproducibility('fgsea', input_gene_ranking_csv, hallmark_gene_sets_gmt, rscript_path)


def test_reproducibility(type, input_gene_ranking_file, hallmark_gene_sets_file, rscript_path=None):

    cachedir = settings.cachedir

    # define paths for first and second files for comparison
    first_test_file_path = os.path.join(cachedir,'first_df.csv')
    second_test_file_path = os.path.join(cachedir,'second_df.csv')

    sce.tl.gsea(input_gene_ranking_file=input_gene_ranking_file, hallmark_gene_sets_file=hallmark_gene_sets_file, type=type, out_dir=first_test_file_path, rscript_path=rscript_path, cache=True)
    sce.tl.gsea(input_gene_ranking_file=input_gene_ranking_file, hallmark_gene_sets_file=hallmark_gene_sets_file, type=type, out_dir=second_test_file_path, rscript_path=rscript_path, cache=True)

    first_gseapy_sample_df = pd.read_csv(first_test_file_path)
    second_gseapy_sample_df = pd.read_csv(second_test_file_path)
    assert first_gseapy_sample_df.equals(second_gseapy_sample_df)
