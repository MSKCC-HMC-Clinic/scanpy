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
    # remove old test file, if applicable
    test_file_path = '_images/test_gsea_barplot.png'
    if os.path.exists(test_file_path):
        os.remove(test_file_path)
    gseapy_df = pd.read_csv('./_data/master_gsea_df.csv')
    sce.tl.gsea_barplot(gseapy_df, score='nes', color=['b','r'], out_dir='_images/test_gsea_barplot.png')
    assert compare_images('_images/test_gsea_barplot.png', '_images/master_gsea_barplot.png', tol=5) is None
