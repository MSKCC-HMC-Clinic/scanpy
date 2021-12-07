import pytest

import scanpy as sc
import pandas as pd
import scanpy.external as sce
from matplotlib.testing.compare import compare_images

# pytest.importorskip("gsea")

def test_gsea_barplot():
    """
    Test that gsea_barplot works
    """
    gseapy_df = pd.read_csv('.\_data\gseapy_sample_df.csv')
    sce.tl.gsea_barplot(gseapy_df, score = 'nes', color = ['b','r'], out_dir = '/_images/test_gsea_barplot.png')
    assert compare_images('.\_images\test_gsea_barplot.png', '.\_images\master_gsea_barplot.png', tol=5) is None
