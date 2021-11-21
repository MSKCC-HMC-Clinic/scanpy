import pytest

import scanpy as sc
import scanpy.external as sce

from functools import partial
from pathlib import Path
import sys
from itertools import repeat, chain, combinations

from matplotlib.testing import setup

import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib.testing.compare import compare_images
import anndata
from anndata import AnnData
import joypy
from scipy.stats import wasserstein_distance

# pytest.importorskip("joyplot")

HERE: Path = Path(__file__).parent.parent
ROOT = HERE / '_images'
FIGS = HERE / 'figures'

def test_joyplot():
    adata = sc.datasets.pbmc68k_reduced()
    sc.tl.leiden(adata)
    #     adata.to_df()
    df = sce.tl.adata_to_gene_expression_df(
        adata, "adata.X", "leiden", ["HES4","TNFRSF4","SSU72"]
    )
    # df
    sce.tl.joyplot(df = df, grouping = "leiden", title = "Gene expression per cluster", 
                            color = ['skyblue','lightpink','red'], x_range = (-4,4), output_file = "joyplottest.png",x_label="Expression range")


    compare_images('joyplottest.png','master_joyplot.png', tol = 10)

def main():
    print("HI")
    test_joyplot()
    print("hello")

if __name__=="__main__":
    main()   
