import pytest
from typing import Optional, Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import scanpy.external as sce
import os
import shutil
import pandas.api.types as ptypes
from anndata import AnnData
import scipy.io as sio
from scipy.sparse import csr_matrix

# pytest.importorskip("scran")


def test_scran_reduced(rscript_path):
    """Test functionality using subsampled adata from pbmc10k.h5ad using
        sc.pp.subsample(adata, fraction=0.10)
    """
    datapath = "./_data/0.1-pbmc10k.h5ad"
    adata = sc.read_h5ad(datapath)
    normed_adata = sce.tl.scran(adata, rscript_path)
    assert 'scran_norm' in normed_adata.layers 
    return normed_adata


def test_scran_layer(rscript_path):
    """Test functionality using adata = sc.datasets.pbmc68k_reduced with input adata.X in non-compressed format
       Test layer input
    """
    datapath = "./_data/0.1-pbmc10k.h5ad"
    adata = sc.read_h5ad(datapath)
    adata.layers['raw_data'] = adata.X.copy()
    normed_adata = sce.tl.scran(adata, rscript_path, layer="raw_data")
    assert 'scran_norm' in normed_adata.layers
    return normed_adata


def test_scran_adata_read_path(rscript_path):
    datapath = "./_data/0.1-pbmc10k.h5ad"
    normed_adata = sce.tl.scran(None, rscript_path, adata_read_path=datapath)
    assert 'scran_norm' in normed_adata.layers    
    return normed_adata
