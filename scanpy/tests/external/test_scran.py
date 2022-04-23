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

def test_scran_X_compressed(rscript_path):
    """Test functionality using adata = sc.datasets.pbmc68k_reduced with input adata.X in compressed format"""
    pbmc = sc.datasets.pbmc68k_reduced()
    pbmc.X = csr_matrix(pbmc.X)
    print(type(pbmc.X))
    normed_adata = sce.tl.scran(pbmc, rscript_path, verbosity=True)
    assert 'scran_norm' in normed_adata.layers
    return normed_adata


def test_scran_X_not_compressed(rscript_path):
    """Test functionality using adata = sc.datasets.pbmc68k_reduced with input adata.X in non-compressed format
       Test that scran() compresses the input automatically
    """
    pbmc = sc.datasets.pbmc68k_reduced()
    normed_adata = sce.tl.scran(pbmc, rscript_path, verbosity=True)
    assert 'scran_norm' in normed_adata.layers
    return normed_adata


def test_scran_layer(rscript_path):
    """Test functionality using adata = sc.datasets.pbmc68k_reduced with input adata.X in non-compressed format
       Test layer input
    """
    pbmc = sc.datasets.pbmc68k_reduced()
    pbmc.layers['raw_data'] = pbmc.X.copy()
    normed_adata = sce.tl.scran(pbmc, rscript_path, layer="raw_data", verbosity=True)
    assert 'scran_norm' in normed_adata.layers
    return normed_adata


def test_scran_adata_read_path(rscript_path):
    pbmc = sc.datasets.pbmc68k_reduced()
    pbmc_temp_path = "./_data/pbmc_tmp.h5ad"
    pbmc.write(pbmc_temp_path)
    normed_adata = sce.tl.scran(None, rscript_path, adata_read_path=pbmc_temp_path, verbosity=True)
    assert 'scran_norm' in normed_adata.layers    
    return normed_adata


def test_scran(rscript_path):

    # make sure this is sparse and you need to change adata.X to whichever data the user wants 
    # this or adata.layers['raw_data'] or something like that
    datapath = "./_data/pbmc10k.h5ad"
    adata = sc.read_h5ad(datapath)
    adata.layers['raw_data'] = adata.X.copy()
    normed_adata = sce.tl.scran(adata, rscript_path, verbosity=True, cache=True)

    assert 'scran_norm' in normed_adata.layers 
    return normed_adata


def test_scran_reduced(rscript_path):
    datapath = "./_data/0.1-pbmc10k.h5ad"
    adata = sc.read_h5ad(datapath)
    normed_adata = sce.tl.scran(adata, rscript_path, cache=True)
    print("in test, past scran")
    assert 'scran_norm' in normed_adata.layers 
    return normed_adata
