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

# pytest.importorskip("scran")

def test_scran(rscript_path):
    pbmc = sc.datasets.pbmc3k_processed()
    sce.tl.scran(pbmc, rscript_path, verbosity=True)

def test_scran_copy(rscript_path):
    datapath = "./_data/pbmc10k.h5ad"
    adata = sc.read_h5ad(datapath)
    adata.layers['raw_data'] = adata.X.copy()
    adata_copy = sc.AnnData(X = np.asarray(adata.layers['raw_data'].todense()).copy(), 
                 obs = pd.DataFrame(index = adata.obs.index),
                 var = pd.DataFrame(index = adata.var.index))
    
    sce.tl.scran(adata_copy, rscript_path, verbosity=True)

