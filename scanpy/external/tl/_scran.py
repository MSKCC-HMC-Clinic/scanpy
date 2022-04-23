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
import scipy
from scipy.sparse import csr_matrix, issparse
from ... import settings

HERE = Path(__file__).parent

def scran(
    adata: AnnData = None,
    rscript_path: 'str' = None,
    layer: Optional[str] = None,
    adata_read_path: Optional['str'] = None,
    adata_write_path: Optional['str'] = 'scran_adata',
    verbosity: Optional[bool] = False,
    cache: Optional[bool] = False
) -> AnnData:

    """\

    Normalize by a varying library size using scran  (Lun ATL, McCarthy DJ, Marioni JC, 2016) where the data is first clustered based on their raw counts.
    For each cluster, it estimates the optimal normalizing factor and uses it to normalize the cells in the cluster.

    Citation:
    Lun ATL, McCarthy DJ, Marioni JC (2016). “A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor.” F1000Res., 5, 2122. doi: 10.12688/f1000research.9501.2.

    ----------
    adata
        Annotated data matrix
    rscript_path
        Required for fgsea, user's path to Rscript.exe
    layer
        If provided, use adata.layers[layer] for expression values instead of adata.X
    adata_read_path
        Optional file path to input annotated data matrix
    norm_adata_out_dir
        Optional file path to output annotated data matrix
    verbosity
        Passed in to execute_r_script. Default false
    cache
        default false
    
    Returns
    -------
    AnnData with normed counts layer 'scran_norm' in layers['scran_norm']

    Example
    -------
    >>> adata_norm = scran(adata)
    >>> scran_norm = adata_norm.layers['scran_norm']
    """

    # try:
    #     import anndata2ri
    #     import poetry
    # except ImportError:
    #     raise ImportError('Please install anndata2ri and poetry, ie `pip install anndata2ri`.')

    try:
        if adata is not None:
            assert adata_read_path is None
        if adata_read_path is not None:
            assert adata is None
    except AssertionError:
        raise AssertionError('Only one of adata and adata_read_path can be an input at a time')

    # create cache directory for temporary files
    cachedir = settings.cachedir
    sce.tl.create_cache()
    
    if adata is None:
        adata = sc.read_h5ad(adata_read_path)

    X = adata.layers[layer] if layer is not None else adata.X

    if not scipy.sparse.issparse(X):
        print("not sparse! making it so")
        X = csr_matrix(X)

    mtx_input_path = os.path.join(cachedir, "input_mtx.mtx")
    mtx_output_path = os.path.join(cachedir, "output_mtx.mtx")

    sio.mmwrite(mtx_input_path, X) 

    args = [mtx_input_path, mtx_output_path]
    sce.tl.execute_r_script(rscript_path, 'scran.R', args, verbosity=verbosity)

    norm_X = sio.mmread(mtx_output_path)
    transpose_norm_X = norm_X.transpose()
    adata.layers['scran_norm'] = transpose_norm_X
    
    if adata_write_path is not None:
        adata.write(adata_write_path)

    if not cache:
        sce.tl.clear_cache()
    
    return adata




