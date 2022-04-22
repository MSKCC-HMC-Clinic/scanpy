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
from ... import settings

HERE = Path(__file__).parent

def scran(
    adata: AnnData = None,
    rscript_path: 'str' = None,
    adata_path: Optional['str'] = None,
    norm_adata_out_dir: Optional['str'] = 'scran_adata',
    verbosity: Optional[bool] = False,
    cache: Optional[bool] = False
) -> AnnData:

    """\
    Normalize by a varying library size using scran, where the data is first clustered based on their raw counts.
    For each cluster, it estimates the optimal normalizing factor and uses it to normalize the cells in the cluster.
    
    ----------
    adata
        Annotated data matrix
    rscript_path
        Required for fgsea, user's path to Rscript.exe
    adata_path
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

    try:
        import anndata2ri
        import poetry
    except ImportError:
        raise ImportError('Please install anndata2ri and poetry, ie `pip install anndata2ri`.')

    try:
        if adata is not None:
            assert adata_path is None
        if adata_path is not None:
            assert adata is None
    except AssertionError:
        raise AssertionError('Only one of adata and adata_path can be an input at a time')

    if adata is not None:
        # create cache directory for temporary files
        cachedir = settings.cachedir
        sce.tl.create_cache()
        adata_path = os.path.join(cachedir, "tmp_adata")
        adata.write(adata_path)
    # else:
    #     adata = sc.read_h5ad(adata_path)

    args = [adata_path, norm_adata_out_dir]
    sce.tl.execute_r_script(rscript_path, 'scran.R', args, verbosity=verbosity)

    # adata = sc.read_h5ad(norm_adata_out_dir)

    # print("in python", adata)
    # if not cache:
    #     sce.tl.clear_cache()
    
    # return adata



