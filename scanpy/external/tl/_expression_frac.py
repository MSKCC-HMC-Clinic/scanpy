import anndata
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

from typing import Optional, Union, Mapping  # Special
from typing import Sequence  # ABCs
from typing import Tuple  # Classes
from anndata import AnnData

def expression_frac(
    adata: AnnData,
    adata_version: Optional["str"] = "adata.X",
    genes: Optional["list"] = ["MT-"],
    xlabel: Optional['str'] = "%MT-Content",
    output_file: Optional['str'] = None,
    dpi: Optional[float] = 300
):
    """\
    
    
    Params
    ------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
        Either adata or adata.layers['norm_log'].
    genes
        User-specified list of genes to explore (either specific gene set or prefix (["MT-"], ["RP-"])).
    xlabel 
        X axis label for figures. Default %MT-Content.
    output_file
        Name of output plot image (include ".png" extension).
    dpi
        Dots per inch specification of output plot image.
    
    Returns
    -------
    expression_frac_df
        A dataframe of rows as cell indices with their corresponding expression fraction in the columns. 

    Example
    --------
    expression_frac(adata)
    expression_frac(adata, genes = ["RP"], xlabel="%RP-Content")

  """
    sc.pp.calculate_qc_metrics(adata, inplace = True)

    if len(genes)==1:
      genes = adata.var_names[adata.var_names.str.startswith(genes[0])]      
    print(f"Gene list: {genes}.")
    index_genes = [adata.var_names.get_loc(j) for j in genes]

    if adata_version == "adata.X":
        frac = np.asarray(np.sum(adata.X[:, index_genes], axis = 1)/np.sum(adata.X, axis = 1)).squeeze() * 100
    elif "layers" in adata_version:
        layer = adata_version[adata_version.find('[')+1:adata_version.find(']')]
        frac = np.asarray(np.sum(adata.layers[layer][:, index_genes], axis = 1)/np.sum(adata.X, axis = 1)).squeeze() * 100

    expression_frac_array = np.asarray(np.sum(adata.X[:, index_genes], axis = 1)/np.sum(adata.X, axis = 1)).squeeze() * 100
    expression_frac_df = pd.DataFrame(expression_frac_array, index = adata.obs.index, columns = ["Expression fraction"])

    fig = plt.figure(figsize = (8*3, 6*1))
    ax = fig.add_subplot(1, 3, 1)
    ax.hist(expression_frac_array, 100);
    ax.set_xlabel(xlabel, fontsize = 14)
    ax.set_ylabel('Frequency', fontsize = 14)

    ax = fig.add_subplot(1, 3, 2)
    ax.scatter(adata.obs['log1p_total_counts'], expression_frac_array);
    ax.set_xlabel('Log library size', fontsize = 14)
    ax.set_ylabel(xlabel, fontsize = 14)

    ax = fig.add_subplot(1, 3, 3)
    ax.scatter(adata.obs['log1p_n_genes_by_counts'], expression_frac_array);
    ax.set_xlabel('Log num. genes per cell', fontsize = 14)
    ax.set_ylabel(xlabel, fontsize = 14)


    fig = plt.figure(figsize = (8*3, 6*1))
    plt.hist(expression_frac_array, 100);
    plt.xlabel(xlabel, fontsize = 14)
    plt.ylabel('Frequency', fontsize = 14)
    if output_file:
      plt.savefig(output_file, format = output_file.split(".")[-1], dpi = dpi)
    
    return expression_frac_df
