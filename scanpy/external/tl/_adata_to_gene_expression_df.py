from anndata.utils import asarray_sparse_dataset
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData

from typing import Optional, Union, Mapping  # Special
from typing import Sequence  # ABCs
from typing import Tuple  # Classes

def adata_to_gene_expression_df(
    adata: AnnData, 
    x: str,
    y: str,
    genes: str
):
    """\
    Generate dataframe of grouping column (leiden, louvain, etc.) with corresponding gene expression columns.
    Gene expression can come from adata.X or transformed gene expressions in adata.obsm or adata.layers.
    Grouping comes from adata.obsm. Examples include leiden, louvain clustering.
    Function creates a dataframe that can be used with joyplot function.
    
    Params
    ------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    x
        The chosen gene expression format. 
        Either adata.X, adata.obsm.[choice here], or adata.layers.[choice here]. 
    y
        Chosen grouping. 
        Must come from adata.obs. 
        Examples: "leiden", "louvain"
    genes
        List of genes whose expression user wants to include in dataframe.
        Must come from available genes in adata.var_names.
    
    Returns
    -------
    Returns dataframe with chosen grouping and expression per genes as columns.
    
    Example
    --------
    >>> adata_to_gene_expression_df(adata, "adata.X", "leiden", ["TNFRSF4","CPSF3L"])
            leiden	TNFRSF4	CPSF3L
        0	0	-0.171470	-0.280812
        1	2	-0.214582	-0.372653
        2	0	-0.376888	-0.295085
        3	4	-0.285241	-0.281735
        4	5	-0.256484	-0.220394
        
    >>> adata_to_gene_expression_df(adata, "adata.obsm['imputed_data']", "Clusters", ["TNFRSF4","CPSF3L"] )
            Clusters	TNFRSF4	CPSF3L
        0	19	-3.287377	-3.100823
        1	19	-3.283743	-3.146124
        2	19	-3.269338	-3.161905
        3	19	-3.212159	-2.647825
        4	19	-3.199698	-2.576992
    >>> adata_to_gene_expression_df(adata, "adata.layers['norm_count']", "Clusters", ["TNFRSF4","CPSF3L"] )
            Clusters	TNFRSF4	CPSF3L
        0	19	0.0	1.171362
        1	19	0.0	0.000000
        2	19	0.0	0.000000
        3	19	0.0	0.000000
        4	19	0.0	0.000000
    """
    adata.var_names_make_unique()
    
    if len(genes) > 8: 
        print("WARNING: Consider using fewer genes for better performance.")
    
    genes_of_interest = [e for e in list(adata.var_names) if e in genes]
    
    x_input = x.replace('[', '.').replace(']', '.').replace("'", "").split(".")
    if "X" in x_input: 
        df = sc.get.obs_df(
                adata,
                keys=[y, *genes_of_interest]
            )
    elif "obsm" in x_input:
        obsm_key = x_input[-2] #will be imputed_data
        var_locs = []
        for gene in genes:
            var_locs.append(adata.var_names.get_loc(str(gene)))
    
        tuple_list = []
        for var_loc in var_locs:
            tuple_list.append((obsm_key,var_loc))
        
        df = sc.get.obs_df(
                adata,
                keys=[y],
                obsm_keys=[*tuple_list] #takes obsm by index
            )
        df.columns = [df.columns[0]]+genes
        
    elif "layers" in x_input:
        layer_key = x_input[-2] #will be something like 'norm_count'
        df = pd.DataFrame(adata.obs[y]).reset_index(drop=True)
        for gene in genes:
            add = pd.DataFrame(adata.layers[layer_key][:, adata.var_names.get_loc(gene)])
            add.columns = [gene]
            df = pd.concat([df,add], axis = 1).reset_index(drop=True)
    else: 
        print("WARNING: Please provide a valid clustering input, e.g. adata.X.")
    return df.sort_values(by=[y]).reset_index(drop=True) #df of gene expression with grouping

