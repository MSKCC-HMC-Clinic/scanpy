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
    y: str,
    genes: str,
    layer: Optional['str'] = None
):
    """\
    Generate dataframe of grouping column (leiden, louvain, etc.) with corresponding gene expression columns.
    Gene expression can come from adata.X or transformed gene expressions in adata.layers. None for layer indicates use of adata.X.
    Grouping comes from adata.obsm. Examples include leiden, louvain clustering.
    Function creates a dataframe that can be used with joyplot function.
    
    Params
    ------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
    layer
        Whether to use adata.X or an anndata layer. 
        Either None or the layer as it appears in adata.layers.[choice here]. 
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
    >>> adata_to_gene_expression_df(adata, "leiden", ["TNFRSF4","CPSF3L"])
            leiden	TNFRSF4	CPSF3L
        0	0	-0.171470	-0.280812
        1	2	-0.214582	-0.372653
        2	0	-0.376888	-0.295085
        3	4	-0.285241	-0.281735
        4	5	-0.256484	-0.220394
    
    >>> adata_to_gene_expression_df(adata, "Clusters", ["TNFRSF4","CPSF3L"], "norm_count")
            Clusters	TNFRSF4	CPSF3L
        0	19	0.0	1.171362
        1	19	0.0	0.000000
        2	19	0.0	0.000000
        3	19	0.0	0.000000
        4	19	0.0	0.000000
    """
#     adata.var_names_make_unique()
    # x_input = x.replace('[', '.').replace(']', '.').replace("'", "").split(".")

    if len(genes) > 8: 
        print("WARNING: Consider using fewer genes for better performance.")
    
    if layer:
        df = pd.DataFrame(adata.obs[y]).reset_index(drop=True)
        for gene in genes:
            add = pd.DataFrame(adata.layers[layer][:, adata.var_names.get_loc(gene)])
            add.columns = [gene]
            df = pd.concat([df,add], axis = 1).reset_index(drop=True)
    elif not layer:
        genes_of_interest = set([e for e in list(adata.var_names) if e in genes]) 
        df = sc.get.obs_df(
                adata,
                keys=[y, *genes_of_interest]
            )
    # elif "obsm" in x_input:
    #     obsm_key = x_input[-2] #will be imputed_data
    #     var_locs = []
    #     for gene in genes:
    #         var_locs.append(adata.var_names.get_loc(str(gene)))
    
    #     tuple_list = []
    #     for var_loc in var_locs:
    #         tuple_list.append((obsm_key,var_loc))
        
    #     df = sc.get.obs_df(
    #             adata,
    #             keys=[y],
    #             obsm_keys=[*tuple_list] #takes obsm by index
    #         )
    #     df.columns = [df.columns[0]]+genes
    else: 
        print("WARNING: Please provide a valid clustering input, e.g. adata.X.")
    return df.sort_values(by=[y]).reset_index(drop=True) #df of gene expression with grouping

