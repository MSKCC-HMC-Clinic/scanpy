import numpy as np
import pandas as pd
import scanpy as sc
import joypy
from anndata import AnnData
from scipy.stats import wasserstein_distance

from typing import Optional, Union, Mapping  # Special
from typing import Sequence  # ABCs
from typing import Tuple  # Classes

def adata_to_gene_expression_df(
    adata: AnnData, 
    x: str,
    y: str,
    genes: Optional[Sequence[str]] = None
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
    if genes:
        genes_of_interest = [e for e in list(adata.var_names) if e in genes]
    else: 
        genes_of_interest = list(adata.var_names)
    
    x_input = x.replace('[', '.').replace(']', '.').replace("'", "").split(".")
    if "X" in x_input: 
        df = sc.get.obs_df(
                adata,
                keys=[y, *genes_of_interest]
            )
    if "obsm" in x_input:
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
        
    if "layers" in x_input:
        layer_key = x_input[-2] #will be something like 'norm_count'
        df = pd.DataFrame(adata.obs[y]).reset_index(drop=True)
        for gene in genes:
            add = pd.DataFrame(adata.layers[layer_key][:, adata.var_names.get_loc(gene)])
            add.columns = [gene]
            df = pd.concat([df,add], axis = 1).reset_index(drop=True)
    return df.sort_values(by=[y]).reset_index(drop=True) #df of gene expression with grouping

def gene_expression_joyplot(
    df: pd.DataFrame, 
    grouping: str,
    x_label: Optional['str'] = None,
    y_label: Optional['str'] = None,
    figsize: Optional[Tuple[float, float]] = (10,8), 
    title: Optional['str'] = None,
    output_file: Optional['str'] = 'joypy_output.png',
    dpi: Optional[float] = 300,
    x_range: Optional[Tuple[float, float]] = (-3,3), 
    color: Optional['str'] = None,
    alpha: Optional[float] = 0.5,
    compute_wasserstein: Optional[bool] = True,
    view_y_axis: Optional[bool] = True
):
    """\
    Generate joyplot given a dataframe of grouping and numeric categories (e.g. gene expression) as columns. 
    Dataframe likely generated from adata_to_gene_expression_df() function.
    
    Params
    ------
    df
        Dataframe with grouping and numeric categories (e.g. gene expression) as columns. 
    grouping
        Chosen grouping. 
        Must be exact name of column as it appears in dataframe. 
        Examples: "leiden", "louvain"
    x_label
        Optional label for x axis.
    y_label
        Optional label for y axis.
    figsize
        Figure size, height by width.
    title
        Plot title.
    output_file
        Name of output plot image (include ".png" extension).
    dpi
        Dots per inch specification of output plot image.
    x_range 
        Range of x values for plot
    color
        matplotlib color for plot. See choices here: https://matplotlib.org/stable/gallery/color/named_colors.html
        Can be a list of colors (e.g. ['b','g']) with length of list corresponding to number of numeric columns compared.
    alpha
        Value between 0 and 1 that determines transparency of density plot.
    compute_wasserstein
        Boolean indicating whether or not to compute pairwise wasserstein distance between 
        provided numeric columns.
    view_y_axis
        Boolean indicating whether or not to view y axis.
        
    Returns
    -------
    Returns joyplot with gene expression distribution per chosen grouping ("leiden", "louvain", etc.).
    
    Example
    --------
    >>> gene_expression_joyplot(df = df_leiden_two_genes, grouping = "leiden", title = "Gene expression per leiden cluster", color = ['skyblue','lightpink'])
    """
    if grouping in df.columns: 
        if len(list(df.columns))>10:
            print("WARNING: Please compare fewer numeric columns for improved plot readability.")
            return

        group_list = list(set(df[grouping]))
        labels=[f"{grouping} {group}" for group in sorted(group_list)]

        all_cols = list(df.columns)
        all_cols.remove(grouping)
        numeric_columns = all_cols
        column_names= list(numeric_columns)
        
        fig, axes = joypy.joyplot(df, by=grouping, column = numeric_columns, labels=labels, overlap = 0, 
                                  linewidth=0.5, figsize=figsize, background='w', 
                                  legend = False, title = title, 
                                  x_range = x_range, color = color, alpha = alpha, linecolor=None)
        legend = axes[0].legend(loc="upper right", fontsize= 'small', facecolor="w")
        
        ax = axes[-1]
        ax.set_xlabel(x_label)
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(y_label, rotation=270)
        ax.yaxis.set_visible(True)
        ax.yaxis.set_ticks([])
        
        x_range_center = (x_range[0]+x_range[1])/(2)
        x_offset = (x_range[1]-x_range[0])/(10) 
        x_tiny_offset = (x_range[1]-x_range[0])/(200)     

        if view_y_axis:
            for i in range(len(axes)-1):
                ax = axes[i]
                y_min, y_max = ax.get_ylim()
                ax.vlines(x_range[0], 0, y_max, color="black") 
                ax.annotate(f"{int(round(y_max))}",(x_range[0]+x_tiny_offset, y_max), fontsize=10)

        if compute_wasserstein==True:
            groups_with_numeric_columns = {}
            for group in sorted(group_list):
                groups_with_numeric_columns[group] = []
                df_subset = df[df[grouping] == group] #subset by group

                for column in numeric_columns:
                    groups_with_numeric_columns[group].append(df_subset[column]) 
#             print("groups_with_numeric_columns", groups_with_numeric_columns)
            
            if len(numeric_columns) == 2:
                wds = []
                for numeric_columns in groups_with_numeric_columns.values():
                    wds.append(round(wasserstein_distance(*numeric_columns),3))
                for i in range(len(axes)-1):
                    wd = wds[i] #grabs wasserstein distance for correct group/subset/cluster
                    ax = axes[i]
                    y_min, y_max = ax.get_ylim()
                    ax.annotate(f"w.d. {wd}",(x_range_center+x_offset,y_max/4), fontsize=10)
            else: 
                print("PAIRWISE WASSERSTEIN DISTANCES BY GROUP")
                for group, numeric_columns in groups_with_numeric_columns.items():
                    print("Group:", group)
                    for i in range(len(numeric_columns)-1):
                        for j in range(i+1, len(numeric_columns)):
                            print(f"  {column_names[i]} - {column_names[j]} : {round(wasserstein_distance(numeric_columns[i], numeric_columns[j]),3)}")
                            
        fig.savefig(output_file, format = output_file.split(".")[-1], dpi = dpi)
    else:
        print("WARNING: Please make sure your dataframe contains the grouping (leiden, louvain, etc.) you are interested in as a column name.")


