import numpy as np
import pandas as pd
import scanpy as sc
import joypy
from anndata import AnnData
from scipy.stats import wasserstein_distance
import scanpy.external as sce

from typing import Optional, Union, Mapping  # Special
from typing import Sequence  # ABCs
from typing import Tuple  # Classes

def joyplot(
    df: pd.DataFrame, 
    grouping: str,
    x_label: Optional['str'] = None,
    y_label: Optional['str'] = None,
    figsize: Optional[Tuple[float, float]] = (10,8), 
    title: Optional['str'] = None,
    output_file: Optional['str'] = None,
    dpi: Optional[float] = 300,
    x_range: Optional[Tuple[float, float]] = (-3,3), 
    color: Optional['str'] = None,
    alpha: Optional[float] = 0.5,
    compute_wasserstein: Optional[bool] = True,
    view_y_axis: Optional[bool] = False
):
    """\
    Generate joyplot given a dataframe of grouping and numeric categories (e.g. gene expression) as columns. 
    Dataframe likely generated from adata_to_gene_expression_df() module.
    
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
        if output_file:                    
            fig.savefig(output_file, format = output_file.split(".")[-1], dpi = dpi)
    else:
        print("WARNING: Please make sure your dataframe contains the grouping (leiden, louvain, etc.) you are interested in as a column name.")

def df_to_joyplot(
    adata: AnnData, 
    y: str,
    genes: str,
    layer: Optional['str'] = None,
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
    view_y_axis: Optional[bool] = False
):
    """\
    Generate a dataframe of gene expression per chosen grouping and then create corresponding joyplot.
    
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
    df = sce.tl.adata_to_gene_expression_df(
        adata, y, genes, layer
    )
    joyplot(df, grouping=y, x_label, y_label, figsize, title, output_file, dpi, x_range, color, alpha, compute_wasserstein, view_y_axis)
 

