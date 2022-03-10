from typing import Union, Optional, List, Tuple, Sequence

import numpy as np
import pandas as pd
import pandas as pd
from anndata import AnnData
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.colors import Colormap, Normalize

import phenograph
import scipy


def visualize_gene_signature(
    adata: AnnData,
    gene_signature: Sequence[str],
    basis: str, #Literal['pca', 'tsne', 'umap', 'diffmap', 'draw_graph_fr']
    visualization_type: str = 'zscore', # ['zscore', 'medoid']
    # clusters: int = 1, # for 'medoid' visualization_type
    cmap: Union[Colormap, str, None] = 'viridis',
    figsize: Optional[Tuple[float, float]] = (6.5,6),
    title: Optional[str] = None,
    # save: Union[bool, str, None] = None, # TODO: follow scanpy save vs show logic
    out_pth: Optional[str] = None,
    return_fig: Optional[bool] = None,
) -> Union[Figure, Axes, None]:

    """\
   
    Visualize a gene signature using the average z-score expression on the provided basis


    Parameters
    ----------
    adata
        Annotated data matrix.
    gene_signature
        List of gene names
    basis
        String that denotes a plotting tool that computed coordinates.
    visualization_type
        String that denotes the values to plot for the color mapping. Options: 'zscore' (default), 'medoid'
    title
        Figure title
    cmap
        Colormap name or sequence of strings
    return_fig
        Boolean return
    
    Returns
    -------
    RIGHT NOW: Matplotlib axes object

    TODO: If show==False a Axes or list of it    

    Example
    -------
    >>> gene_signature = ['F13A1', 'VCAN', 'CCR2', 'ITGAL', 'CX3CR1', 'CD300E', 'ACE', 'S100A8', 'S100A9', 'CSF3R', 
                     'IL1B', 'CLEC10A', 'CD74', 'FLT3', 'CSF1R']
    >>> visualize_gene_signature(adata=adata, gene_signature=gene_signature, basis='umap', return_fig=True)

    """

    if basis not in ['pca', 'tsne', 'umap']:
      raise ValueError('basis must be either: `pca`, `tsne`, or `umap`')

    if visualization_type not in ['zscore', 'medoid']:
      raise ValueError('visualization type must be either: `zscore` or `medoid`')

    basis_id = 'X_' + basis
    gene_ids = [adata.var_names.get_loc(j) for j in gene_signature]

    if visualization_type == 'zscore':
      # first get zscores of each cell by genes in gene signature
      avg_zs = np.mean(scipy.stats.zscore(adata.X[:, gene_ids]), axis = 1) # mean across rows of zscores

      # axs = ax if ax else plt.figure(figsize=figsize)
      c = avg_zs
    # elif visualization_type == 'medoid':
    #   kmedoids = KMedoids(n_clusters=clusters, random_state=0).fit(points)
      
    ax = plt.scatter(adata.obsm[basis_id][:, 0], adata.obsm[basis_id][:, 1], s = 2, c = c, 
                    cmap = cmap) 
    
    title = title if title else basis + ' colored by ' + visualization_type
    plt.title(title, fontsize = 14)
    plt.colorbar(ax)
    plt.axis('off')
    plt.savefig(out_pth, figsize=figsize)

    if return_fig:
      return ax