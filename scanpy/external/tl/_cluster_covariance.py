import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import anndata
import seaborn as sns
from anndata import AnnData
import pysal
from libpysal.weights import lat2W
from esda.moran import Moran
from splot.esda import moran_scatterplot,plot_moran

from typing import Optional, Union, Mapping  # Special
from typing import Sequence  # ABCs
from typing import Tuple  # Classes

def cluster_covariance(
    adata: AnnData,
    cluster_id: float,
    adata_version: Optional["str"] = "adata.X",
    threshold: Optional[float] = 0.0,
    adata_cluster_label: Optional['str'] = 'cluster',
    title: Optional['str'] = None,
    output_file: Optional['str'] = None,
    dpi: Optional[float] = 300,
    gene_choice: Optional['str'] = 'sum', 
    morans_i: Optional[bool] = False
):
    """\
    Generates and saves a covariance matrix across genes within a cluster. 
    Reports percent of genes with covariance greater than chosen threshold (default = 0).
    Optionally reports Moran's I calculation on covariance matrix and Moran's I scatterplot.
    
    Params
    ------
    adata
        The annotated data matrix of shape `n_obs` × `n_vars`.
        Rows correspond to cells and columns to genes.
        Either adata or adata.layers['norm_log'].
    cluster_id
        Chosen cluster. 
    adata_version
        Use either adata.X or adata.layers['norm_log']. Input norm_log or default X.
    adata_cluster_label
        The label for the clusters groupings in adata, default = 'clusters'
    title
        Plot title.
    output_file
        Name of output plot image (include ".png" extension).
    dpi
        Dots per inch specification of output plot image.
    gene_choice
        How to select the genes of interest.
        Possible values: 'sum' (choose n [300] most highly expressed)
        'variance' (choose n [300] genes with highest variance)
        'highly_variable' (choose highly variable genes: 
        comes from adata.var['highly_variable'] after calling scanpy.pp.highly_variable_genes)
    morans_i
        Whether or not to perform Morans I calculation and output corresponding plots.
    
    Returns
    -------
    A measure of correlation of genes within a cluster. 

    Examples
    --------
    cluster_covariance(adata, adata_version = "adata.layers['norm_log']", cluster_id = 0, gene_choice = 'highly_variable', morans_i = True)
    cluster_covariance(adata1, cluster_id = 0, gene_choice = 'sum', morans_i = True)
    cluster_covariance(adata1, cluster_id = 0, gene_choice = 'variance', morans_i = True)

  """
    # Identify the index of cells that belong to the nth cluster
    id_cells = np.where(adata.obs[adata_cluster_label] == cluster_id)[0]

    if adata_version == "adata.X":
      data_subset = adata.X[id_cells, :]
    elif adata_version == "adata.layers['norm_log']":
      data_subset = adata.layers['norm_log'][id_cells, :]

    # Number of genes to investigate correlation structure of
    n_genes = 300

    if gene_choice == 'sum':
      # Order genes based no their expression in descending order
      gene_expr_order = np.argsort(np.sum(data_subset, axis = 0))[::-1]
    elif gene_choice == 'variance':
      gene_expr_order = np.argsort(np.var(data_subset, axis = 0)**2)[::-1]
    elif gene_choice == 'highly_variable':
      gene_expr_order = np.array([i for i,e in enumerate(adata.var['highly_variable']) if e == True])

    data_subset_ordered = data_subset[:, gene_expr_order[0:n_genes]]
    
    indices_to_delete = []
    num_cells, num_genes = data_subset_ordered.shape
    for i in range(num_genes):
      if np.var(data_subset_ordered[:,i]) == 0: 
        indices_to_delete.append(i)
    data_subset_ordered = np.delete(data_subset_ordered,indices_to_delete,1)
    print(f"WARNING: Performed covariance calculation on {n_genes-len(indices_to_delete)} of {n_genes} genes.")

    cov = np.corrcoef(data_subset_ordered.T)

    tmp = sns.clustermap(cov,cmap='bwr', yticklabels="",xticklabels="", vmax=.5, center=0, vmin=-.5, 
                        figsize = (6, 6))
    plt.title('Cluster ' + str(cluster_id))
    if output_file:
      plt.savefig(output_file, format = output_file.split(".")[-1], dpi = dpi)

    all_cov_entries = np.concatenate(cov,axis=0)
    print(f"Percentage genes greater than correlation of {threshold}: {round((len(np.where(all_cov_entries>threshold)[0])/len(all_cov_entries))*100,3)}%")

    if morans_i:
      # Create the matrix of weigthts 
      w = lat2W(cov.shape[0], cov.shape[1])

      # Crate the pysal Moran object 
      mi = Moran(cov, w)

      # Verify Moran's I results 
      print(f"Moran's I Calculation: {round(mi.I,3)}") 

      # "Spatial lag" is a weighted sum or a weighted average of the neighboring 
      # values for the variable "(normalized) gene expression" in our case
      moran_scatterplot(mi)
      plt.xlabel("Covariance values")
      plt.show()