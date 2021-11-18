from typing import Optional

import numpy as np
import warnings
from scipy.sparse import issparse, coo_matrix, csr_matrix
from anndata import AnnData
#from .. import settings

def nearest_neighbors_formatted(
    adata: AnnData,
    key_added: Optional[str] = None,
    copy = False) -> Optional[AnnData]:
    """
    Cleans up and better displays nearest neighbor information in
    k x N arrays.

    Reformats nearest neighbor from sparse pairwise matrix to two 
    k x N arrays for k neighbors and N cells. One array contains 
    the cell ids of the k nearest neighbors in order of distance. 
    The next contains the distances between the given cell (which row of the array) 
    and the k neighbors, also in order of distance.

    This requires having run :func:‘~scanpy.pp.neighbors’

    Parameters
    ----------
    adata
        The annotated data matrix.
    key_added
        Key under which to add the arrays of k nearest neighbors.
    copy
        Copy adata or modify it inplace.

    Returns
    -------
    :obj:‘None‘
        By default (``copy=False``), updates ``adata`` with the following fields:
        ``adata.uns[‘neighbors_formatted’]`` (:class:`numpy.ndarray`, dtype ``int32``)
            Two arrays of dim (number of cells) by k-nearest-neighbors that store:
            in [`k_nearest`] the numbers of each of the k-nearest cells in 
            ascending order by distance. In [`k_nearest_dists`] the distances to
            each of the nearest k cells, in correspondence with the cell numbers 
            in [`k_nearest`]. 

    :class:`~anndata.Anndata`
        When ``copy=True`` is set, a copy of ``adata`` with those fields is 
        returned.
    """
    # Make copy of adata if param is to copy
    adata = adata.copy() if copy else adata

    # Throw error if neighbors has not yet been run
    if 'neighbors' not in adata.uns:
        raise ValueError(
        'You need to run `pp.neighbors` first to compute a neighborhood graph.')
        
    # Add custom key if key_added
    if key_added is None:
        key_added = 'neighbors_formatted'
        cell_nums_key = 'k_nearest'
        cell_dists_key = 'k_nearest_dists'
    else:
        cell_nums_key = key_added + '_k_nearest'
        cell_dists_key = key_added + '_k_nearest_dists'

    # Initialize dictionary to store outputs   
    adata.uns[key_added] = {}
    formatted_dict = adata.uns[key_added]
    
    # Get data from scanpy.pp.neighbors call 
    dists = adata.obsp['distances']     # sparse matrix of k nearest distances
    neighs_info = adata.uns['neighbors']     # get k from metadata
    num_neighbors = neighs_info['params']['n_neighbors']   

    # first list = row num, second list = column num, where non-zero
    dist_lists = dists.nonzero()    # scipy sparse matrix function
    
    # Distance output gives us k nearest cells in array of arrays
    # cell numbers, not the distances
    k_nearest_cells = [np.concatenate((np.array([dist_lists[0][i]]), dist_lists[1][i:i + num_neighbors - 1]))
                  for i in range(0, adata.n_obs*(num_neighbors-1), num_neighbors-1)]  # adata.n_obs = num cells (observations) in data
    
    # X coordinate from dist_lists[0], Y coordinate from dist_lists[1]
    # Access X,Y coordinate in dists
    k_distances = [dists[dist_lists[0][i], dist_lists[1][i]] 
                   for i in range(adata.n_obs*(num_neighbors-1))]
    k_distances_output = [np.concatenate((np.array([0]), k_distances[i: i + num_neighbors - 1])) 
                          for i in range(0, adata.n_obs*(num_neighbors-1), num_neighbors-1)] # Ensure self included by adding 0 as first distance  
    
     # Sort arrays by nearest
    k_nearest_sorted = np.array([k_nearest_cells[i][k_distances_output[i].argsort()] 
                        for i in range(adata.n_obs)])   # k nearest neighbors in order
    k_distances_output_sorted = np.array([k_distances_output[i][k_distances_output[i].argsort()] 
                                 for i in range(adata.n_obs)])     # distances to ordered k nearest neighbors
    
    formatted_dict[cell_nums_key] = k_nearest_sorted
    formatted_dict[cell_dists_key] = k_distances_output_sorted
    
    return adata if copy else None