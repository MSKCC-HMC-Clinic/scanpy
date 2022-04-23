from typing import Optional, Union

import pandas as pd
import numpy as np
from anndata import AnnData
import scanpy as sc


def compute_entropy(
    adata: AnnData,
	sample_ids: Optional[pd.DataFrame] = None,
    adata_obs_index: Optional[str] = None,
    neighbors_index: Optional[str] = 'neighbors',
    k: Optional[int] = 30,
    n_pcs: Optional[int] = None,
    inplace: Optional[bool] = True,
    copy: bool = False,
    key_added: Optional[str] = 'entropy',
) -> Optional[Union[AnnData, pd.DataFrame]]:
    """ Computes the batch/sample entropy of each cell.
 
    For data that spans multiple batches or samples, this computes the entropy of each
    cell as quantified by Shannon Entropy (Shannon, C E 1948). Gives a measure of how 
    mixed the data is across samples. 
    
    Parameters 
    ----------
    adata
        AnnData object, which may or may not have neighbors computed already.
    sample_ids
        A dataframe indexed by cell index, must be length # cells, storing a batch or
        sample ID of each cell. Not necessary if adata_obs_index is provided.
    adata_obs_index
        Index in `adata.obs` where sample/batch IDs may be stored for each cell. For
        example, if sample labels stored in adata.obs[‘Timepoints’], input ‘Timepoints’
        here. Not necessary if user provides sample_ids.
    neighbors_index
        Index in adata.uns to check if Scanpy neighbors has been run on the data.
        Defaults to `neighbors`, the default index in Scanpy function. Optional, can 
        pass in data where neighbors has not yet been computed. 
    k
        Number of k-nearest neighbors to use in Scanpy neighbors computation if neighbors
        has not yet been computed.
    n_pcs
        Number of principal components to use in Scanpy neighbors computation if neighbors
        has not yet been computed on the adata. Default is None, matching Scanpy neighbors
        function default.
    inplace
        If not making a copy of the adata to return, this indicates whether to modify adata
        by adding computed entropy DataFrame to adata.obs. Default to True. If copy == True,
        this value must be True, else function will throw an error. If False (and copy False),
        will return pd.DataFrame with indices of cell IDs, and values of entropy.   
    copy
        Whether to return a copy of the adata with entropy stored in adata.obs. If True,
        inplace must be True as well.
    key_added 
        Defaults to `entropy`, allows user to store entropy data in 
        adata.obs[key_added + '_' + n_neighbors] if modifying adata or making copy 
        to return. 

    Returns
    -------
    :obj:‘None‘
	    By default (``copy=False``), updates ``adata`` with the following field:
	    ``adata.obs[‘key_added’ + `_` + k]`` (:class:`pd.Series`, dtype ``float64``)
		A panda series indexed by cell, with entropy computed on the k nearest neighbors 
        for the value. 

    :class:`~anndata.Anndata`
	    When ``copy=True`` is set, a copy of ``adata`` with those fields is returned.

    :obj:`~pandas.Series`
	    When inplace = False, returns the `pd.Series` alone, with index of cell ids, and datatype float64
    """
    # Preprocessing checks:
    if copy:
        if not inplace:  # must be inplace if copy is true
            raise ValueError("`copy=True` cannot be used with `inplace=False`.")
        adata = adata.copy()

    if neighbors_index not in adata.uns:  # if needed, run neighbors
        sc.pp.neighbors(adata, n_neighbors = k, n_pcs = n_pcs)
        neighbor_arr = adata.obsp['distances']  # sparse matrix of distances
    else:
        neighbor_arr = adata.obsp[neighbors_index + '_distances']

    if sample_ids is not None:  # if inputed cell -> batch id mapping
        if len(sample_ids) != neighbor_arr.shape[0]:
            raise ValueError('The number of sample/batch ids do not match the number of cells')
    elif adata_obs_index not in adata.obs.columns:
        raise ValueError('No batch ids provided and no obs key index provided')
    else:
        sample_ids = adata.obs[adata_obs_index]  # contains batch ids

    # Build n x k arr of cells x k-neighbor indices
    n_neighbors = adata.uns[neighbors_index]['params']['n_neighbors']
    num_cells = adata.shape[0]
    neighbor_inds = np.reshape(neighbor_arr.indices, newshape=(num_cells, n_neighbors-1))  # n x k-1 neighbors-doesn't include self
    cell_ints = np.arange(num_cells, dtype=int)  # add self to matrix of k nearest neighbors
    neighbor_inds = np.column_stack((cell_ints, neighbor_inds))

    # Map batch ids to ints
    batches = sample_ids.unique()
    batch_map = {}  # dict of batch, int representing id in order of batches
    for batch_ind in range(len(batches)):
        batch_map[batches[batch_ind]] = batch_ind
    cell_to_batch = np.array([batch_map[sample_id] for sample_id in sample_ids])  # n cells -> batch id int

    # Compute entropy
    output_df = pd.Series(index = adata.obs.index, dtype='float64')
    for i in range(num_cells):
        # replace cell index w batch index, get unique+counts
        _, cts = np.unique(cell_to_batch[neighbor_inds[i]], return_counts=True)
        batch_freq = cts/n_neighbors
        cell_entropy = -(sum((batch_freq)*np.log2(batch_freq)))  # Shannon's Entropy
        output_df[i] = cell_entropy

    # Save as part of anndata if inplace or copy
    if inplace or copy:
        key_k = key_added + "_" + str(n_neighbors)
        adata.obs[key_k] = output_df

    if copy:
        return adata
    elif not inplace:
        return output_df