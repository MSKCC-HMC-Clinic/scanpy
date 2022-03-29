from typing import Union, Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import coo_matrix, csr_matrix

from .._compat import Literal

_Metric = Literal[
    'euclidean',
    'manhattan',
    'chebyshev',
    'minkowski',
    'wminkowski',
    'seuclidean',
    'mahalanobis',
    'hamming',
    'canberra',
    'braycurtis'
]

def mutual_nearest_neighbors(
  arr1: np.ndarray,   # not clear on what input formats will be valid
  arr2: np.ndarray,
  k1: int = 30,
  k2: int = 30,
  dist_metric: _Metric = 'euclidean',
  connectivities: bool = False,
  n_jobs: int = -2
  ):
  """ Computes the mutual nearest neighbors between two datasets  

  Takes in two subsets of data or datasets across which to compute
  mutual nearest neighbors. Can compute different numbers of neighbors
  in each direction. Returns connectivities sparse matrix and distances
  sparse matrix, not guaranteed to have any k non-zero values for any
  row/column. 

  Parameters
  ----------
  arr1
    First subset of the data to compute mutual nearest neighbors of, often
    in PCA space. For n1 number of cells, and k PCA dimensions, this will
    have shape n1 x k
  arr2
    Second subset of the data for the mutual nearest neighbor computation,
    which should be in the same dimensionality space as arr1 (whether PCA
    or not). For n2 number of cells, and the same k PCA dimensions, this
    will have shape n2 x k
  k1
    For each cell in arr2, find the k1 nearest neighbors from the arr1 data.
    Default value 30.
  k2
    For each cell in arr1, find the k2 nearest neighbors from the arr2 data.
    Default value is 30.
  dist_metric
    Distance metric to compute neighbors, from sklearn. Default is ``euclidean``.
    Valid metrics are: ``euclidean``, ``manhattan``, ``chebyshev``, ``minkowski``,
    ``wminkowski``, ``seuclidean``, ``mahalanobis``, ``hamming``, ``canberra``,
    ``braycurtis``
  connectivities
    Whether to return the connectivities matrix, otherwise distance. Default `False`
  n_jobs
    Number of parallel jobs to run for sklearn neighbors function. Default value is -2.

  Returns
  -------
  :class: `scipy.sparse.csr.csr_matrix`
    If connectivities, a matrix of dimension n1 x n2 with 1 in cells (i,j) if n1_i and
    n2_j are mutual nearest neighbors, 0 elsewhere. 
    If connectivities is false, a matrix of dimension n1 x n2 where cell (i,j) is the
    distance between n1_i and n2_j if they are mutual nearest neighbors, otherwise 0.
"""
  # sklearn NearestNeighbors objects for each input array
  nbrs1 = NearestNeighbors(n_neighbors=k1, metric=dist_metric, n_jobs=n_jobs)
  nbrs2 = NearestNeighbors(n_neighbors=k2, metric=dist_metric, n_jobs=n_jobs)

  # Fit nbrs1 object to arr1, query k1 neighbors from arr2
  nbrs1.fit(arr1)
  t1_nbrs = nbrs1.kneighbors_graph(arr2, mode='distance')

  # Fit nbrs2 object to arr2, query k2 neighbors from arr1
  nbrs2.fit(arr2)
  t2_nbrs = nbrs2.kneighbors_graph(arr1, mode='distance')

  # Mututally nearest neighbors
  mnn = t2_nbrs.multiply(t1_nbrs.T)  # anywhere not mutual will have value 0
  mnn = mnn.sqrt()  # correct distances post-multiplication

  if connectivities:  # return connectivity matrix instead of distance matrix
    n_rows, n_cols = mnn.shape
    # calculate connectivities:
    row_connectivity, col_connectivity = mnn.nonzero()
    data = np.ones(len(row_connectivity))
    connectivities_matrix = coo_matrix((data, (row_connectivity, col_connectivity)), shape=(n_rows, n_cols))
    return connectivities_matrix

  return mnn
