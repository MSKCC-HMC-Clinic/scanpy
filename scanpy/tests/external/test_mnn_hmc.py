import pytest

import numpy as np
import pandas as pd
from anndata import AnnData

import scanpy as sc
import scanpy.external as sce

# Import neighbor testing functionality
from scanpy.external.pp import mnn_hmc

# pytest.importorskip("mnn_hmc")

# input data
## Manhattan distance example 
# x _ _ x _ 
# _ _ o _ x
# o _ x _ o
sample_data_x = np.array([[0,0],[0,3],[1,4],[2,2]])
sample_data_o = np.array([[1,2], [2,0], [2,4]])
k1_neighbors_1 = 2  # 2 x neighbors for each o pt
k1_neighbors_2 = 1
k2_neighbors = 1  # 1 o neighbor for each x pt


# distances
distances_manhattan_1 = np.array([
    [0, 2, 0],
    [2, 0, 0],
    [0, 0, 1],
    [1, 0, 0]
])

# cell connectivity manhattan
connectivity_manhattan_1 = np.array([
    [0, 1, 0],
    [1, 0, 0],
    [0, 0, 1],
    [1, 0, 0]
])

# connectivity for euclidean, k=1 each direction
connectivity_euclidean = np.array([
    [0., 1., 0.],
    [0., 0., 0.],
    [0., 0., 1.],
    [1., 0., 0.]
])

def test_nn_formatted():
    # Compute neighbors, neighbors_formatted
    mutual_nearest_1 = sce.pp.mnn_hmc(arr1=sample_data_x, arr2=sample_data_o, k1=k1_neighbors_1,
            k2=k2_neighbors, dist_metric='manhattan')
    mutual_nearest_1_c = sce.pp.mnn_hmc(arr1=sample_data_x, arr2=sample_data_o, k1=k1_neighbors_1,
            k2=k2_neighbors, dist_metric='manhattan', connectivities=True)
    mutual_nearest_2 = sce.pp.mnn_hmc(arr1=sample_data_x, arr2=sample_data_o, k1=k1_neighbors_2,
            k2=k2_neighbors, dist_metric='euclidean', connectivities=True)

    assert np.allclose(mutual_nearest_1.toarray(), distances_manhattan_1)
    assert np.allclose(mutual_nearest_1_c.toarray(), connectivity_manhattan_1)
    assert np.allclose(mutual_nearest_2.toarray(), connectivity_euclidean)
