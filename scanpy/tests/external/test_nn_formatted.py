import pytest

import numpy as np
import pandas as pd
from anndata import AnnData

import scanpy as sc
import scanpy.external as sce

# Import neighbor testing functionality
from scanpy import Neighbors
from scanpy.external.pp import nearest_neighbors_formatted

# pytest.importorskip("nn_formatted")

# the input data
X = [[1, 0], [3, 0], [5, 6], [0, 4]]
n_neighbors = 3  # includes data points themselves

# Construct Anndata obj
adata = AnnData(np.array(X))

# distances
distances_euclidean = [
    [0.0, 2.0, 4.123105525970459],
    [0.0, 2.0, 5.0],
    [0.0, 5.385164737701416, 6.324555397033691],
    [0.0, 4.123105525970459, 5.0],
]

# cell nearest
nearest_cells = [
    [0, 1, 3],
    [1, 0, 3],
    [2, 3, 1],
    [3,0,1],
]


def test_nn_formatted():
    # Compute neighbors, neighbors_formatted
    sc.pp.neighbors(adata, method='umap', n_neighbors=n_neighbors)
    sce.pp.nearest_neighbors_formatted(adata)

    k_nearest = adata.uns['neighbors_formatted']['k_nearest']
    k_nearest_dists = adata.uns['neighbors_formatted']['k_nearest_dists']

    assert np.allclose(k_nearest_dists, np.array(distances_euclidean))
    assert np.allclose(k_nearest, np.array(nearest_cells))
