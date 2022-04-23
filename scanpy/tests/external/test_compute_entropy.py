import pytest

import numpy as np
import pandas as pd

import scanpy as sc
import scanpy.external as sce

from scanpy.external.pp import compute_entropy

# pytest.importorskip("compute_entropy")

def test_entropy():
    '''
        Testing compute entropy method with scanpy built in dataset `moignard15`
    '''
    # Use scanpy built in moignard dataset
    test_adata = sc.datasets.moignard15()
    sce.pp.compute_entropy(test_adata, adata_obs_index = 'exp_groups', inplace = True, 
                           copy = False, key_added = 'entropy')

    # first five values for entropy 
    test_five_vals = [1.616537, 1.597417, 1.877468, 1.711080, 0.766510]
    adata_five_vals = test_adata.obs['entropy_30'][:5]
    assert np.allclose(test_five_vals, adata_five_vals)

    # test entropy per batch id
    target_group_entropy = [1.3101538182305053, 1.4328202732516984, 1.3817733417572298, 0.10271282611216405, 1.1822043071815793]
    group_entropy = []
    for group in test_adata.obs['exp_groups'].cat.categories:
        adata_group = np.where(test_adata.obs['exp_groups']==group)
        adata_temp = test_adata.copy()[adata_group]  # adata with only data from group
        group_entropy.append(np.mean(adata_temp.obs['entropy_30']))  # mean entropy of batch

    assert np.allclose(group_entropy, target_group_entropy)