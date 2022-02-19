import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import seaborn as sns
import umap
import phenograph
from sklearn.neighbors import NearestNeighbors
from typing import Union, Optional, List, Tuple, Sequence
import scipy
from scipy.sparse import find, csr_matrix
from scipy.sparse.linalg import eigs
from anndata import AnnData

from matplotlib.cm import get_cmap
cols = get_cmap('tab20').colors

def rank_gene(
    adata: AnnData,
    cell_group_by: str,
    n_components: Optional[int] = 10,
    knn: Optional[int] = 30,
    metric: Optional[str] = 'euclidean',
    n_jobs: Optional[int] = -1,
    zscore: Optional[bool] = False,
    self_autocorrelate: Optional[bool] = True,
    # normalization_type: Optional[str] = None,
) -> pd.DataFrame:

    """\

    Given an anndata object, a list of genes, and a factor matrix, calculate the
    spatial autocorrelation for each factor by gene expression

    Parameters
    ----------
    adata: AnnData object
        Annotated data matrix with adata.obsm[‘X_pca’] and adata.obsm[‘X_factors’]
    cell_group_by: string
        A cell type or cluster within adata.obs
    n_components: int
        Number of diffusion components in PCA dimension reduction; default is 10
    knn: int
        Number of nearest neighbors to calculate in computing transition matrix; default is 30
    metric: string
        Chosen distance metric to compute nearest neighbor distances; default is ‘euclidean’
    n_jobs: int
        Number of cores to use in parallelizing NN calculation; default to -1 (all cores in use)
    zscore: boolean
	    Boolean indicating whether to normalize using z-score
    self_autocorrelate: boolean
	    Boolean indicating whether to include values of a cell to itself when calculating spatial autocorrelation; default true
    normalization_type: str
        TODO: test different normalization techniques on spatial autocorrelation output

    Returns
    -------    
    DataFrame of cell neighborhood types by factors

    Notes
    -------   

    Example
    -------


    """

    def run_diffusion_maps(pca_projections, n_components, knn, metric, n_jobs):
        """Run Diffusion maps using the adaptive anisotropic kernel
        :param pca_projections: PCA projections of the data
        :param n_components: Number of diffusion components
        :return: Diffusion components, corresponding eigen values and the diffusion operator
        """

        # Determine the kernel

        # Construct the distance matrix
        print('Determing nearest neighbor graph using ' + metric + ' distance metric')
        nbrs = NearestNeighbors(int(knn), metric,
                                n_jobs).fit(pca_projections.values)
        kNN = nbrs.kneighbors_graph(pca_projections.values, mode='distance')

        # Adaptive k
        adaptive_k = int(np.floor(knn / 3))
        nbrs = NearestNeighbors(int(adaptive_k),
                                metric, n_jobs).fit(pca_projections.values)
        adaptive_std = nbrs.kneighbors_graph(
            pca_projections.values, mode='distance').max(axis=1)
        adaptive_std = np.ravel(adaptive_std.todense())

        # Find the distance similarity via kernel

        # Kernel 
        print('Determing kernel')
        N = pca_projections.shape[0] # number of dim we've reduced it down to
        x, y, dists = find(kNN) # x - row indices, y - col indices, dists - non zero-values 

        # X, y specific stds
        dists = dists / adaptive_std[x] 
        W = csr_matrix((np.exp(-dists), (x, y)), shape=[N, N])
        W.setdiag(1) 
        
        # Diffusion components
        kernel = W + W.T 

        # Markov or Transition matrix
        D = np.ravel(kernel.sum(axis=1))
        D[D != 0] = 1 / D[D != 0]
        T = csr_matrix((D, (range(N), range(N))), shape=[N, N]).dot(kernel) 
        
        # Eigen value decomposition
        D, V = eigs(T, n_components, tol=1e-4, maxiter=1000) 
        D = np.real(D) # eigenvalues 
        V = np.real(V) # eigenvectors
        inds = np.argsort(D)[::-1]
        D = D[inds]
        V = V[:, inds]

        # Normalize vectors
        for i in range(V.shape[1]):
            V[:, i] = V[:, i] / np.linalg.norm(V[:, i])

        # Create results dictionary
        # T - transition matrix
        res = {'T': T, 'EigenVectors': V, 'EigenValues': D, 'kernel': kernel}
        res['EigenVectors'] = pd.DataFrame(
            res['EigenVectors'], index=pca_projections.index)
        res['EigenValues'] = pd.Series(res['EigenValues'])

        return res

    def feature_autocorrelation(feature, sim_graph, cells_of_interest, zscore, self_autocorrelate):
        # feature - cells by one factor
        # sim_graph - labeled transition matrix, see below
        # cells of interest - all cells relating to provided cell neighborhood/group (ie celltype): Question: could there be multiples of same cell id?
        # zscore
        # check for input
        
        # condition 1: sim_graph should be provided as a pandas data frame
        if isinstance(sim_graph, pd.DataFrame):
            cells_num_row = len(set(cells_of_interest).difference(set(sim_graph.index)))
            cells_num_col = len(set(cells_of_interest).difference(set(sim_graph.columns)))

            # condition 2: cells of interest must be present in sim_graph index and columns
            if (cells_num_row != 0) & (cells_num_col != 0):
                raise ValueError('Some cells of interest not in similarity graph')
        else:
            raise ValueError('sim_graph must be pandas dataframe with cell names as index and columns')
            
        # condition 3: feature must also be provided as pandas dataframe
        if isinstance(feature, pd.DataFrame):
            cells_num = len(set(cells_of_interest).difference(set(feature.index)))
            if (cells_num != 0):
                raise ValueError('Some cells of interest not in provided feature')
        else:
            raise ValueError('feature must be pandas dataframe with cell names as index')
        
        if zscore:
            feature_zs = scipy.stats.zscore(feature)
        else:
            feature_zs = feature
        
        # isolate the similarity graph for the cells of interest
        sim_graph_interest = sim_graph.loc[cells_of_interest][cells_of_interest]
        
        # isolate the feature values for the cells of interest
        feature_interest = feature_zs.loc[cells_of_interest]
        
        # compute outer product for faster computation between feature values
        feature_outer = np.outer(feature_interest.values, feature_interest.values)
        
        # product of weight with feature values
        weight_dot_feature = sim_graph_interest.values * feature_outer

        if self_autocorrelate:
            k = 0
        else:
            k = 1    
        # compute score
        # summing time: avoiding the j=i diagonal by summing the upper triangle
        score = np.sum(np.triu(weight_dot_feature, k))
        
        # TODO: play with normalization possibilities
        # normalize score in some way <- here you can play with different techniques
        #score_norm = score/len(cells_of_interest)        
        score_norm = score/np.sum(np.triu(sim_graph_interest.values, k))

        # return
        return score, score_norm
        
    dm_res = run_diffusion_maps(pd.DataFrame(adata.obsm['X_pca'], index = adata.obs_names), 
                                    n_components,
                                    knn,
                                    metric,
                                    n_jobs)
    sim_graph = pd.DataFrame(dm_res['T'].toarray(), index = adata.obs_names, columns = adata.obs_names)

    df_temp = pd.DataFrame(index = np.unique(adata.obs[cell_group_by]), 
                        columns = adata.obsm['X_factors'].columns, dtype = 'float32') # this specifically is celltype by factors

    for feature_item in adata.obsm['X_factors'].columns: # for every factor
        feature = pd.DataFrame(adata.obsm['X_factors'][feature_item], index = adata.obs.index) # this is every cell by single factor
        for item in df_temp.index: # for every celltype/condition/neighborhood notion
            cells_of_interest = adata.obs.index[adata.obs[cell_group_by] == item] # select cells that are of the right cell neighborhood
            df_temp.loc[item][feature_item] = feature_autocorrelation(feature, sim_graph, 
                                                                cells_of_interest = cells_of_interest, zscore=zscore, self_autocorrelate=self_autocorrelate)[0] # update output dataframe of the given cell neighborhood

    return df_temp

