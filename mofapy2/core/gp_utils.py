"""
Utility funcions for the GP calculation of the Sigma Node

"""

import scipy as s
import scipy.spatial as SS
import numpy as np


def covar_rescaling_factor(C):
    """
    Returns the rescaling factor for the Gower normalizion on covariance matrix C
    the rescaled covariance matrix has sample variance of 1 s.t. one obtaines an unbiased estimate of the variance explained
    (based on https://github.com/PMBio/limix/blob/master/limix/utils/preprocess.py - covar_rescaling_factor_efficient;
    https://limix.readthedocs.io/en/stable/api/limix.qc.normalise_covariance.html)
    """
    n = C.shape[0]
    P = s.eye(n) - s.ones((n,n))/float(n) # Gower’s centering matrix
    CP = C - C.mean(0)[:, s.newaxis]
    trPCP = s.sum(P * CP) # trace of doubly centered covariance matrix
    r = (n-1) / trPCP
    return r

def covar_to_corr(C):
    """
    Transforms the covariance matrix into a correlation matrix
    """
    Cdiag = np.diag(C)
    Ccor = np.diag(1/np.sqrt(Cdiag)) @ C @ np.diag(1/np.sqrt(Cdiag))
    return Ccor

def SE(X, l):
    """
    squared exponential covariance function on input X with lengthscale l
    """
    if l == 0:
        return np.eye(X.shape[0])
    noise_var = 1e-3
    tmp = SS.distance.pdist(X,'euclidean')**2.
    tmp = SS.distance.squareform(tmp)
    cov = (1-noise_var) * np.exp(-tmp/ (2 * l** 2.))
    cov += noise_var * np.eye(X.shape[0]) # avoid singularities

    return cov

def PE(X, l):
    """
    periodic covariance function on input X with lengthscale l
    """
    if l == 0:
        return np.eye(X.shape[0])
    noise_var = 1e-3
    tmp = SS.distance.pdist(X,'euclidean')
    tmp = SS.distance.squareform(tmp)
    cov = (1-noise_var) * np.cos( np.pi/l  * tmp)
    cov += noise_var * np.eye(X.shape[0]) # avoid singularities

    return cov

def BlockSE(X, clust, l):
    """
    covariance function yielding block matrix with sqaured exponential kernel per block
    """
    if l == 0:
        return np.eye(X.shape[0])

    cov = np.zeros([X.shape[0], X.shape[0]])
    for clust_ix in np.unique(clust):
        cells_ix = np.where(clust == clust_ix)[0]
        X_tmp = X[cells_ix,:]
        se_tmp = SE(X_tmp, l)
        cov[np.ix_(cells_ix, cells_ix)] = se_tmp

    return cov

def BlockInv(mat, clust):
    """
    calculate inverse of block matrix with blocks given by clust
    """
    assert mat.shape[0] == mat.shape[1], "non-squared matrix"

    inv_mat = np.zeros([mat.shape[0], mat.shape[0]])

    for clust_ix in np.unique(clust):
        cells_ix = np.where(clust == clust_ix)[0]
        mat_tmp = mat[np.ix_(cells_ix, cells_ix)]
        inv_mat[np.ix_(cells_ix, cells_ix)] = s.linalg.inv(mat_tmp) # TODO speed up with cholesky

    return inv_mat


def get_l_limits(X):
    """
    Get boundaries for the grid of lengthscales to optimize over (as implemented in spatialDE) 
    Boundaries of the grid are the shortest observed distance, divided by 2, and the longest observed distance multiplied by 2
    """
    tmp = SS.distance.pdist(X,'euclidean')**2.
    tmp = SS.distance.squareform(tmp)
    tmp_vals = np.unique(tmp.flatten())
    tmp_vals = tmp_vals[tmp_vals > 1e-8]

    l_min = np.sqrt(tmp_vals.min()) / 2.
    l_max = np.sqrt(tmp_vals.max()) * 2.

    return l_min, l_max

def get_l_grid(X, n_grid = 5):
    """
    Function to get points in a logarithmic grid for lengthscales
    """
    l_min, l_max = get_l_limits(X)
    return np.logspace(np.log10(l_min), np.log10(l_max), n_grid)
