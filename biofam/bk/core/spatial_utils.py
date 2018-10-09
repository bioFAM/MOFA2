#Copied from SpatialFA project

"""
Module with functions to initialise the model
runSingleTrial: run a single trial
runMultipleTrial: run multiple trials, optionally in Parallel (not implemented)
"""

import scipy as s
import scipy.spatial as SS
from sys import path
from time import time,sleep
import pandas as pd
import numpy as np
#from joblib import Parallel, delayed

# from spatialfa.run.init_nodes import *
# from biofam.core.BayesNet import BayesNet
# from biofam.core.utils import *

def covar_rescaling_factor(C):
    """
    Returns the rescaling factor for the Gower normalizion on covariance matrix C
    the rescaled covariance matrix has sample variance of 1
    """
    n = C.shape[0]
    P = s.eye(n) - s.ones((n,n))/float(n)
    CP = C - C.mean(0)[:, s.newaxis]
    trPCP = s.sum(P * CP)
    r = (n-1) / trPCP
    return r

def SE(X, l):
    if l == 0:
        return np.eye(X.shape[0])

    tmp = SS.distance.pdist(X,'euclidean')**2.
    tmp = SS.distance.squareform(tmp)
    cov = np.exp(-tmp/l**2.)
    # cov*=covar_rescaling_factor(cov)  # should be done outside only
    cov += 1e-3 * np.eye(X.shape[0])

    return cov

def BlockSE(X, clust, l):
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
    assert mat.shape[0] == mat.shape[1], "non-squared matrix"

    inv_mat = np.zeros([mat.shape[0], mat.shape[0]])

    for clust_ix in np.unique(clust):
        cells_ix = np.where(clust == clust_ix)[0]
        mat_tmp = mat[np.ix_(cells_ix, cells_ix)]
        import pdb; pdb.set_trace()
        inv_mat[np.ix_(cells_ix, cells_ix)] = s.linalg.inv(mat_tmp)

    return inv_mat

def buildSpatialCov(K_spatial, K_non_spatial, X):
    """
    Function to build a list of K N*N covariance function based on a SE kernel with
    multiple length scales
    """
    tmp = SS.distance.pdist(X,'euclidean')**2.
    tmp = SS.distance.squareform(tmp)
    l_grid = get_l_grid(X)
    l_all = np.tile(l_grid, K_spatial)

    k_count = 0
    cov_list = np.zeros([5*K_spatial + K_non_spatial, X.shape[0], X.shape[0]])

    count = 0
    for l in l_all:
        cov = np.exp(-tmp/l**2.)
        cov*=covar_rescaling_factor(cov)
        cov += 1e-3 * np.eye(X.shape[0])
        cov_list[count,:,:] = cov
        count+=1

    diag_cov = np.eye(X.shape[0])
    for i in range(K_non_spatial):
        cov_list[count,:,:] = diag_cov
        count+=1

    return np.array(cov_list)

def get_l_limits(X):
    Xsq = np.sum(np.square(X), 1)
    R2 = -2. * np.dot(X, X.T) + (Xsq[:, None] + Xsq[None, :])
    R2 = np.clip(R2, 0, np.inf)
    R_vals = np.unique(R2.flatten())
    R_vals = R_vals[R_vals > 1e-8]

    l_min = np.sqrt(R_vals.min()) / 2.
    l_max = np.sqrt(R_vals.max()) * 2.

    return l_min, l_max

def get_l_grid(X, n_grid = 5):
    l_min, l_max = get_l_limits(X)
    return np.logspace(np.log10(l_min), np.log10(l_max), n_grid)