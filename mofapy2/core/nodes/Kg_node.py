from __future__ import division
import numpy as np

from mofapy2.core.nodes.variational_nodes import *
from mofapy2.core.gp_utils import *
import scipy as s
import time
from mofapy2.core import gpu_utils
import pandas as pd
from scipy.spatial.distance import euclidean
import copy

class Kg_Node(Node):
    """
    Sigma node to model the covariance structure K_c along the given covariate.
    This node constructs the covariate kernel.

    PARAMETERS
    ----------
    dim: dimensionality of the node (= number of covariate observations x number of latent factors)
    covariates: unique covariates for construction of the covariate part of the covariance matrix
    n_grid: number of grid points for the lengthscale parameter
    """

    def __init__(self, dim, groups, rank, sigma = 0.5, sigma_const=True):
        super().__init__(dim)
        self.groups = groups                                # covariate
        self.G = dim[1]                                             # number of observations for covariate
        self.K = dim[0]                                             # number of latent processes
        self.rank = rank
        self.sigma_const = sigma_const
        if sigma_const:
            self.sigma = np.array([sigma] * self.K)
        else:
            self.sigma = [np.array([sigma] * self.G)] * self.K

        # initialize components in node
        self.compute4init()


    def compute4init(self):
        """
        Function to initialize kernel matrix
        """
        # initialize x by full connectedness of groups (matrix of 1's)
        self.x = np.sqrt(np.ones([self.K, self.rank, self.G]) * 1/self.rank)

        # initialise kernel matrix
        self.Kmat = np.zeros([self.K, self.G, self.G])

        # initialise spectral decomposition
        self.V = np.zeros([self.K, self.G, self.G])  # eigenvectors of kernel matrix
        self.D = np.zeros([self.K, self.G])  # eigenvalues of kernel matrix

        self.compute_kernel()

    def compute_kernel(self):
        """
        Function to compute kernel matrix for all lengthscales
        """
        for k in range(self.K):
            self.compute_kernel_k(k)


    def compute_kernel_k(self, k):
        # build kernel matrix based on low-rank approximation
        if self.sigma_const:
            self.Kmat[k, :, :] = np.dot(self.x[k,:,:].transpose(), self.x[k,:,:]) + self.sigma[k] * np.eye(self.G)
        else:
            self.Kmat[k, :, :] = np.dot(self.x[k,:,:].transpose(), self.x[k,:,:]) +  np.diag(self.sigma[k])

        # compute spectral decomposition
        # Sigma = VDV^T with V^T V = I
        # important to use eigh and not eig to obtain orthogonal eigenvector (which always exist for symmetric real matrices)
        self.D[k, :], self.V[k, :, :] = s.linalg.eigh(self.Kmat[k, :, :])


    def removeFactors(self, idx, axis=1):
        self.updateDim(0, self.dim[0] - len(idx))
        self.K = self.K - 1

    def get_kernel_components_k(self, k):
        """
        Method to fetch components of kernel
        """
        return self.V[k, :, :], self.D[k, :]

    def set_parameters(self, x, sigma, k):
        self.set_x(x,k)
        self.set_sigma(sigma, k)
        self.compute_kernel_k(k)

    def set_x(self, x, k):
        self.x[k,:,:] = x

    def set_sigma(self, sigma, k):
        self.sigma[k] = sigma

    def get_x(self):
        return self.x


    def get_sigma(self):
        return self.sigma
