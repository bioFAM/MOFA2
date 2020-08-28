from __future__ import division
import numpy as np

from mofapy2.core.nodes.variational_nodes import *
from mofapy2.core.gp_utils import *
import scipy as s
from mofapy2.core import gpu_utils


class Kc_Node(Node):
    """
    Sigma node to model the covariance structure K_c along the given covariate(s)
    This node constructs the covariate kernel matrix on a grid of lengthscale parameters and
    keeps track of the ELBO-optimal idx.

    By default a squared exponential kernel is used:
    KC_{ij} = SE(i,j) = exp(-0.5 * (i-j)^2 / (l**2))

    PARAMETERS
    ----------
    dim: dimensionality of the node (= number of latent factors x number of covariate observations)
    covariates: unique covariates for construction of the covariate part of the covariance matrix
    n_grid: number of grid points for the lengthscale parameter
    """

    def __init__(self, dim, covariates, n_grid=10):
        super().__init__(dim)
        self.covariates = covariates                                # covariate
        self.C = dim[1]                                             # number of observations for covariate
        self.K = dim[0]                                             # number of latent processes
        self.n_grid = n_grid                                        # number of grid points to optimize lengthscale on
        self.gridix = np.zeros(self.K, dtype = np.int8)             # index of the grid values for lengthscale selected per factor


        # initialize components in node
        self.compute4init()

    def compute4init(self):
        """
        Function to initiailize the grid, kernel matrix and spectral decomposition
        """
        self.l_grid = get_l_grid(self.covariates, n_grid = self.n_grid)

        # add the diagonal covariance (lengthscale 0) to the grid
        self.l_grid = np.insert(self.l_grid, 0, 0)
        self.n_grid += 1

        # initialise kernel matrix
        self.Kmat = np.zeros([self.n_grid, self.C, self.C])  # kernel matrix on lengthscale grid

        # initialise spectral decomposition
        self.V = np.zeros([self.n_grid, self.C, self.C])    # eigenvectors of kernel matrix on lengthscale grid
        self.D = np.zeros([self.n_grid, self.C])            # eigenvalues of kernel matrix on lengthscale grid

        # compute for each lengthscale the kernel matrix
        self.compute_kernel()

    def compute_kernel(self):
        """
        Function to compute kernel matrix for all lengthscales
        """
        for i in range(self.n_grid):
            self.compute_kernel_at_gridpoint(i)

    def compute_kernel_at_gridpoint(self, i):

        # build kernel matrix based on given covariance function
        self.Kmat[i, :, :] = SE(self.covariates, self.l_grid[i], zeta=0)
        # self.Kmat[i, :, :] = Cauchy(self.sample_cov_transformed, self.l_grid[i], zeta=0)

        # compute spectral decomposition
        # Kc = VDV^T with V^T V = I
        self.D[i, :], self.V[i, :, :] = s.linalg.eigh(self.Kmat[i, :,:])

    def get_ls(self):
        """
        Method to fetch ELBO-optimal length-scales
        """
        ls = np.array([self.l_grid[i] for i in self.gridix])
        return ls

    def get_kernel_components_k(self, k):
        """
        Method to ELBO optimal components of kernel for given factor k
        """
        best_ls_idx = self.gridix[k]
        return self.V[best_ls_idx, :, :], self.D[best_ls_idx,:]


    def removeFactors(self, idx, axis=1):
        self.gridix = s.delete(self.gridix, axis=0, obj=idx)
        self.updateDim(0, self.dim[0] - len(idx))
        self.K = self.K - 1


    def set_gridix(self, lidx, k):
        self.gridix[k] = lidx
        # no recomputation required as stored on grid

    def eval_at_newpoints_k(self, new_cov, k):

        Kc_new = SE(new_cov, self.l_grid[self.gridix[k]], zeta=0)

        return Kc_new
