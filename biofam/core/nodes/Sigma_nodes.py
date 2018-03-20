#Copied from SpatialFA project

# from __future__ import division
import numpy as np

from biofam.core.nodes.variational_nodes import *
from biofam.core.spatial_utils import *
import scipy as s

# TODO we need to deal with covariates too
class SigmaGrid_Node(Node):
    # dim should be the number of latent variables
    def __init__(self, dim, X, start_opt=20, n_grid=10, n_diag=0):
        super(SigmaGrid_Node,self).__init__(dim)
        self.X = X
        self.N = X.shape[0]
        self.start_opt = start_opt
        self.n_grid = n_grid
        self.iter = 0
        self.ix = [0] * dim[0]
        #self.ix = np.zeros(dim[0])  # index of the grid values to use
        self.spatial_sig = np.zeros(dim[0])
        self.precompute()
        self.compute_cov()

        self.n_diag = n_diag  # only used for sampling

    def precompute(self):
        self.l_grid = get_l_grid(self.X, n_grid=self.n_grid)

        self.l_grid = self.l_grid[1:-1]
        self.n_grid-=2

        self.l_grid = np.insert(self.l_grid, 0, 0)
        self.n_grid += 1

        self.grid_cov = np.zeros([self.n_grid, self.N, self.N])
        self.grid_cov_inv = np.zeros([self.n_grid, self.N, self.N])
        self.grid_cov_inv_diag = np.zeros([self.n_grid, self.N])

    def compute_cov(self):
        for i in range(self.n_grid):
            # intialisation of covariance
            self.grid_cov[i,:,:] = SE(self.X, self.l_grid[i])
            self.grid_cov[i,:,:] *= covar_rescaling_factor(self.grid_cov[i,:,:])

            # intialisation of inverse covariance
            self.grid_cov_inv[i,:,:] = s.linalg.inv(self.grid_cov[i,:,:])  # TODO speed up

            # initialisation of diagonal terms
            self.grid_cov_inv_diag[i,:] = s.diag(self.grid_cov_inv[i,:,:])

    def getExpectations(self):
        cov = np.array([self.grid_cov[i,:,:] for i in self.ix])
        inv = np.array([self.grid_cov_inv[i,:,:] for i in self.ix])
        inv_diag = np.array([self.grid_cov_inv_diag[i,:] for i in self.ix])
        return {'cov':cov, 'inv': inv, 'inv_diag':inv_diag, 'E':cov}

    def get_ls(self):
        ls = np.array([self.l_grid[i] for i in self.ix])
        return ls

    def getParameters(self):
        ls = self.get_ls()
        return {'l':ls, 'sig': self.spatial_sig, 'X':self.X}

    def removeFactors(self, idx, axis=1):
        self.ix = s.delete(self.ix, axis=0, obj=idx)
        self.spatial_sig = s.delete(self.spatial_sig, axis=0, obj=idx)
        self.updateDim(0, self.dim[0] - 1)

    def optimise(self):
        # as the multiple latent variables are independent in the elbo term,
        # no need to test for all possible length scale combinations, optimise one at a time
        child_node = self.markov_blanket.keys()[0]
        child = self.markov_blanket[child_node]
        K = child.dim[1]
        assert K == len(self.ix), 'problem in dropping factor'

        # use grid search to optimise hyperparameters
        for k in range(K):
            best_i = -1
            best_elbo = -np.Inf
            for i in range(self.n_grid):
                self.ix[k] = i
                elbo = child.calculateELBO_k(k)
                if elbo > best_elbo:
                    best_elbo = elbo
                    best_i = i
                if i == 0:
                    elbo0 = elbo
            self.spatial_sig[k] = best_elbo - elbo0
            self.ix[k] = best_i

    def updateParameters(self):
        self.iter += 1
        if self.iter >= self.start_opt:
            self.optimise()

    # probably already defined
    def calculateELBO(self):
        return 0

    def sample(self, dist='P'):
        ix = s.random.choice(range(1, self.n_grid), self.dim[0], replace=True)
        i0 = s.random.choice(range(self.dim[0]), self.n_diag, replace=False)
        ix[i0] = 0

        self.samp = ix   # saving locally only the grid indices
        return np.array([self.grid_cov[i,:,:] for i in ix])

# TODO add clusters as parameters
class BlockSigmaGrid_Node(SigmaGrid_Node):
    def __init__(self, dim, X, clusters, start_opt=20, n_grid=10, n_diag=0):
        self.clusters = clusters
        super(BlockSigmaGrid_Node, self).__init__(dim, X, start_opt, n_grid, n_diag)

    def compute_cov(self):
        # compute SE block covariance based on clusters
        for i in range(self.n_grid):
            # intialisation of covariance
            self.grid_cov[i,:,:] = BlockSE(self.X, self.clusters, self.l_grid[i])
            self.grid_cov[i,:,:] *= covar_rescaling_factor(self.grid_cov[i,:,:])

            # intialisation of inverse covariance
            self.grid_cov_inv[i,:,:] = BlockInv(self.grid_cov[i,:,:], self.clusters)  # TODO speed up with cholesky

            # initialisation of diagonal terms
            self.grid_cov_inv_diag[i,:] = s.diag(self.grid_cov_inv[i,:,:])