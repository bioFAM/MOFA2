from __future__ import division
import numpy as np

from mofapy2.core.nodes.variational_nodes import *
from mofapy2.core.gp_utils import *
import scipy as s
import time
from mofapy2.core import gpu_utils


class Sigma_Node(Node):
    pass

class SigmaGrid_Node(Node):
    """
    Sigma node to opitmize the lengthscale parameter for each factor, not a variational node
    This is outside of the variational updates and done by grid search.
    
    PARAMETERS
    ----------
    dim: dimensionality of the node (= number of latent factors)
    sample_cov: covariates for construction of the covariance matrix from distances
    start_opt: in which iteration to start with optimizing the length scale parameters
    n_grid: number of grid points to optimize over
    mv_Znode: whether a multivariate variational is modelled for the node with prior Sigma covariance
    idx_inducing: Index of inducing points (default None - use the full model)
    """
    def __init__(self, dim, sample_cov, start_opt=20, n_grid=10, mv_Znode = False, idx_inducing = None, smooth_all = False):
        super().__init__(dim)
        self.AO = False and mv_Znode #TODO get rid of AO
        self.mini_batch = None
        self.sample_cov = sample_cov
        self.N = sample_cov.shape[0]
        self.start_opt = start_opt
        self.n_grid = n_grid
        self.mv_Znode = mv_Znode
        self.iter = 0                                       # counter of iteration to keep track when to optimize lengthscales
        self.gridix = np.zeros(dim[0], dtype = np.int8)     # index of the grid values to use
        self.struct_sig = np.zeros(dim[0])                  # store improvments compared to diagonal covariance
        self.idx_inducing = idx_inducing
        self.smooth_all = smooth_all
        if not self.idx_inducing is None:
            self.Nu = len(idx_inducing)
        self.compute4init()
        # TODO make this node more memory efficient (avoid loading full covariance matrix into memeory + enable stoch updates - only needs sample_cov and Sigam_U)

    def compute4init(self):
        """
        Function to get the lengthscale grid and covariance matrices on initilisation of the Sigma node
        """
        # get grid of lengthscales
        self.l_grid = get_l_grid(self.sample_cov, n_grid=self.n_grid)
        # add the diagonal covariance (lengthscale 0) to the grid

        if not self.smooth_all:
	        self.l_grid = np.insert(self.l_grid, 0, 0)
	        self.n_grid += 1

        # initialise covariance matrix, inverse and diagonal of inverse
        self.Sigma = np.zeros([self.n_grid, self.N, self.N])
        if self.idx_inducing is None:
            self.Simga_inv = np.zeros([self.n_grid, self.N, self.N])
            self.Simga_inv_diag = np.zeros([self.n_grid, self.N])
            self.Simga_inv_logdet = np.zeros([self.n_grid])
        else:
            self.Simga_inv = np.zeros([self.n_grid, self.Nu, self.Nu])
            self.Simga_inv_diag = np.zeros([self.n_grid, self.Nu])
            self.Simga_inv_logdet = np.zeros([self.n_grid])

        # compute for each lengthscale
        self.compute_cov()
        
    def precompute(self, options):
        gpu_utils.gpu_mode = options['gpu_mode'] # currently not in use for Sigma

    def compute_cov(self):
        """
        Function to compute covariance matrices for all lengthscales
        """
        for i in range(self.n_grid):
            
            # covariance (scaled by Gower factor)
            self.Sigma[i,:,:] = SE(self.sample_cov, self.l_grid[i])
            if not self.mv_Znode:
                self.Sigma[i,:,:] *= covar_rescaling_factor(self.Sigma[i,:,:])

            if self.idx_inducing is None:
                # compute inverse (without Cholesky)
                self.Simga_inv[i,:,:] = s.linalg.inv(self.Sigma[i,:,:])
                # diagonal of inverse
                self.Simga_inv_diag[i, :] = s.diag(self.Simga_inv[i, :, :])
                # determinant of inverse
                self.Simga_inv_logdet[i] = np.linalg.slogdet(self.Simga_inv[i, :, :])[1]

                # compute inverse using Cholesky decomposition (slower) TODO fix
                # L = s.linalg.cholesky(self.Sigma[i,:,:], lower=True)
                # # Li = s.linalg.inv(L)
                # Li = s.linalg.solve_triangular(L, s.eye(self.Sigma.shape[1]), lower = True)
                # self.Simga_inv[i,:,:] = Li.transpose().dot(Li)
                # # determinant of inverse
                # self.Simga_inv_logdet[i] = 2.0 * s.sum(s.log(s.diag(Li)))  # using cholesky decomp and det. of triangular matrix
                # # diagonal of inverse
                # self.Simga_inv_diag[i,:] = s.diag(self.Simga_inv[i,:,:])

            else:
                # compute inverse (without Cholesky)
                self.Simga_inv[i,:,:] = s.linalg.inv(self.Sigma[i][self.idx_inducing,:][:,self.idx_inducing])
                # diagonal of inverse
                self.Simga_inv_diag[i, :] = s.diag(self.Simga_inv[i, :, :])
                # determinant of inverse
                self.Simga_inv_logdet[i] = np.linalg.slogdet(self.Simga_inv[i, :, :])[1]


    # def define_mini_batch(self, ix):
    #     """
    #     Method to define a mini-batch (only for stochastic inference)
    #     """
    #     tmp = self.getExpectations()
    #     # note that the inverse of a subsample will not correspond to the submatrix of the inverse
    #     cov_mini = np.array([tmp['cov'][i][ix][:,ix] for i in range(tmp['cov'].shape[0])])
    #     inv_mini = np.array([tmp['inv'][i][ix][:,ix] for i in range(tmp['inv'].shape[0])])
    #     inv_diag_mini = np.array([tmp['inv_diag'][i][ix] for i in range(tmp['inv'].shape[0])])
    #     self.mini_batch = {'cov' : cov_mini, 'inv': inv_mini, 'inv_diag' : inv_diag_mini, 'E' : cov_mini}
        
    def get_mini_batch(self):
        """ 
        Method to fetch minibatch 
        """
        if self.mini_batch is None:
            return self.getExpectations()
        else:
            return self.mini_batch
    
    def getExpectation(self):
        """ 
        Method to fetch ELBO-optimal covariance matrix per factor (only used upon saving a model) 
        """
        cov = np.array([self.Sigma[i,:,:] for i in self.gridix])
        return cov
    
    def getExpectations(self):
        """ 
        Method to fetch ELBO-optimal covariance matrix, its  inverse and the diagonal of the inverse per factor
        """
        cov = np.array([self.Sigma[i,:,:] for i in self.gridix])
        if not self.AO:
            inv = np.array([self.Simga_inv[i,:,:] for i in self.gridix])
            inv_diag = np.array([self.Simga_inv_diag[i,:] for i in self.gridix])
            cov_inv_logdet = np.array([self.Simga_inv_logdet[i] for i in self.gridix])
            return {'cov':cov, 'inv': inv, 'inv_diag':inv_diag, 'E':cov, 'inv_logdet' : cov_inv_logdet}
        else:
            return{'cov' : cov}
    
    def get_ls(self):
        """ 
        Method to fetch ELBO-optimal length-scales
        """
        ls = np.array([self.l_grid[i] for i in self.gridix])
        return ls

    def getParameters(self):
        """ 
        Method to fetch ELBO-optimal length-scales, improvements compared to diagonal covariance prior and structural positions
        """
        ls = self.get_ls()
        return {'l':ls, 'sig': self.struct_sig, 'sample_cov':self.sample_cov}

    def removeFactors(self, idx, axis=1):
        """
        Method to remove factors 
        """
        self.gridix = s.delete(self.gridix, axis=0, obj=idx)
        self.struct_sig = s.delete(self.struct_sig, axis=0, obj=idx)
        self.updateDim(0, self.dim[0] - len(idx))

    def optimise(self):
        """
        Method to find for each factor the lengthscale parameter that optimises the ELBO of the factor.
        The optimization can be carried out on a per-factor basis (not required on all combinations) as latent variables are independent in the elbo
        """
        if not self.idx_inducing is None:
            var = self.markov_blanket['U']
        else:
            var = self.markov_blanket['Z']
        K = var.dim[1]
        assert K == len(self.gridix), 'problem in dropping factor'

        # use grid search to optimise hyperparameters
        for k in range(K):
            best_i = -1
            best_elbo = -np.Inf
            for i in range(self.n_grid):
                self.gridix[k] = i
                elbo = var.calculateELBO_k(k)
                if elbo > best_elbo:
                    best_elbo = elbo
                    best_i = i
                if i == 0:
                    elbo0 = elbo
            self.struct_sig[k] = best_elbo - elbo0
            self.gridix[k] = best_i
        print('Sigma node has been optimised: Lengthscales =', self.l_grid[self.gridix])
        # print('Sigma node has been optimised: struct_sig =', self.struct_sig)

    def updateParameters(self, ix=None, ro=1.):
        """
        Public method to update the nodes parameters
        Optional arguments for stochastic updates are:
            - ix: list of indices of the minibatch
            - ro: step size of the natural gradient ascent
        Stochastic updates have no effect here yet as ELBO is calculated and optimized on all samples
        """
        self.iter += 1
        if self.iter >= self.start_opt:
            self.optimise()

    def calculateELBO(self):
        return 0