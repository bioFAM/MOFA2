from __future__ import division
import numpy as np

from mofapy2.core.nodes.variational_nodes import *
from mofapy2.core.gp_utils import *
import scipy as s
import time
from mofapy2.core import gpu_utils
import pandas as pd
# from fastdtw import fastdtw
from dtw import dtw # note this is dtw-python not dtw
from scipy.spatial.distance import euclidean
import copy

class Sigma_Node(Node):
    pass


# TODO
#  - make this node more memory efficient for sparse GP (avoid loading full covariance matrix into memeory - only needs sample_cov and Sigam_U)
# - implement warping for more than one covariate
# - do not construct Sigma and Sigma stuff but instead pass U D and zeta directly to the Z node and get rid of Sigma completely?

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
    def __init__(self, dim, sample_cov, groups, start_opt=20, n_grid=10, idx_inducing = None,
                 warping = False, warping_freq = 20, warping_ref = 0, warping_open_begin = True, warping_open_end = True,
                 opt_freq = 10):
        super().__init__(dim)
        self.zeta_opt = False
        self.setscaletoone = False # TODO: set back to True?
        self.mini_batch = None
        self.sample_cov = sample_cov
        self.sample_cov_transformed = copy.copy(sample_cov) # keep original covariate in place
        self.groups = groups
        self.groupsidx = pd.factorize(self.groups)[0]
        self.n_groups = len(np.unique(self.groups))
        self.N = sample_cov.shape[0]
        self.start_opt = start_opt
        self.n_grid = n_grid
        self.mv_Znode = True
        self.iter = 0                                       # counter of iteration to keep track when to optimize lengthscales
        self.n_factors = dim[0]
        self.zeta = np.ones(self.n_factors) * 0.5
        self.gridix = np.zeros(self.n_factors, dtype = np.int8)     # index of the grid values to use per factor
        # self.shift = np.zeros(self.n_groups)                  # group-sepcific offset of covariates
        # self.scaling = np.ones(self.n_groups)                 # group specific scaling of covariates
        self.struct_sig = np.zeros(self.n_factors)                  # store improvments compared to diagonal covariance
        self.idx_inducing = idx_inducing
        self.opt_freq = opt_freq
        self.warping = warping
        assert warping_ref < self.n_groups, "Reference group not correctly specified, exceeds the number of groups."
        self.reference_group = warping_ref
        self.warping_freq = warping_freq
        self.warping_open_begin = warping_open_begin
        self.warping_open_end = warping_open_end
        if not self.idx_inducing is None:
            self.Nu = len(idx_inducing)
        # self.transform_sample_cov()
        self.compute4init()


    def compute4init(self):
        """
        Function to get the lengthscale grid and covariance matrices on initilisation of the Sigma node
        """
        # get grid of lengthscales
        if self.warping:
            idx = np.where(self.groupsidx == self.reference_group)[0]
        else:
            idx = None
        self.l_grid = get_l_grid(self.sample_cov, n_grid=self.n_grid, idx = idx)

        # add the diagonal covariance (lengthscale 0) to the grid
        if not self.zeta_opt:
            # for ls optimization:
            self.l_grid = np.insert(self.l_grid, 0, 0)
            self.n_grid += 1
        else:
            # for zeta optimization:
            self.l_grid = np.insert(self.l_grid, self.n_grid, 0)

        # initialise kernel matrix
        self.K = np.zeros([self.n_grid, self.N, self.N])        # kernel matrix on lengthscale grid

        # initialise covariance matrix
        self.Sigma = np.zeros([self.n_factors, self.N, self.N])

        # initialise inverse and diagonal of inverse and spectral decomposition (only required for inducing points)
        if self.idx_inducing is None:
            self.V = np.zeros([self.n_grid, self.N, self.N])  # eigenvectors of kernel matrix on lengthscale grid
            self.D = np.zeros([self.n_grid, self.N])  # eigenvalues of kernel matrix on lengthscale grid
            self.Sigma_inv = np.zeros([self.n_factors, self.N, self.N])
            self.Sigma_inv_diag = np.zeros([self.n_factors, self.N])
            self.Sigma_inv_logdet = np.zeros([self.n_factors])
        else:
            self.V = np.zeros([self.n_grid, self.Nu, self.Nu])
            self.D = np.zeros([self.n_grid, self.Nu])
            self.Sigma_inv = np.zeros([self.n_factors, self.Nu, self.Nu])
            self.Sigma_inv_diag = np.zeros([self.n_factors, self.Nu])
            self.Sigma_inv_logdet = np.zeros([self.n_factors])

        # compute for each lengthscale the kernel matrix and spectral decomposition
        self.compute_kernel()
        self.compute_cov()

    def precompute(self, options):
        gpu_utils.gpu_mode = options['gpu_mode'] # currently not in use for Sigma

    def compute_kernel(self):
        """
        Function to compute covariance matrices for all lengthscales
        """
        for i in range(self.n_grid):
            self.compute_kernel_at_gridpoint(i)

    def compute_kernel_at_gridpoint(self, i):
        # covariance (scaled by Gower factor)
        self.K[i, :, :] = SE(self.sample_cov_transformed, self.l_grid[i], zeta=0) # zeta is added later
        # self.K[i, :, :] = Cauchy(self.sample_cov_transformed, self.l_grid[i], zeta=0)

        # if not self.mv_Znode:
        #     self.K[i, :, :] *= covar_rescaling_factor(self.K[i, :, :])

        # compute spectral decomposition
        if self.idx_inducing is None:
            self.D[i,:], self.V[i,:,:] = s.linalg.eigh( self.K[i,:,:]) # Sigma = VDV^T with V^T * V = I # important to use eigh and not eig to obtain orthogonal eigenvector (which always exist for symmetric real matrices)
        else:
            self.D[i, :], self.V[i, :, :] = s.linalg.eigh(self.K[i][self.idx_inducing, :][:, self.idx_inducing])

    def compute_cov(self):
        for k in range(self.n_factors):
            self.compute_cov_k(k)

    def compute_cov_k(self,k):
        self.Sigma[k, :, :] = (1 - self.zeta[k]) * self.K[self.gridix[k], :, :] + self.zeta[k] * s.eye(self.N)

        if self.idx_inducing is None:
            # compute inverse using spectral decomposition of K
            if self.zeta[k] != 1:
                self.Sigma_inv[k, :, :] = 1 / (1 - self.zeta[k]) * gpu_utils.dot(gpu_utils.dot(self.V[self.gridix[k], :, :],
                                                                                               s.diag(1 / (self.D[self.gridix[k],:] +
                                                                                                           self.zeta[k] / (1 - self.zeta[k])))),
                                                                                 self.V[self.gridix[k], :, :].transpose())
                self.Sigma_inv_diag[k, :] = s.diag(self.Sigma_inv[k, :, :])
                det_fac = self.N
                self.Sigma_inv_logdet[k] = - det_fac * s.log(1 - self.zeta[k]) - s.log(
                    (self.D[self.gridix[k], :] + self.zeta[k] / (1 - self.zeta[k]))).sum()
            else:
                self.Sigma_inv[k, :, :] = s.eye(self.N)
                self.Sigma_inv_diag[k, :] = s.diag(self.Sigma_inv[k, :, :])
                self.Sigma_inv_logdet[k] = 1

        else:
            if self.zeta[k] != 1:
                self.Sigma_inv[k, :, :] = 1 / (1 - self.zeta[k]) * gpu_utils.dot(gpu_utils.dot(self.V[self.gridix[k], :, :],
                                                                                               s.diag(1 / (self.D[self.gridix[k],:] + self.zeta[k] / (1 - self.zeta[k])))),
                                                                                 self.V[self.gridix[k], :, :].transpose())
                self.Sigma_inv_diag[k, :] = s.diag(self.Sigma_inv[k, :, :])
                det_fac = self.Nu
                self.Sigma_inv_logdet[k] = - det_fac * s.log(1 - self.zeta[k]) -s.log(
                     (self.D[self.gridix[k], :] + self.zeta[k] / (1 - self.zeta[k]))).sum()
            else:
                self.Sigma_inv[k, :, :] = s.eye(self.Nu)
                self.Sigma_inv_diag[k, :] = s.diag(self.Sigma_inv[k, :, :])
                self.Sigma_inv_logdet[k] = 1


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
        Method to fetch ELBO-optimal covariance matrix per factor
        """
        cov = self.Sigma
        return cov
    
    def getExpectations(self):
        """ 
        Method to fetch ELBO-optimal covariance matrix, its  inverse and the diagonal of the inverse per factor
        """
        cov = self.Sigma
        inv = self.Sigma_inv
        inv_diag = self.Sigma_inv_diag
        cov_inv_logdet = self.Sigma_inv_logdet
        return {'cov':cov, 'inv': inv, 'inv_diag':inv_diag, 'E':cov, 'inv_logdet' : cov_inv_logdet}

    def get_ls(self):
        """ 
        Method to fetch ELBO-optimal length-scales
        """
        ls = np.array([self.l_grid[i] for i in self.gridix])
        return ls

    def get_zeta(self):
        """
        Method to fetch ELBO-optimal length-scales
        """
        return self.zeta

    def getParameters(self):
        """ 
        Method to fetch ELBO-optimal length-scales, improvements compared to diagonal covariance prior and structural positions
        """
        ls = self.get_ls()
        scale = 1 - self.get_zeta()
        return {'l':ls, 'scale': scale, 'sig': self.struct_sig, 'sample_cov':self.sample_cov_transformed}

    def removeFactors(self, idx, axis=1):
        """
        Method to remove factors 
        """
        self.gridix = s.delete(self.gridix, axis=0, obj=idx)
        self.zeta = s.delete(self.zeta, axis=0, obj=idx)
        self.Sigma  = s.delete(self.Sigma, axis=0, obj=idx)
        self.Sigma_inv  = s.delete(self.Sigma_inv, axis=0, obj=idx)
        self.Sigma_inv_diag  = s.delete(self.Sigma_inv_diag, axis=0, obj=idx)
        self.Sigma_inv_logdet  = s.delete(self.Sigma_inv_logdet, axis=0, obj=idx)
        self.struct_sig = s.delete(self.struct_sig, axis=0, obj=idx)
        self.updateDim(0, self.dim[0] - len(idx))
        self.n_factors = self.n_factors - 1


    def calc_neg_elbo_k(self, zeta, var, k):
        self.zeta[k] = zeta
        self.compute_cov_k(k)
        elbo = var.calculateELBO_k(k)
        return -elbo

    def optimise(self):
        """
        Method to find for each factor the lengthscale parameter that optimises the ELBO of the factor.
        The optimization can be carried out on a per-factor basis (not required on all combinations) as latent variables are independent in the elbo
        """

        if not self.idx_inducing is None:
            var = self.markov_blanket['U']
            Zvar = self.markov_blanket['U'].markov_blanket['Z']  # use all factor values for alignment #TODO clean this
        else:
            var = self.markov_blanket['Z']
            Zvar = var

        K = var.dim[1]
        assert K == len(self.gridix), 'problem in dropping factor'
        assert K == len(self.zeta), 'problem in dropping factor'
        assert K == self.n_factors, 'problem in dropping factor'

        # perform DTW to align groups
        if self.warping and self.n_groups > 1 and self.iter % self.warping_freq == 0:
            self.align_sample_cov_dtw(Zvar.getExpectation())
            print("Covariates were aligned between groups.")
            # self.align_sample_cov_linear()
            # print("Covariates were aligned between groups: shift:", self.shift, ", scaling:", self.scaling)

        # optimise hyperparamters of GP
        if self.iter % self.opt_freq == 0:
            for k in range(K):
                best_i = -1
                best_zeta = -1
                best_elbo = -np.Inf
                # use grid search to optimise lengthscale hyperparameters
                for i in range(self.n_grid):
                    self.gridix[k] = i
                    res = s.optimize.minimize(self.calc_neg_elbo_k, args=(var, k), x0 = self.zeta[k], bounds=[(1e-10, 1-1e-10)]) # L-BFGS-B
                    elbo = -res.fun
                    # print("ELBO", elbo, ", i:", self.gridix[k], ", zeta:", res.x)
                    if elbo > best_elbo:
                        best_elbo = elbo
                        best_i = i
                        best_zeta = res.x
                    if i == 0 and not self.zeta_opt and self.setscaletoone: # for i = 0 (lengthscale of 0)
                        elbo0 = elbo
                        best_zeta = 1 # TODO this is not the same (striped vs diagonal)
                    if self.zeta_opt:
                        if best_zeta > 1-1e-7: # for best_zeta = 1, lengthscale is irrelevant (set to zero)
                            best_zeta = 1
                            best_i = self.n_grid # last grid point at position self.n_grid + 1 corresponds to zero (not included in K grid as 0 * K)

                # self.struct_sig[k] = best_elbo - elbo0
                self.struct_sig[k] = np.nan # could measure best_elbo - calc_neg_elbo_k(1, var, k)
                self.gridix[k] = best_i
                self.zeta[k] = best_zeta

            self.compute_cov()
            print('Sigma node has been optimised: Lengthscales =', self.l_grid[self.gridix], ', Scales =',  1- self.zeta)
            # print('Sigma node has been optimised: struct_sig =', self.struct_sig)


    def align_sample_cov_dtw(self, Z):
        """
        Method to perform DTW between groups in the factor space
        """

        paths = []
        for g in range(self.n_groups):
            if g is not self.reference_group:
                # reorder by covariate value to ensure monotonicity constrains are correctly placed
                idx_ref_order = np.argsort(self.sample_cov[self.groupsidx == self.reference_group,0])
                idx_query_order = np.argsort(self.sample_cov[self.groupsidx == g,0])
                # allow for partial matching(no corresponding end and beginning)
                step_pattern = "asymmetric" if self.warping_open_begin or self.warping_open_end else "symmetric2"
                alignment = dtw(Z[self.groupsidx == g, :][idx_query_order,:], Z[self.groupsidx == self.reference_group, :][idx_ref_order,:],
                                open_begin = self.warping_open_begin, open_end = self.warping_open_end, step_pattern=step_pattern)
                query_idx = alignment.index1 # dtw-python
                ref_idx = alignment.index2
                # alignment = dtw(Z[self.groupsidx == g, :][idx_query_order,:], Z[self.groupsidx == self.reference_group,:][idx_ref_order,:], dist = euclidean) #dtw
                # query_idx = alignment[3][0] #dtw
                # ref_idx = alignment[3][1]
                ref_val = self.sample_cov[self.groupsidx == self.reference_group, 0][idx_ref_order][ref_idx]
                idx = np.where(self.groupsidx == g)[0][idx_query_order][query_idx]
                self.sample_cov_transformed[idx, 0] = ref_val

        # recompute covariance matrices
        self.compute_kernel()
        self.compute_cov()

                # distance, path = fastdtw(Z[self.groupsidx == self.reference_group,:][idx_ref_order,:],  Z[self.groupsidx == g, :][idx_query_order,:], dist=euclidean) # fastdtw
                # for i in range(len(path)):
                #     ref_val = self.sample_cov[self.groupsidx == self.reference_group, :][idx_ref_order,:][path[i][0]]
                #     idx = np.where(self.groupsidx == g)[0][idx_query_order][path[i][1]]
                #     self.sample_cov_transformed[idx, : ] = ref_val


    # def transform_sample_cov_linear(self):
    #     self.sample_cov_transformed =  self.scaling[self.groupsidx,None] * self.sample_cov + self.shift[self.groupsidx, None]


    # def align_sample_cov_linear(self):
    #     """
    #     Method to find for a linear transformation of covariates to align across groups by optimising the ELBO
    #     """
    #     start_val = np.hstack([self.scaling[1:], self.shift[1:]])
    #     bounds = np.tile([0, None, None, None], self.n_groups - 1)
    #     res = s.optimize.minimize(self.calc_ELBO_in_warping, start_val,
    #                               bounds= s.optimize.Bounds(lb = np.repeat([0, -np.inf],self.n_groups - 1), ub =np.repeat(np.inf,2*self.n_groups - 2)))
    #
    #     self.scaling = np.insert(res['x'][:self.n_groups -1], 0, 1)
    #     self.shift = np.insert(res['x'][self.n_groups-1:], 0, 0)
    #
    #     # if self.iter == 20: # for debugging purposes - plot objective function
    #     #     import matplotlib.pyplot as plt
    #     #     plt.figure(5)
    #     #     a = np.linspace(0.05,5)
    #     #     b = np.linspace(-5,5)
    #     #     l = np.zeros([len(a), len(b)])
    #     #     for i in range(len(a)):
    #     #         for j in range(len(b)):
    #     #             l[i, j] = self.calc_ELBO_in_warping([a[i], b[j]])
    #     #     aa, bb = np.meshgrid(a, b)
    #     #     plt.pcolormesh(aa, bb, l.transpose())
    #     #     plt.xlabel("scale")
    #     #     plt.ylabel("shift")
    #     #     plt.axvline(x=0.2, color = "black")
    #     #     plt.axhline(y=0, color = "black")
    #     #     plt.axvline(x=self.scaling[1], color="red")
    #     #     plt.axhline(y=self.shift[1], color="red")
    #     #     plt.colorbar()
    #     #     plt.show()
    #
    #     self.transform_sample_cov() # transform sample cov
    #     self.compute_cov()  # recalculate all Sigma matrices in grid
    #
    # def calc_ELBO_in_warping(self, x):
    #     """"
    #     Method to calculate the ELBO terms in the 2 warping parameters per group
    #     """
    #     a = x[:self.n_groups -1]
    #     b = x[self.n_groups-1:]
    #
    #     if not self.idx_inducing is None:
    #         var = self.markov_blanket['U']
    #     else:
    #         var = self.markov_blanket['Z']
    #     self.scaling = np.insert(a,0,1)
    #     self.shift = np.insert(b,0,0)
    #     self.transform_sample_cov()
    #     for i in np.unique(self.gridix):
    #         self.compute_cov_at_gridpoint(i) # recompute covariance matrix at current lengthscale parameters (might differ per factor)
    #     elb = var.calculateELBO()
    #
    #     return -elb


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