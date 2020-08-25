from __future__ import division
import numpy as np

from mofapy2.core.nodes.variational_nodes import *
from mofapy2.core.nodes.Kc_node import Kc_Node
from mofapy2.core.nodes.Kg_node import Kg_Node
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

     #TODO set option to makse K_g a dense matrix of ones (all groups connected)
class kronSigma_Node(Node):
    """
    Sigma node to optimises the GP hyperparameters parameter for each factor and
    perform alignment of covariates per group.
    
    PARAMETERS
    ----------
    dim: dimensionality of the node (= number of latent factors)
    sample_cov: covariates for construction of the covariance matrix from distances (array of length G x C)
    groups: group label of each observation (array of length G)
    start_opt: in which iteration to start with optimizing the GP hyperparameters
    n_grid: number of grid points to optimize the lengthscale on
    idx_inducing: Index of inducing points (default None - models the full kernel matrix)
    warping: Boolean, whether to perform warping of covariates across groups in the latent space
    warping_freq: how often to perform wraping, every n-th iteration
    warping_ref: reference group for warping
    warping_open_begin: Allow free beginning for the warped covariates in each group?
    warping_open_end:  Allow free end for the warped covariates in each group?
    opt_freq: how often to hyperparameter optimization, every n-th iteration
    """
    def __init__(self, dim, sample_cov, groups, start_opt=20, n_grid=10, idx_inducing = None,
                 warping = False, warping_freq = 20, warping_ref = 0, warping_open_begin = True,
                 warping_open_end = True, opt_freq = 10, rankx = 2):
        super().__init__(dim)

        # dimensions and inputs
        self.sample_cov = sample_cov
        self.sample_cov_transformed = copy.copy(sample_cov)         # keep original covariate in place
        self.group_labels = groups

        # covariate kernel is initliazed after first wraping
        self.N = sample_cov.shape[0]                                # total number of observation (C*G in the complete case, otherwise less)
        self.K = dim[0]                                             # number of factors

        # hyperparameter optimization
        self.start_opt = start_opt
        self.n_grid = n_grid
        self.iter = 0                                               # counter of iteration to keep track when to optimize lengthscales
        self.zeta = np.ones(self.K) * 0.5                           # noise hyperparameter
        self.gridix = np.zeros(self.K, dtype = np.int8)             # index of the lengthscale grid values to use per factor
        self.struct_sig = np.zeros(self.K)                          # store ELBO improvements compared to diagonal covariance
        self.opt_freq = opt_freq
        assert self.start_opt % self.opt_freq == 0,\
            "start_opt should be a multiple of opt_freq"  # to ensure in the first opt. step optimization is performed

        # initalize group kernel
        self.initKg(self.group_labels, rankx)

        # initialize covariate kernel if no wraping, otherwise after first warping
        if not warping:
            self.initKc(self.sample_cov_transformed)

        # warping
        self.warping = warping
        assert warping_ref < self.G,\
            "Reference group not correctly specified, exceeds the number of groups."
        self.reference_group = warping_ref
        self.warping_freq = warping_freq
        self.warping_open_begin = warping_open_begin
        self.warping_open_end = warping_open_end
        assert self.warping_freq % self.opt_freq == 0,\
            "start_opt should be a multiple of opt_freq" # to ensure in the first opt. step alignment is performed

        # sparse GPs
        self.idx_inducing = idx_inducing
        if not self.idx_inducing is None:
            self.Nu = len(idx_inducing)

        # initialize Sigma terms (unstructured)
        self.Sigma_inv = np.zeros([self.K, self.N, self.N])
        for k in range(self.K):
            self.Sigma_inv[k, :, :] = np.eye(self.N)
        self.Sigma_inv_logdet = np.ones(self.K)


    def initKc(self, transformed_sample_cov):
        """
        Method to initialize the components required for the covariate kernel
        """
        self.covariates = np.unique(transformed_sample_cov, axis=0)  # distinct covariate values
        self.covidx = [np.where((self.covariates == transformed_sample_cov[j, :]).all(axis=1))
                       for j in range(self.N)]  # for each sample gives the idx in covariates
        self.C = self.covariates.shape[0]  # number of distinct covariate values

        # set covariate kernel
        self.Kc = Kc_Node(dim=(self.K, self.C), covariates = self.covariates, n_grid = self.n_grid)


    def initKg(self, group_labels, rank):
        """
        Method to initialize the components required for the group kernel
        """
        self.groupsidx = pd.factorize(group_labels)[0]                    # for each sample gives the idx in groups
        self.groups = np.unique(group_labels)                             # distinct group labels
        self.G = len(self.groups)                                         # number of groups

        # set group kernel
        self.Kg = Kg_Node(dim=(self.K, self.G), groups = self.groups, rank = rank)


    def precompute(self, options):
        gpu_utils.gpu_mode = options['gpu_mode'] # TODO: implement, this is currently not in use

    
    def get_components(self,k):
        """
        Method to fetch ELBO-optimal covariance matrix components for a given factor k
        """
        Vc, Dc = self.Kc.get_kernel_components_k(k)
        Vg, Dg = self.Kg.get_kernel_components_k(k)
        zeta = self.get_zeta()[k]

        return {'Vc':Vc, 'Vg': Vg, 'Dc':Dc, 'Dg':Dg, 'zeta' : zeta}

    def calc_sigma_inverse(self):
        """
        Method to compute the inverse of sigma and its log determinant based on the spectral decomposition
         of the kernel matrices for all factors
        """
        for k in range(self.K):
            self.calc_sigma_inverse_k(k)

    def calc_sigma_inverse_k(self, k):
        """
        Method to compute the inverse of sigma and its log determinant based on the spectral decomposition
         of the kernel matrices for a given factor k
        """
        assert self.N == self.G * self.C and (np.bincount(self.groupsidx) == self.C).all(), \
            "Data has no kronecker structure - fill in missing values and shut off alignment" # TODO
        components = self.get_components(k)
        term1 = gpu_utils.dot(np.kron(components['Vc'], np.eye(self.G)), np.kron(np.eye(self.C), components['Vg']))
        term2diag = 1/ (np.repeat(components['Dc'], self.G) * np.tile(components['Dg'], self.C) + self.zeta[k] / (1-self.zeta[k]))
        assert (np.diag(term2diag) == np.linalg.inv(np.kron(np.diag(components['Dc']), np.diag(components['Dg'])) + self.zeta[k]/(1-self.zeta[k]) * np.eye(self.C * self.G))).all(),\
            "Calculation of diagonal went wrong..."
        term3 = gpu_utils.dot(np.kron(np.eye(self.C), components['Vg'].transpose()), np.kron(components['Vc'].transpose(), np.eye(self.G)))

        self.Sigma_inv[k, :, :] = 1 / (1-self.zeta[k]) * gpu_utils.dot(term1, gpu_utils.dot(np.diag(term2diag), term3))
        self.Sigma_inv_logdet[k] =  - self.N * s.log(1 - self.zeta[k]) + s.log(term2diag).sum()

        # TODO non-Kronecker strcutre - some gourps can have multiple samples with same covariate value, some groups might be missing samples with covariate value

    def get_sigma_inv_Terms_k(self, k):
        """ 
         Method to fetch ELBO-optimal inverse covariance matrix and its determinant for a given factor k
        """
        return {'inv': self.Sigma_inv[k,:,:], 'inv_logdet':  self.Sigma_inv_logdet[k]}


    def get_sigma_inv_Terms(self):
        """ 
        Method to fetch ELBO-optimal inverse covariance matrix and its determinant for all factors
        """
        return {'inv': self.Sigma_inv, 'inv_logdet': self.Sigma_inv_logdet}


    def get_ls(self):
        """
        Method to fetch ELBO-optimal length-scales
        """
        ls = self.Kc.get_ls()
        return ls


    def get_x(self):
        """
        Method to get low rank group covariance matrix
        """
        x = self.Kg.get_x()
        return x

    def get_sigma(self):
        """
        Method to get sigma hyperparametr
        """
        sigma = self.Kg.get_sigma()
        return sigma

    def getExpectation(self):
        """
        Method to get sigma hyperparametr
        """
        Sigma = np.zeros([self.K, self.N, self.N])
        assert self.N == self.G * self.C and (np.bincount(self.groupsidx) == self.C).all(), \
            "Data has no kronecker structure - fill in missing values and shut off alignment" # TODO
        for k in range(self.K):
            components = self.get_components(k)
            term1 = gpu_utils.dot(np.kron(components['Vc'], np.eye(self.G)), np.kron(np.eye(self.C), components['Vg']))
            term2diag = (np.repeat(components['Dc'], self.G) * np.tile(components['Dg'], self.C) + self.zeta[k])
            term3 = gpu_utils.dot(np.kron(np.eye(self.C), components['Vg'].transpose()), np.kron(components['Vc'].transpose(), np.eye(self.G)))

            Sigma[k, :, :] = gpu_utils.dot(term1, gpu_utils.dot(np.diag(term2diag), term3))

        return Sigma

    def get_zeta(self):
        """
        Method to fetch noise parameter
        """
        return self.zeta


    def getParameters(self):
        """ 
        Method to fetch ELBO-optimal length-scales, improvements compared to diagonal covariance prior and structural positions
        """
        ls = self.get_ls()
        zeta = self.get_zeta()
        x = self.get_x()
        sigma = self.get_sigma()
        return {'l':ls, 'scale': 1-zeta, 'sig': self.struct_sig, 'sample_cov':self.sample_cov_transformed,  'x': x, 'sigma' : sigma}


    def removeFactors(self, idx, axis=1):
        """
        Method to remove factors 
        """
        self.Kg.removeFactors()
        self.Kc.removeFactors()
        self.zeta = s.delete(self.zeta, axis=0, obj=idx)
        self.struct_sig = s.delete(self.struct_sig, axis=0, obj=idx)
        self.updateDim(0, self.dim[0] - len(idx))
        self.K = self.K - 1


    def calc_neg_elbo_k(self, par, lidx, k, var):
        self.zeta[k] = par[0]
        self.Kc.set_gridix(lidx, k)
        sigma = par[1]
        x = par[2:]
        assert len(x) == self.Kg.rank * self.G,\
            "Length of x incorrect: Is %s, should be  %s * %s" % (len(x), self.Kg.rank, self.G)
        x = x.reshape(self.Kg.rank, self.G)
        self.Kg.set_parameters(x=x, sigma=sigma, k=k)
        self.calc_sigma_inverse_k(k)
        elbo = var.calculateELBO_k(k)
        return -elbo


    def optimise(self):
        """
        Method to find for each factor the lengthscale parameter that optimises the ELBO of the factor.
        The optimization can be carried out on a per-factor basis (not required on all combinations) as latent variables are independent in the elbo
        """

        # get Z node of Markov blanket
        Zvar = self.markov_blanket['Z']
        ZE = Zvar.getExpectation()
        K = Zvar.dim[1]
        assert K == len(self.get_ls()) and K == len(self.zeta) and K == self.K,\
            'problem in dropping factor'

        # perform DTW to align groups
        if self.warping and self.G > 1 and self.iter % self.warping_freq == 0:
            self.align_sample_cov_dtw(ZE)
            print("Covariates were aligned between groups.")

        # optimise hyperparamters of GP
        if self.iter % self.opt_freq == 0:
            for k in range(K):

                best_i = -1
                best_zeta = -1
                best_elbo = -np.Inf
                par_init = np.hstack([self.zeta[k], self.get_sigma()[k], self.get_x()[k,:,:].flatten()])

                # use grid search to optimise lengthscale hyperparameters
                for lidx in range(self.n_grid):
                    # par = (zeta, x, lidx), loop over lidx
                    bounds = [(1e-10, 1-1e-10), (1e-10, 1-1e-10)]
                    bounds = bounds + [(-1,1)] * self.G *  self.Kg.rank
                    # make sure initial parameters are in admissible region
                    par_init = np.max(np.vstack([par_init, [bounds[k][0] for k in range(len(bounds))]]), axis = 0)
                    par_init = np.min(np.vstack([par_init, [bounds[k][0] for k in range(len(bounds))]]), axis = 0)

                    # [bounds[k][0] for k in range(len(bounds))]
                    # [bounds[k][1] for k in range(len(bounds))]
                    res = s.optimize.minimize(self.calc_neg_elbo_k, args=(lidx, k, Zvar), x0 = par_init, bounds=bounds) # L-BFGS-B # TODO suitable bounds?
                    elbo = -res.fun
                    # print("ELBO", elbo, ", i:", self.gridix[k], ", zeta:", res.x)
                    if elbo > best_elbo:
                        best_elbo = elbo
                        best_lidx = lidx
                        best_zeta = res.x[0]
                        best_sigma = res.x[1]
                        best_x = res.x[2:].reshape(self.Kg.rank, self.G)

                # save optimized kernel paramters
                self.struct_sig[k] = np.nan # TODO: sth like best_elbo - self.calc_neg_elbo_k(var, k, zeta = 0) but no 1-a to a interparationa any more as first term can have any size
                self.Kc.set_gridix(best_lidx, k)
                self.Kg.set_parameters(x=best_x, sigma=best_sigma, k=k)
                self.zeta[k] = best_zeta

            self.calc_sigma_inverse()
            print('Sigma node has been optimised: Lengthscales =', self.get_ls(), ', Scale =',  1-self.get_zeta(), ', sigma =',  self.get_sigma())
            for k in range(self.K):
                print('Sigma node has been optimised: x_%s =' % k , np.dot(self.get_x()[k,:,:].transpose(), self.get_x()[k,:,:]))
            # print('Sigma node has been optimised: Smoothness =', self.struct_sig)


    def align_sample_cov_dtw(self, Z):
        """
        Method to perform DTW between groups in the factor space.
        The set of possible values for covaraites cannot be expaned (all need to be contained in the reference group)
        Thus, this does not requrie an update of Kc but only of indices mapping samples to covariates.
        """
        paths = []
        for g in range(self.G):
            if g is not self.reference_group:
                # reorder by covariate value to ensure monotonicity constrains are correctly placed
                idx_ref_order = np.argsort(self.sample_cov[self.groupsidx == self.reference_group,0])
                idx_query_order = np.argsort(self.sample_cov[self.groupsidx == g,0])
                # allow for partial matching (no corresponding end and beginning)
                step_pattern = "asymmetric" if self.warping_open_begin or self.warping_open_end else "symmetric2"
                alignment = dtw(Z[self.groupsidx == g, :][idx_query_order,:], Z[self.groupsidx == self.reference_group, :][idx_ref_order,:],
                                open_begin = self.warping_open_begin, open_end = self.warping_open_end, step_pattern=step_pattern)
                query_idx = alignment.index1 # dtw-python
                ref_idx = alignment.index2
                ref_val = self.sample_cov[self.groupsidx == self.reference_group, 0][idx_ref_order][ref_idx]
                idx = np.where(self.groupsidx == g)[0][idx_query_order][query_idx]
                self.sample_cov_transformed[idx, 0] = ref_val

        # adapt covaraite kernel to warped covariates
        if self.iter == self.start_opt:
            self.initKc(self.sample_cov_transformed)
        else:
            # only adapt indices
            self.covidx = [np.where((self.covariates == self.sample_cov_transformed[j, :]).all(axis=1))
                        for j in range(self.N)]  # for each sample gives the idx in covariates

        # TODO need to restore kroncker strucutre here

    def updateParameters(self, ix, ro):
        """
        Public method to update the nodes parameters
        """
        self.iter += 1
        if self.iter >= self.start_opt:
            self.optimise()

    def calculateELBO(self): # no contribution to ELBO
        return 0