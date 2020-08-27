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
from dtw import dtw  # note this is dtw-python not dtw
from scipy.spatial.distance import euclidean
import copy
import gpytorch
import torch

class Gpytorch_Sigma_Node(Node):
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
                 warping_open_end = True, opt_freq = 10, rankx = 2, torch_seed = 2346):
        super().__init__(dim)

        # dimensions and inputs
        self.sample_cov = sample_cov
        self.sample_cov_transformed = copy.copy(sample_cov)         # keep original covariate in place
        self.group_labels = groups

        # covariate kernel is initliazed after first wraping
        self.N = sample_cov.shape[0]                                # total number of observation (C*G in the complete case, otherwise less)
        self.K = dim[0]                                             # number of factors

        # groups
        self.groupsidx = pd.factorize(self.group_labels)[0]                    # for each sample gives the idx in groups
        self.groups = np.unique(self.group_labels)                             # distinct group labels
        self.G = len(self.groups)                                              # number of groups
        self.rankx = rankx

        # hyperparameter optimization
        self.start_opt = start_opt
        self.n_grid = n_grid
        self.iter = 0                                               # counter of iteration to keep track when to optimize lengthscales
        likelihood = gpytorch.likelihoods.MultitaskGaussianLikelihood(num_tasks=2)
        self.struct_sig = np.zeros(self.K)                          # store ELBO improvements compared to diagonal covariance
        self.opt_freq = opt_freq
        assert self.start_opt % self.opt_freq == 0, \
            "start_opt should be a multiple of opt_freq"  # to ensure in the first opt. step optimization is performed

        # warping
        self.warping = warping
        assert warping_ref < self.G, \
            "Reference group not correctly specified, exceeds the number of groups."
        self.reference_group = warping_ref
        self.warping_freq = warping_freq
        self.warping_open_begin = warping_open_begin
        self.warping_open_end = warping_open_end
        assert self.warping_freq % self.opt_freq == 0, \
            "start_opt should be a multiple of opt_freq" # to ensure in the first opt. step alignment is performed

        # sparse GPs
        self.idx_inducing = idx_inducing
        if not self.idx_inducing is None:
            self.Nu = len(idx_inducing)

        # initialize Sigma terms (unstructured) (corresponds to post hoc fitting of GPs)
        self.Sigma =  np.zeros([self.K, self.N, self.N])
        self.Sigma_inv = np.zeros([self.K, self.N, self.N])
        for k in range(self.K):
            self.Sigma[k, :, :] = np.eye(self.N)
            self.Sigma_inv[k, :, :] = np.eye(self.N)
        self.Sigma_inv_logdet = np.ones(self.K)

        # initialize list of GP per factor
        self.gp = [None] * self.K
        self.likelihood = [gpytorch.likelihoods.MultitaskGaussianLikelihood(num_tasks=self.G)] * self.K
        self.noise = [np.nan] * self.K
        self.sigma = [np.nan] * self.K
        self.ls = [np.nan] * self.K
        self.Gmat = [np.nan] * self.K
        torch.manual_seed(torch_seed)

    def optimise(self, gp_iter = 300):
        """
        Method to find for each factor the lengthscale parameter that optimises the ELBO of the factor.
        The optimization can be carried out on a per-factor basis (not required on all combinations) as latent variables are independent in the elbo
        """

        # get Z node of Markov blanket
        Zvar = self.markov_blanket['Z']
        ZE = Zvar.getExpectation()
        K = Zvar.dim[1]
        assert K == len(self.get_ls()) and K == len(self.noise) and K == self.K, \
            'problem in dropping factor'

        # perform DTW to align groups
        if self.warping and self.G > 1 and self.iter % self.warping_freq == 0:
            self.align_sample_cov_dtw(ZE)
            print("Covariates were aligned between groups.")

        # optimise hyperparamters of GP
        ZE = torch.as_tensor(ZE, dtype=torch.float32)
        if self.iter % self.opt_freq == 0:
            for k in range(K):
                ytrain = torch.stack([ZE[self.groupsidx == g,k] for g in range(self.G)])
                xtrain = torch.as_tensor(self.sample_cov_transformed[self.groupsidx == 0], dtype=torch.float32) # TODO only works for kroncker structure
                self.gp[k] = MultitaskGPModel(train_x =xtrain, train_y= ytrain, likelihood=self.likelihood[k], # TODO avoid entire re-init, do hot-start
                                              n_tasks=self.G, rank_x=self.rankx) # TODO this does not respect the posterior uncertainty of Z add this ot the lokelihood
                # set bounds on diagonal of task kernel
                # self.gp[k].covar_module.task_covar_module.register_constraint("raw_var", gpytorch.constraints.Interval(0.0,1.0))
                # print(f'Variance constraint: {self.gp[k].covar_module.task_covar_module.raw_var_constraint}')

                training_iterations = gp_iter

                # Find optimal model hyperparameters
                self.gp[k].train()
                self.likelihood[k].train()

                # Use the adam optimizer
                optimizer = torch.optim.Adam([
                    {'params': self.gp[k].parameters()},  # Includes GaussianLikelihood parameters
                ], lr=0.1)

                # "Loss" for GPs - the marginal log likelihood
                mll = gpytorch.mlls.ExactMarginalLogLikelihood(self.likelihood[k], self.gp[k])

                for i in range(training_iterations):
                    optimizer.zero_grad()
                    output = self.gp[k](xtrain)
                    loss = -mll(output, ytrain)
                    loss.backward()
                    print('Iter %d/%d - Loss: %.3f' % (i + 1, training_iterations, loss.item()))
                    optimizer.step()

                self.gp[k].eval()
                self.likelihood[k].eval()

                print('Sigma node for factor %s has been optimised: ' %k)
                # Show fitted parameters for one example model (note that these are transformed for the actual values):
                self.noise[k] = (self.gp[k].likelihood.noise).detach().numpy()
                self.ls[k] = (self.gp[k].covar_module.data_covar_module.lengthscale).detach().numpy()
                self.sigma[k] = (self.gp[k].covar_module.task_covar_module.var).detach().numpy()
                B = (self.gp[k].covar_module.task_covar_module.covar_factor).detach().numpy()
                self.Gmat[k] = np.dot(B, B.transpose()) + np.diag(self.sigma[k])
            print('Sigma node has been optimised: Lengthscales =', self.get_ls(), ', noise =', self.get_noise())
            for k in range(self.K):
                print('Sigma node has been optimised: Gmat_%s =' % k, self.Gmat[k])

            self.calc_sigma_terms() #TODO avoid recomputing this, extract form object?

    def make_predictions(self, k, test_x):
        # Make predictions
        with torch.no_grad(), gpytorch.settings.fast_pred_var():
            predictions = self.likelihood[k](self.gp[k](test_x))
            mean = predictions.mean
            lower, upper = predictions.confidence_region()

        return mean, lower, upper

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
                idx_ref_order = np.argsort(self.sample_cov[self.groupsidx == self.reference_group, 0])
                idx_query_order = np.argsort(self.sample_cov[self.groupsidx == g, 0])
                # allow for partial matching (no corresponding end and beginning)
                step_pattern = "asymmetric" if self.warping_open_begin or self.warping_open_end else "symmetric2"
                alignment = dtw(Z[self.groupsidx == g, :][idx_query_order, :],
                                Z[self.groupsidx == self.reference_group, :][idx_ref_order, :],
                                open_begin=self.warping_open_begin, open_end=self.warping_open_end,
                                step_pattern=step_pattern)
                query_idx = alignment.index1  # dtw-python
                ref_idx = alignment.index2
                ref_val = self.sample_cov[self.groupsidx == self.reference_group, 0][idx_ref_order][ref_idx]
                idx = np.where(self.groupsidx == g)[0][idx_query_order][query_idx]
                self.sample_cov_transformed[idx, 0] = ref_val

        # adapt covaraite kernel to warped covariates
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

    def calculateELBO(self):  # no contribution to ELBO
        return 0


    def calc_sigma_terms(self):
        """
         Method to save inverse covariance matrix and its determinant for a given factor k
        """
        for k in range(self.K):
            # output = self.gp[k](torch.as_tensor(self.sample_cov_transformed[self.groupsidx == 0], dtype=torch.float32)) # TODO this is the covariance of the predicitve distribution not of the marginal!
            # sigma_mixed = output.covariance_matrix.detach().numpy()
            # self.Sigma[k] = np.vstack([np.hstack(
            #     [sigma_mixed[np.arange(self.N) % self.G == g][:, np.arange(self.N) % self.G == h]
            #     for h in range(self.G)]
            # ) for g in range(self.G)])

            # todo: doublecheck - is there a direct way to obtain the maringal covariance from the fitted gp model
            KT = SE(self.sample_cov_transformed[self.groupsidx == 0], self.ls[k], zeta=0)
            self.Sigma[k] = np.kron(self.Gmat[k], KT) + self.noise[k] * np.eye(self.N)

            # uncomment the following to propagate back to Z updates
            self.Sigma_inv[k] = np.linalg.inv(self.Sigma[k]) # TODO re-use stuff from gp model
            self.Sigma_inv_logdet[k] = np.linalg.slogdet(self.Sigma_inv[k])[1] # TODO re-use stuff from gp model


    def get_sigma_inv_Terms_k(self, k):
        """
         Method to fetch ELBO-optimal inverse covariance matrix and its determinant for a given factor k
        """
        return {'inv': self.Sigma_inv[k, :, :], 'inv_logdet': self.Sigma_inv_logdet[k]}

    def get_sigma_inv_Terms(self):
        """
        Method to fetch ELBO-optimal inverse covariance matrix and its determinant for all factors
        """
        return {'inv': self.Sigma_inv, 'inv_logdet': self.Sigma_inv_logdet}

    def getExpectation(self):
        """
        Method to get covariance matrix
        """
        return self.Sigma

    def get_ls(self):
        """
        Method to fetch ELBO-optimal length-scales
        """
        return self.ls


    def get_Gmat(self):
        """
        Method to get sigma hyperparametr
        """
        return self.Gmat

    def get_noise(self):
        """
        Method to fetch noise parameter
        """
        return self.noise

    def getParameters(self):
        """
        Method to fetch ELBO-optimal length-scales, improvements compared to diagonal covariance prior and structural positions
        """
        ls = self.get_ls()
        noise = self.get_noise()
        Gmat = self.get_Gmat()
        return {'l': ls, 'noise': noise, 'sig': self.struct_sig, 'sample_cov': self.sample_cov_transformed, 'Gmat': Gmat}

    def removeFactors(self, idx, axis=1):
        """
        Method to remove factors
        """
        self.noise = s.delete(self.zeta, axis=0, obj=idx)
        self.ls = s.delete(self.ls, axis=0, obj=idx)
        self.Gmat = s.delete(self.Gmat, axis=0, obj=idx)
        self.struct_sig = s.delete(self.struct_sig, axis=0, obj=idx)
        self.updateDim(0, self.dim[0] - len(idx))
        self.K = self.K - 1
        self.gp = s.delete(self.gp, axis=0, obj=idx)
        self.likelihood = s.delete(self.likelihood, axis=0, obj=idx)

