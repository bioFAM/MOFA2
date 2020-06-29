from __future__ import division
import numpy as np
import numpy.ma as ma
from copy import deepcopy
from mofapy2.core import gp_utils
from mofapy2.core import gpu_utils
from mofapy2.core.distributions import *
import scipy as s

# Import manually defined functions
from .variational_nodes import MultivariateGaussian_Unobserved_Variational_Node

# Z_GP_Node_mv
class Z_GP_Node_mv(MultivariateGaussian_Unobserved_Variational_Node):
    """
    Z node with a Multivariate Gaussian prior and variational distribution
    """
    def __init__(self, dim, pmean, pcov, qmean, qcov, weight_views, qE=None):
        super().__init__(dim=dim, pmean = pmean, pcov =pcov, qmean =qmean, qcov =qcov, qE=qE)

        self.mini_batch = None
        self.factors_axis = 1
        self.struct = None
        self.length_scales = None
        self.N = self.dim[0]
        self.weight_views = weight_views

        # Precompute terms (inverse covariance ant its determinant for each factor) to speed up computation
        tmp = self.P.params['cov']
        self.p_cov_inv = np.array([s.linalg.inv(cov) for cov in tmp]) # TODO speed-ups
        self.p_cov_inv_diag = np.array([s.diag(c) for c in self.p_cov_inv]) # TODO speed-ups

    def add_characteristics(self, struct=None, length_scales=None):
        self.struct = np.array(struct)
        self.length_scales = np.array(length_scales)

    def precompute(self, options):
        """ Method to precompute terms to speed up computation """
        gpu_utils.gpu_mode = options['gpu_mode']

    def removeFactors(self, idx, axis=0):
        super().removeFactors(idx, axis)
        self.p_cov_inv = s.delete(self.p_cov_inv, axis=0, obj=idx)
        self.p_cov_inv_diag = s.delete(self.p_cov_inv_diag, axis=0, obj=idx)
        if not self.length_scales is None:
            self.length_scales = s.delete(self.length_scales, obj=idx)
        if not self.struct is None:
            self.struct = s.delete(self.struct, obj=idx)

    def get_mini_batch(self):
        """ Method to fetch minibatch """
        if self.mini_batch is None:
            return self.getExpectations()
        else:
            return self.mini_batch
        
    def updateParameters(self, ix=None, ro=1.):
        
        # Get expectations from other nodes
        W = self.markov_blanket["W"].getExpectations()
        Y = self.markov_blanket["Y"].get_mini_batch()
        tau = self.markov_blanket["Tau"].get_mini_batch()
        mask = [self.markov_blanket["Y"].nodes[m].getMask() for m in range(len(Y))]
      
        if "Sigma" in self.markov_blanket:
            Sigma = self.markov_blanket['Sigma'].get_mini_batch()
            p_cov_inv = Sigma['inv']
            p_cov_inv_diag = Sigma['inv_diag']
        else:
            p_cov_inv = self.p_cov_inv
            p_cov_inv_diag = self.p_cov_inv_diag
            if ix is not None: 
                p_cov_inv = np.array([p_cov_inv[i][ix][:,ix] for i in range(p_cov_inv.shape[0])])
                p_cov_inv_diag = np.array([p_cov_inv_diag[i][ix] for i in range(p_cov_inv_diag.shape[0])])

        # Get variational parameters of current node
        Q = self.Q.getParameters()
        Qmean, Qcov = Q['mean'], Q['cov']
        # TODO implement inducing points --- stochasticVI
        # if ix is not None:
        #     self.mini_batch = {}
        #     Qmean = Qmean[ix,:]
        #     Qvar = Qvar[ix,:]
        
        par_up = self._updateParameters(Y, W, tau, Qmean, Qcov, p_cov_inv, mask)
    
        # Update parameters
        # if ix is None:
        Q['mean'] = par_up['Qmean']
        Q['cov'] = par_up['Qcov']
        # else:
        #     self.mini_batch['E'] = par_up['Qmean']
        #     self.mini_batch['E2'] = s.square(par_up['Qmean']) + par_up['Qvar']

        #     Q['mean'][ix,:] = par_up['Qmean']
        #     Q['var'][ix,:] = par_up['Qvar']

        self.Q.setParameters(mean=Q['mean'], cov=Q['cov'])  # NOTE should not be necessary but safer to keep for now
    
    
    
    def _updateParameters(self, Y, W, tau, Qmean, Qcov, p_cov_inv, mask):
        """ Hidden method to compute parameter updates """

        N = Y[0].shape[0]  # this is different from self.N for minibatch
        M = len(Y)
        K = self.dim[1]

        # Masking
        for m in range(M):
            tau[m][mask[m]] = 0.

        if self.weight_views:
            #scaling factor according to number of features per view m
            weights = [N / Y[m].shape[1] for m in range(M)]
        else:
            weights = [1] * M

        # Precompute terms to speed up GPU computation
        foo = gpu_utils.array(s.zeros((N,K)))
        precomputed_bar = gpu_utils.array(s.zeros((N,K)))
        for m in range(M):
            tau_gpu = gpu_utils.array(tau[m])
            foo += weights[m] * gpu_utils.dot(tau_gpu, gpu_utils.array(W[m]["E2"]))
            bar_tmp1 = gpu_utils.array(W[m]["E"])
            bar_tmp2 = tau_gpu * gpu_utils.array(Y[m])
            precomputed_bar += weights[m] * gpu_utils.dot(bar_tmp2, bar_tmp1)
        foo = gpu_utils.asnumpy(foo)

        # Calculate variational updates
        for k in range(K):
            bar = gpu_utils.array(s.zeros((N,)))
            tmp_cp1 = gpu_utils.array(Qmean[:, s.arange(K) != k])
            for m in range(M):
                tmp_cp2 = gpu_utils.array(W[m]["E"][:, s.arange(K) != k].T)

                bar_tmp1 = gpu_utils.array(W[m]["E"][:,k])
                bar_tmp2 = gpu_utils.array(tau[m])*(-gpu_utils.dot(tmp_cp1, tmp_cp2))

                bar += weights[m] * gpu_utils.dot(bar_tmp2, bar_tmp1)
            bar += precomputed_bar[:,k]
            bar = gpu_utils.asnumpy(bar)

            Qcov[k,:,:] = np.linalg.inv(np.eye(N) * foo[:,k] +  p_cov_inv[k,:,:]) # TODO speed-ups using diagonal structure of first term,
            Qmean[:, k] = gpu_utils.dot(Qcov[k,:,:], bar) 

        # Save updated parameters of the Q distribution
        return {'Qmean': Qmean, 'Qcov':Qcov}

# Elbo calculation per factor - required for grid search for optimal hyperparameter in sigma per factor
    def calculateELBO_k(self, k):
        # Collect parameters and expectations of current node
        Qpar, Qexp = self.Q.getParameters(), self.Q.getExpectations()
        Qmean, Qcov = Qpar['mean'], Qpar['cov']
        
        QE  = Qexp['E']

        if 'Sigma' in self.markov_blanket:
            Sigma = self.markov_blanket['Sigma'].getExpectations()
            p_cov = Sigma['cov']
            p_cov_inv = Sigma['inv']
            p_cov_inv_logdet = Sigma['inv_logdet']
        else:
            p_cov = self.P.params['cov']
            p_cov_inv = self.p_cov_inv
            p_cov_inv_logdet = np.linalg.slogdet(self.p_cov_inv)[1]
            
        # compute term from the exponential in the Gaussian
        tmp1 = -0.5 * (np.trace(gpu_utils.dot(p_cov_inv[k,:,:], Qcov[k,:,:])) +  gpu_utils.dot(QE[:,k].transpose(), gpu_utils.dot(p_cov_inv[k,:,:], QE[:,k])))# expectation of quadratic form

        # compute term from the precision factor in front of the Gaussian
        tmp2 = 0.5 * p_cov_inv_logdet[k]
        lb_p = tmp1 + tmp2
        
        lb_q = -0.5 * np.linalg.slogdet(Qcov[k,:,:])[1] # term -N*(log(2* np.pi)) cancels out between p and q term; -N/2 is added below

        return lb_p - lb_q
    
    # sum up individual EBLO calculations
    def calculateELBO(self):
        elbo = 0
        for k in range(self.dim[1]):
            elbo += self.calculateELBO_k(k)         
        elbo-= self.dim[0] * self.dim[1] / 2. 

        return elbo
    
    def sample(self):
        mu = self.P.params['mean']
        if "Sigma" in self.markov_blanket:
            Sigma = self.markov_blanket['Sigma']
            cov = Sigma.sample()
        else:
            cov = self.P.params['cov']

        samp_tmp = [s.random.multivariate_normal(mu[:,i], cov[i,:,:]) for i in range(self.dim[1])]
        self.samp = s.array([tmp-tmp.mean() for tmp in samp_tmp]).transpose()
        return self.samp
