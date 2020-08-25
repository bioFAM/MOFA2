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
class Z_GP_Node_mv_kron(MultivariateGaussian_Unobserved_Variational_Node):
    """
    Z node with a Multivariate Gaussian prior and variational distribution
    """
    def __init__(self, dim, pmean, pcov, qmean, qcov, qE=None, weight_views = False):
        super().__init__(dim=dim, pmean = pmean, pcov =pcov, qmean =qmean, qcov =qcov, qE=qE)

        self.mini_batch = None
        self.factors_axis = 1
        self.struct = None
        self.length_scales = None
        self.N = self.dim[0]
        self.K = self.dim[1]
        self.weight_views = weight_views

    def precompute(self, options):
        """ Method to precompute terms to speed up computation """
        gpu_utils.gpu_mode = options['gpu_mode']

    def removeFactors(self, idx, axis=0):
        super().removeFactors(idx, axis)
        self.K = self.dim[1]

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

        assert "Sigma" in self.markov_blanket,\
            "Sigma node needs to be in Markov blanket if using Z-GP node"
        Sigma = self.markov_blanket['Sigma'].get_sigma_inv_Terms()
        p_cov_inv = Sigma['inv']

        # Get variational parameters of current node
        Q = self.Q.getParameters()
        Qmean, Qcov = Q['mean'], Q['cov']
        
        par_up = self._updateParameters(Y, W, tau, Qmean, Qcov, p_cov_inv, mask)
    
        # Update parameters
        Q['mean'] = par_up['Qmean']
        Q['cov'] = par_up['Qcov']

        self.Q.setParameters(mean=Q['mean'], cov=Q['cov'])  # NOTE should not be necessary but safer to keep for now

    
    def _updateParameters(self, Y, W, tau, Qmean, Qcov, p_cov_inv, mask):
        """ Hidden method to compute parameter updates """

        N = Y[0].shape[0]  # this is different from self.N for minibatch
        M = len(Y)
        K = self.dim[1]

        # Masking
        for m in range(M):
            tau[m][mask[m]] = 0.

        weights = [1] * M
        if self.weight_views and M > 1:
            total_w = np.asarray([Y[m].shape[1] for m in range(M)]).sum()
            weights = np.asarray([total_w / (M * Y[m].shape[1]) for m in range(M)])
            weights = weights / weights.sum() * M

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

            # scale Sigma by GP scale parameters alpha
            Qcov[k,:,:] = np.linalg.inv(np.eye(N) * foo[:,k] +  p_cov_inv[k,:,:]) # TODO speed-ups using diagonal structure of first term,
            Qmean[:, k] = gpu_utils.dot(Qcov[k,:,:], bar)

        # Save updated parameters of the Q distribution
        return {'Qmean': Qmean, 'Qcov':Qcov}

    # ELBO calculation per factor - required for grid search for optimal hyperparameter per factor
    def calculateELBO_k(self, k):
        # Collect parameters and expectations of current node
        Qpar, Qexp = self.Q.getParameters(), self.Q.getExpectations()
        Qmean, Qcov = Qpar['mean'], Qpar['cov']
        
        QE  = Qexp['E']

        assert "Sigma" in self.markov_blanket, \
            "Sigma node needs to be in Markov blanket if using Z-GP node"
        Sigma = self.markov_blanket['Sigma'].get_sigma_inv_Terms_k(k)
        p_cov_inv_k = Sigma['inv']
        p_cov_inv_logdet_k = Sigma['inv_logdet']

        # compute term from the exponential in the Gaussian (given by expectation of quadratic form)
        tmp1 = -0.5 * (np.trace(gpu_utils.dot(p_cov_inv_k, Qcov[k,:,:])) +  gpu_utils.dot(QE[:,k].transpose(), gpu_utils.dot(p_cov_inv_k, QE[:,k])))
        # compute term from the precision factor in front of the Gaussian
        tmp2 = 0.5 * p_cov_inv_logdet_k
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
