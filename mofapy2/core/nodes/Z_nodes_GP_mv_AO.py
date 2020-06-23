from __future__ import division
import numpy as np
import numpy.ma as ma
from copy import deepcopy
from smofapy.core import gp_utils
from smofapy.core import gpu_utils
from smofapy.core.distributions import *
import scipy as s
import sys

# Import manually defined functions
from .variational_nodes import MultivariateGaussian_AO_Unobserved_Variational_Node


# Z_GP_Node_mv
class Z_GP_Node_mv_AO(MultivariateGaussian_AO_Unobserved_Variational_Node):
    """
    Z node with a MultivariateGaussian prior and posterior
    """

    def __init__(self, dim, pmean, pcov, qalpha, qlamb, qE=None):
        super().__init__(dim=dim, pmean=pmean, pcov=pcov, axis_cov=0, qalpha=qalpha, qlamb=qlamb, qE=qE)

        self.mini_batch = None
        self.factors_axis = 1
        self.struct = None
        self.length_scales = None
        self.N = self.dim[0]

    def add_characteristics(self, struct=None, length_scales=None):
        self.struct = np.array(struct)
        self.length_scales = np.array(length_scales)

    def precompute(self, options):
        """ Method to precompute terms to speed up computation """
        gpu_utils.gpu_mode = options['gpu_mode']

    def removeFactors(self, idx, axis=0):
        super().removeFactors(idx, axis)
        self.length_scales = s.delete(self.length_scales, obj=idx)
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
            self.Q.params['K'] = Sigma['cov']
        else:
            self.Q.params['K'] = self.P.params['cov']

        # Get variational parameters of current node
        Q = self.Q.getParameters()
        self.Q.updateExpectations()
        Qmean = self.Q.getExpectation()
        if ix is not None:
            print("Stochastic implemented not implemented in the AO node")
            sys.exit()
        #     self.mini_batch = {}
        #     Qmean = Qmean[ix,:]
        #     Qvar = Qvar[ix,:]

        par_up = self._updateParameters(Y, W, tau, Qmean, self.Q.params['K'], mask)

        # Update parameters
        Q['alpha'] = par_up['alpha']
        Q['lamb'] = par_up['lamb']
        Q['K'] = self.Q.params['K']

        self.Q.setParameters(alpha=Q['alpha'], lamb=Q['lamb'], K = Q['K'])  # NOTE should not be necessary but safer to keep for now

    def _updateParameters(self, Y, W, tau, Qmean, QK, mask):
        """ Hidden method to compute parameter updates """

        N = Y[0].shape[0]  # this is different from self.N for minibatch
        M = len(Y)
        K = self.dim[1]

        # Masking
        for m in range(M):
            tau[m][mask[m]] = 0.

        # Precompute terms to speed up GPU computation
        foo = gpu_utils.array(s.zeros((N, K)))
        precomputed_bar = gpu_utils.array(s.zeros((N, K)))
        for m in range(M):
            tau_gpu = gpu_utils.array(tau[m])
            foo += gpu_utils.dot(tau_gpu, gpu_utils.array(W[m]["E2"]))
            bar_tmp1 = gpu_utils.array(W[m]["E"])
            bar_tmp2 = tau_gpu * gpu_utils.array(Y[m])
            precomputed_bar += gpu_utils.dot(bar_tmp2, bar_tmp1)
        foo = gpu_utils.asnumpy(foo)
        Qlamb = np.sqrt(foo)

        # Calculate variational updates
        # following the approach of Opper & Archambeau (as in https://gpflow.readthedocs.io/en/develop/notebooks/theory/vgp_notes.html)
        Qalpha = gpu_utils.array(s.zeros((N, K)))
        for k in range(K):
            bar = gpu_utils.array(s.zeros((N,)))
            tmp_cp1 = gpu_utils.array(Qmean[:, s.arange(K) != k])
            for m in range(M):
                tmp_cp2 = gpu_utils.array(W[m]["E"][:, s.arange(K) != k].T)
                bar_tmp1 = gpu_utils.array(W[m]["E"][:, k])
                bar_tmp2 = gpu_utils.array(tau[m]) * (-gpu_utils.dot(tmp_cp1, tmp_cp2))

                bar += gpu_utils.dot(bar_tmp2, bar_tmp1)
            bar += precomputed_bar[:, k]
            bar = gpu_utils.asnumpy(bar)
            Qalpha[:, k] = gpu_utils.dot(gpu_utils.dot(np.linalg.inv(QK[k, :, :]), np.linalg.inv(np.diag(Qlamb[:, k]) + QK[k, :, :])), bar)

        # Save updated parameters of the Q distribution
        return {'alpha': Qalpha, 'lamb': Qlamb}

    # Eblo calculation per factor - required for grid search for optimal hyperparameter in sigma per factor
    # following the approach of Opper & Archambeau (as in https://gpflow.readthedocs.io/en/develop/notebooks/theory/vgp_notes.html)
    # KL = 0.5(log |A| + α⊤Kα + tr(A^{−1}) − N)
    # TODO debug
    def calculateELBO_k(self, k):
        # Collect parameters and expectations of current node
        Qpar, Qexp = self.Q.getParameters(), self.Q.getExpectations()
        Qlamb, Qalpha = Qpar['lamb'], Qpar['alpha']

        QE = Qexp['E']

        if 'Sigma' in self.markov_blanket:
            Sigma = self.markov_blanket['Sigma'].getExpectations()
            p_cov = Sigma['cov']
        else:
            p_cov = self.P.params['cov']

        # following Opper & Archambeeu as in https://gpflow.readthedocs.io/en/develop/notebooks/theory/vgp_notes.html and GPflow model SVGP(AO)
        A = gpu_utils.dot(np.diag(Qlamb[:, k]), p_cov[k, :, :])
        A = gpu_utils.dot(A, np.diag(Qlamb[:, k])) + np.eye(self.dim[0])
        L = np.linalg.cholesky(A)
        Li = s.linalg.solve_triangular(L, np.eye(self.N))
        A_logdet = 2.0 * (np.log(np.diag(L))).sum() # https://xcorr.net/2008/06/11/log-determinant-of-positive-definite-matrices-in-matlab/
        trAinv = np.square(Li).sum()
        elb_k = 0.5 * (A_logdet + gpu_utils.dot(gpu_utils.dot(Qalpha[:,k].transpose(), p_cov[k,:,:]), Qalpha[:,k]) + trAinv)

        return elb_k

    # sum up individual EBLO calculations
    def calculateELBO(self):
        elbo = 0
        for k in range(self.dim[1]):
            elbo += self.calculateELBO_k(k)
        elbo -= self.dim[0] * self.dim[1] / 2.

        return elbo

    def sample(self):
        mu = self.P.params['mean']
        if "Sigma" in self.markov_blanket:
            Sigma = self.markov_blanket['Sigma']
            cov = Sigma.sample()
        else:
            cov = self.P.params['cov']

        samp_tmp = [s.random.multivariate_normal(mu[:, i], cov[i, :, :]) for i in range(self.dim[1])]
        self.samp = s.array([tmp - tmp.mean() for tmp in samp_tmp]).transpose()
        return self.samp

