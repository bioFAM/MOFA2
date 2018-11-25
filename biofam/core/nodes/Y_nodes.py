from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
import math

from biofam.core.utils import dotd
from biofam.core import gpu_utils

# Import manually defined functions
from .variational_nodes import Constant_Variational_Node

class Y_Node(Constant_Variational_Node):
    def __init__(self, dim, value):
        Constant_Variational_Node.__init__(self, dim, value)

        # Create a boolean mask of the data to hide missing values
        if type(self.value) != ma.MaskedArray:
            self.mask()

        # Replace NA by zero
        self.value[np.isnan(self.value)] = 0.

    def precompute(self, options=None):
        """ Precompute some terms to speed up the calculations """
        self.N = self.dim[0] - self.mask.sum(axis=0)
        self.D = self.dim[1] - self.mask.sum(axis=1)

        gpu_utils.gpu_mode = options['gpu_mode']

        # Precompute the constant term of the likelihood
        self.likconst = -0.5 * s.sum(self.N) * s.log(2.*s.pi)

    def mask(self):
        """ Generate binary mask """
        self.mask = ma.getmask( ma.masked_invalid(self.value) )

    def getMask(self):
        return self.mask

    def calculateELBO(self):
        """ Calculate evidence lower bound """

        # Collect expectations from nodes
        Y = self.getExpectation()
        Tau = self.markov_blanket["Tau"].getExpectations()
        Wtmp = self.markov_blanket["W"].getExpectations()
        Ztmp = self.markov_blanket["Z"].getExpectations()
        W, WW = Wtmp["E"], Wtmp["E2"]
        Z, ZZ = Ztmp["E"], Ztmp["E2"]

        # Get Mask
        mask = ma.getmask(Y)

        # Compute ELBO
        ZW =  gpu_utils.array(Z).dot(gpu_utils.array(W.T))
        ZW[self.mask] = 0.

        term1 = gpu_utils.square(gpu_utils.array(Y))

        term2 = gpu_utils.array(ZZ).dot(gpu_utils.array(WW.T))
        term2[self.mask] = 0

        term3 = - gpu_utils.dot(gpu_utils.square(gpu_utils.array(Z)),gpu_utils.square(gpu_utils.array(W)).T)
        term3[self.mask] = 0.
        term3 += gpu_utils.square(ZW)

        ZW *= gpu_utils.array(Y)  # WARNING ZW becomes ZWY
        term4 = 2.*ZW

        tmp = 0.5 * (term1 + term2 + term3 - term4)

        Tau["lnE"][self.mask] = 0
        lik = self.likconst + 0.5 * gpu_utils.sum(gpu_utils.array(Tau["lnE"])) -\
              gpu_utils.sum(gpu_utils.array(Tau["E"]) * tmp)

        return lik

    def sample(self, dist='P'):
        # Y does NOT call sample recursively but relies on previous calls
        Z_samp = self.markov_blanket['Z'].samp
        W_samp = self.markov_blanket['W'].samp

        Tau_samp = self.markov_blanket['Tau'].samp
        F = Z_samp.dot(W_samp.transpose())

        # DEPRECATED (tau is expanded inside the node)
        # if Tau_samp.shape != mu.shape:
        #     Tau_samp = s.repeat(Tau_samp.copy()[None,:], self.dim[0], axis=0)
        var = 1./Tau_samp

        if self.markov_blanket['Tau'].__class__.__name__ == "TauN_Node": #TauN
            self.samp = np.array([s.random.normal(F[i, :], math.sqrt(var[i])) for i in range(F.shape[0])])
        else: #TauD
            self.samp = np.array([s.random.normal(F[:, i],math.sqrt(var[i])) for i in range(F.shape[1])]).T

        self.value = self.samp

        return self.samp
