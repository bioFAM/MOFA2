from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
import math

from biofam.core.utils import dotd

# Import manually defined functions
from .variational_nodes import Constant_Variational_Node

class Y_Node(Constant_Variational_Node):
    def __init__(self, dim, value):
        Constant_Variational_Node.__init__(self, dim, value)

        # Create a boolean mask of the data to hide missing values
        if type(self.value) != ma.MaskedArray:
            self.mask()

    def precompute(self, options=None):
        # Precompute some terms to speed up the calculations
        self.N = self.dim[0] - ma.getmask(self.value).sum(axis=0)
        self.D = self.dim[1] - ma.getmask(self.value).sum(axis=1)

        # Precompute the constant term of the likelihood
        self.likconst = -0.5 * s.sum(self.N) * s.log(2.*s.pi)

    def mask(self):
        # Mask the observations if they have missing values
        self.value = ma.masked_invalid(self.value)

    def getMask(self):
        return ma.getmask(self.value)

    def calculateELBO(self):
        # Calculate evidence lower bound
        # Collect expectations from nodes

        Y = self.getExpectation().copy()
        mask = ma.getmask(Y)

        Tau = self.markov_blanket["Tau"].getExpectations()

        Wtmp = self.markov_blanket["W"].getExpectations()
        Ztmp = self.markov_blanket["Z"].getExpectations()

        W, WW = Wtmp["E"], Wtmp["E2"]
        Z, ZZ = Ztmp["E"], Ztmp["E2"]

        # Mask matrices
        Y = Y.data
        Y[mask] = 0.

        ZW = ma.array(s.dot(Z, W.T),mask=mask)

        # term1 = s.square(Y).sum(axis=0).data # not checked numerically with or without mask
        term1 = s.square(Y.astype("float64"))

        # term2 = 2.*(Y*s.dot(Z,W.T)).sum(axis=0).data
        term2 = 2. * (Y * ZW)  # save to modify

        term3 = ma.array(ZZ.dot(WW.T), mask=mask)

        term4 = s.square(ZW) - ma.array(s.dot(s.square(Z),s.square(W).T),mask=mask)

        tmp = 0.5 * (term1 - term2 + term3 + term4)

        Tau["lnE"][mask] = 0
        lik = self.likconst + 0.5 * s.sum(Tau["lnE"]) - s.sum(s.multiply(Tau["E"], tmp))
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
