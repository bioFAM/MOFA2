from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
import math

from biofam.core.utils import dotd

# Import manually defined functions
from .variational_nodes import Constant_Variational_Node

class Y_Node(Constant_Variational_Node):
    def __init__(self, dim, value, transpose_noise):
        tau_d = not transpose_noise
        Constant_Variational_Node.__init__(self, dim, value, {'tau_d': tau_d})

        # Create a boolean mask of the data to hide missing values
        if type(self.value) != ma.MaskedArray:
            self.mask()

        # Precompute some terms
        self.precompute()

    def precompute(self):
        # Precompute some terms to speed up the calculations
        self.N = self.dim[0] - ma.getmask(self.value).sum(axis=0)
        self.D = self.dim[1] - ma.getmask(self.value).sum(axis=1)

        # Precompute the constant depending on the noise dimensions
        if self.opts['tau_d']:
            self.likconst = -0.5 * s.sum(self.N) * s.log(2.*s.pi)
        else:
            self.likconst = -0.5 * s.sum(self.D) * s.log(2.*s.pi)

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

        Tau = {k: v[0, :] for (k, v) in self.markov_blanket["Tau"].getExpectations().items()}

        if "SW" in self.markov_blanket:
            Wtmp = self.markov_blanket["SW"].getExpectations()
            Ztmp = self.markov_blanket["Z"].getExpectations()
        else:
            Wtmp = self.markov_blanket["W"].getExpectations()
            Ztmp = self.markov_blanket["SZ"].getExpectations()
        W, WW = Wtmp["E"], Wtmp["E2"]
        Z, ZZ = Ztmp["E"], Ztmp["E2"]

        # Mask matrices
        Y = Y.data
        Y[mask] = 0.

        # term1 = s.square(Y).sum(axis=0).data # not checked numerically with or without mask
        term1 = s.square(Y.astype("float64")).sum(axis=0)

        # term2 = 2.*(Y*s.dot(Z,W.T)).sum(axis=0).data
        term2 = 2. * (Y * s.dot(Z, W.T)).sum(axis=0)  # save to modify

        term3 = ma.array(ZZ.dot(WW.T), mask=mask).sum(axis=0)

        WZ = ma.array(W.dot(Z.T), mask=mask.T)
        term4 = dotd(WZ, WZ.T) - ma.array(s.dot(s.square(Z), s.square(W).T), mask=mask).sum(axis=0)

        tmp = term1 - term2 + term3 + term4
        tmp /= 2.

        if self.opts['tau_d']:
            lik = self.likconst + 0.5 * s.sum(self.N * (Tau["lnE"])) - s.dot(Tau["E"], tmp)
        else:
            lik = self.likconst + 0.5 * s.sum(self.D * (Tau["lnE"])) - s.dot(Tau["E"], tmp)
            
        return lik

    def sample(self, dist='P'):
        # Y does NOT call sample recursively but relies on previous calls
        if "SW" in self.markov_blanket:
            W_samp = self.markov_blanket['SW'].samp
            Z_samp = self.markov_blanket['Z'].samp
        else:
            Z_samp = self.markov_blanket['SZ'].samp
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
