from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s

# Import manually defined functions
from .variational_nodes import Constant_Variational_Node

class Y_Node(Constant_Variational_Node):
    def __init__(self, dim, value):
        Constant_Variational_Node.__init__(self, dim, value)

        # Create a boolean mask of the data to hide missing values
        if type(self.value) != ma.MaskedArray:
            self.mask()

        # Precompute some terms
        self.precompute()

    def precompute(self):
        # Precompute some terms to speed up the calculations
        self.N = self.dim[0] - ma.getmask(self.value).sum(axis=0)
        self.D = self.dim[1]
        self.likconst = -0.5*s.sum(self.N * self.D)*s.log(2.*s.pi)

    def mask(self):
        # Mask the observations if they have missing values
        self.value = ma.masked_invalid(self.value)

    def getMask(self):
        return ma.getmask(self.value)

    def calculateELBO(self):
        # Calculate evidence lower bound
        # We use the trick that the update of Tau already contains the Gaussian likelihod.
        # However, it is important that the lower bound is calculated after the update of Tau is performed
        if "SW" in self.markov_blanket:
            tauQ_param = {k:v[0,:] for (k, v) in self.markov_blanket["Tau"].getParameters("Q").items()}
            tauP_param = {k:v[0,:] for (k, v) in self.markov_blanket["Tau"].getParameters("P").items()}
            tau_exp = {k:v[0,:] for (k, v) in self.markov_blanket["Tau"].getExpectations().items()}
            lik = self.likconst + 0.5*s.sum(self.N*(tau_exp["lnE"])) - s.dot(tau_exp["E"], tauQ_param["b"]-tauP_param["b"])
        else:
            tauQ_param = {k:v[:,0] for (k, v) in self.markov_blanket["Tau"].getParameters("Q").items()}
            tauP_param = {k:v[:,0] for (k, v) in self.markov_blanket["Tau"].getParameters("P").items()}
            tau_exp = {k:v[:,0] for (k, v) in self.markov_blanket["Tau"].getExpectations().items()}
            lik = self.likconst + 0.5*s.sum(self.D*(tau_exp["lnE"])) - s.dot(tau_exp["E"], tauQ_param["b"]-tauP_param["b"])
        return lik

    def sample(self, dist='P'):
        # Y does NOT call sample recursively but relies on previous calls
        Z_samp = self.markov_blanket['Z'].samp
        W_samp = self.markov_blanket['SW'].samp
        Tau_samp = self.markov_blanket['Tau'].samp

        mu = Z_samp.dot(W_samp.transpose())
        # DEPRECATED (tau is expanded inside the node)
        # if Tau_samp.shape != mu.shape:
        #     Tau_samp = s.repeat(Tau_samp.copy()[None,:], self.dim[0], axis=0)
        var = 1./Tau_samp

        self.samp = s.random.normal(mu, var)
        return self.samp
