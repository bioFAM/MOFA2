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

        Tau["lnE"][self.mask] = 0

        Y_gpu = gpu_utils.array(Y)
        Z_gpu = gpu_utils.array(Z)
        Wt_gpu = gpu_utils.array(W.T)

        # ZW =  Z_gpu.dot(Wt_gpu)
        # term1 = gpu_utils.square(Y_gpu)
        # term2 = gpu_utils.array(ZZ).dot(gpu_utils.array(WW.T))
        # term3 = - gpu_utils.dot(gpu_utils.square(Z_gpu),gpu_utils.square(Wt_gpu))
        # term3 += gpu_utils.square(ZW)
        # ZW *= Y_gpu  # WARNING ZW becomes ZWY
        # term4 = 2.*ZW
        # tmp = 0.5 * (term1 + term2 + term3 - term4)

        ZW = Z_gpu.dot(Wt_gpu)
        tmp = gpu_utils.square(Y_gpu) \
            + gpu_utils.array(ZZ).dot(gpu_utils.array(WW.T)) \
            - gpu_utils.dot(gpu_utils.square(Z_gpu),gpu_utils.square(Wt_gpu)) + gpu_utils.square(ZW) \
            - 2*ZW*Y_gpu 
        tmp *= 0.5

        tmp[self.mask] = 0.

        # lik = self.likconst + 0.5 * gpu_utils.sum(gpu_utils.array(Tau["lnE"])) - gpu_utils.sum(gpu_utils.array(Tau["E"]) * tmp)
        lik = self.likconst + 0.5 * Tau["lnE"].sum() - gpu_utils.sum(gpu_utils.array(Tau["E"]) * tmp)

        return lik