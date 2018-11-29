from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
import scipy.special as special

from biofam.core.utils import *
from biofam.core import gpu_utils

# Import manually defined functions
from .variational_nodes import Gamma_Unobserved_Variational_Node

from biofam.core.distributions import *


class TauD_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, groups, qE=None):

        self.groups = groups
        self.N = len(self.groups)
        self.n_groups = len(np.unique(groups))

        assert self.n_groups == dim[0], "node dimension does not match number of groups"

        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)

    def precompute(self, options):
        """ Method to precompute some terms to speed up the calculations """

        # GPU mode
        gpu_utils.gpu_mode = options['gpu_mode']

        # Constant ELBO terms
        self.lbconst = s.sum(self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']))

        # compute number of samples per group
        self.n_per_group = np.zeros(self.n_groups)
        for c in range(self.n_groups):
            self.n_per_group[c] = (self.groups == c).sum()

        self.mini_batch = None

    def getExpectations(self, expand=True):
        QExp = self.Q.getExpectations()
        if expand:
            expanded_E = QExp['E'][self.groups, :]
            expanded_lnE = QExp['lnE'][self.groups, :]
            return {'E': expanded_E, 'lnE': expanded_lnE}
        else:
            return QExp

    def getExpectation(self, expand=True):
        QExp = self.getExpectations(expand)
        return QExp['E']

    def define_mini_batch(self, ix):
        """ Method to define minibatch of data for all nodes to use """
        QExp = self.Q.getExpectations()

        # expand only the size of the minibatch
        self.groups_batch = self.groups[ix]
        expanded_E = QExp['E'][self.groups_batch, :]
        # expanded_lnE = s.repeat(QExp['lnE'][None, :], len(ix), axis=0)
        self.mini_batch = expanded_E

    def get_mini_batch(self):
        if self.mini_batch is None:
            return self.getExpectation()
        return self.mini_batch

    def updateParameters(self, ix=None, ro=1.):
        """
        Public function to update the nodes parameters
        Optional arguments for stochastic updates are:
            - ix: list of indices of the minibatch
            - ro: step size of the natural gradient ascent
        """

        # Get expectations from other nodes
        Y = self.markov_blanket["Y"].get_mini_batch()
        mask = self.markov_blanket["Y"].getMask()
        Wtmp = self.markov_blanket["W"].getExpectations()
        Ztmp = self.markov_blanket["Z"].get_mini_batch()
        W, WW = Wtmp["E"], Wtmp["E2"]
        Z, ZZ = Ztmp["E"], Ztmp["E2"]

        # Collect parameters from the P distributions of this node
        P = self.P.getParameters()
        Pa, Pb = P['a'], P['b']

        # subset mini-batch
        if ix is None:
            groups = self.groups
        else:
            groups = self.groups[ix]

        # compute the updates
        Qa, Qb = self._updateParameters(Y, W, WW, Z, ZZ, Pa, Pb, mask, ro, groups)

        self.Q.setParameters(a=Qa, b=Qb)

    def _updateParameters(self, Y, W, WW, Z, ZZ, Pa, Pb, mask, ro, groups):
        """ Hidden method to compute parameter updates """
        Q = self.Q.getParameters()
        Qa, Qb = Q['a'], Q['b']

        # Move matrices to the GPU
        # Y_gpu = gpu_utils.array(Y)
        # Z_gpu = gpu_utils.array(Z)
        # W_gpu = gpu_utils.array(W).T

        # Calculate terms for the update
        ZW =  gpu_utils.array(Z).dot(gpu_utils.array(W.T))
        ZW[mask] = 0.

        term1 = gpu_utils.square(gpu_utils.array(Y))

        term2 = gpu_utils.array(ZZ).dot(gpu_utils.array(WW.T))
        term2[mask] = 0

        term3 = - gpu_utils.dot(gpu_utils.square(gpu_utils.array(Z)),gpu_utils.square(gpu_utils.array(W)).T)
        term3[mask] = 0.
        term3 += gpu_utils.square(ZW)

        ZW *= gpu_utils.array(Y)  # WARNING ZW becomes ZWY
        term4 = 2.*ZW

        tmp = gpu_utils.asnumpy(term1 + term2 + term3 - term4)

        Qa *= (1-ro)
        Qb *= (1-ro)
        for g in range(self.n_groups):
            g_mask = (groups == g)

            n_batch = g_mask.sum()
            if n_batch == 0: continue
            n_total = self.n_per_group[g]
            coeff = n_total/n_batch
            # if ro < 1: import pdb; pdb.set_trace()
            Qb[g,:] += ro * (Pb[g,:] + 0.5*coeff*tmp[g_mask,:].sum(axis=0))

            # TO VERIFY
            mask_tmp = mask[g_mask, :]
            # Qa[g,:] += ro * (Pa[g,:] + (mask.shape[0] - coeff*mask_tmp.sum(axis=0))/2)
            Qa[g,:] += ro * (Pa[g,:] + 0.5*coeff*(mask.shape[0] - mask_tmp.sum(axis=0)))


        return Qa, Qb

    def calculateELBO(self):
        # Collect parameters and expectations from current node
        P, Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.expectations['E'], self.Q.expectations['lnE']

        # Do the calculations
        lb_p = self.lbconst + s.sum((Pa-1.)*QlnE) - s.sum(Pb*QE)
        lb_q = s.sum(Qa*s.log(Qb)) + s.sum((Qa-1.)*QlnE) - s.sum(Qb*QE) - s.sum(special.gammaln(Qa))

        return lb_p - lb_q
