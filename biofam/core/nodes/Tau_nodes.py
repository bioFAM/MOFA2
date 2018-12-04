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
    def __init__(self, dim, pa, pb, qa, qb, groups, groups_dic, qE=None):

        self.groups = groups
        self.group_names = groups_dic
        self.N = len(self.groups)
        self.n_groups = len(np.unique(groups))

        assert self.n_groups == dim[0], "node dimension does not match number of groups"

        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)

    def precompute(self, options):
        # self.N = self.dim[0]
        self.lbconst = s.sum(self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']))
        gpu_utils.gpu_mode = options['gpu_mode']

        # update of Qa
        Y = self.markov_blanket["Y"].getExpectation()
        self.mask = self.markov_blanket["Y"].getMask()
        # mask = ma.getmask(Y)
        # Y = Y.data

        self.Qa_pre = self.P.getParameters()['a'].copy()

        for g in range(self.n_groups):
            g_mask = (self.groups == g)
            Y_tmp = Y[g_mask, :]
            mask_tmp = self.mask[g_mask, :]
            self.Qa_pre[g,:] += (Y_tmp.shape[0] - mask_tmp.sum(axis=0))/2.

    def getExpectations(self, expand=True):
        QExp = self.Q.getExpectations()
        if expand:
            expanded_E = QExp['E'][self.groups, :]
            expanded_lnE = QExp['lnE'][self.groups, :]
            return {'E': expanded_E, 'lnE': expanded_lnE}
        else:
            return QExp

    def getExpectation(self, expand=True):
        QExp = self.Q.getExpectations()
        if expand:
            return QExp['E'][self.groups, :]
        else:
            return Qexp

    def updateParameters(self):
        # Collect expectations from other nodes

        Y = self.markov_blanket["Y"].getExpectation()

        Wtmp = self.markov_blanket["W"].getExpectations()
        Ztmp = self.markov_blanket["Z"].getExpectations()

        W, WW = Wtmp["E"], Wtmp["E2"]
        Z, ZZ = Ztmp["E"], Ztmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P, Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']
        Qb = Q['b']

        # Copy matrices to GPU
        Y_gpu = gpu_utils.array(Y)
        Z_gpu = gpu_utils.array(Z)
        W_gpu = gpu_utils.array(W).T

        ZW = Z_gpu.dot(W_gpu)
        tmp = gpu_utils.square(Y_gpu) \
            + gpu_utils.array(ZZ).dot(gpu_utils.array(WW.T)) \
            - gpu_utils.dot(gpu_utils.square(Z_gpu),gpu_utils.square(W_gpu)) + gpu_utils.square(ZW) \
            - 2*ZW*Y_gpu 
        tmp[self.mask] = 0.

        # Perform updates of the Q distribution
        for g in range(self.n_groups):
            g_mask = (self.groups == g)
            Qb[g,:] = Pb[g,:] + tmp[g_mask,:].sum(axis=0) / 2.

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=self.Qa_pre, b=Qb)

    def calculateELBO(self):
        # Collect parameters and expectations from current node
        P, Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.expectations['E'], self.Q.expectations['lnE']

        # Do the calculations
        lb_p = self.lbconst + s.sum((Pa-1.)*QlnE) - s.sum(Pb*QE)
        lb_q = s.sum(Qa*s.log(Qb)) + s.sum((Qa-1.)*QlnE) - s.sum(Qb*QE) - s.sum(special.gammaln(Qa))

        return lb_p - lb_q
