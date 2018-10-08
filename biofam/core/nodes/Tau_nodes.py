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
        mask = ma.getmask(Y)
        Y = Y.data

        self.Qa_pre = self.P.getParameters()['a'].copy()

        for g in range(self.n_groups):
            g_mask = (self.groups == g)
            Y_tmp = Y[g_mask, :]
            mask_tmp = mask[g_mask, :]
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
        QExp = self.getExpectations(expand)
        return QExp['E']

    # @profile
    def updateParameters(self):
        # Collect expectations from other nodes
        # TODO sum(axis = 0) to change
        Y = self.markov_blanket["Y"].getExpectation()
        mask = ma.getmask(Y)

        Wtmp = self.markov_blanket["W"].getExpectations()
        Ztmp = self.markov_blanket["Z"].getExpectations()

        W, WW = Wtmp["E"], Wtmp["E2"]
        Z, ZZ = Ztmp["E"], Ztmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P, Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']
        Qb = Q['b']

        # Mask matrices
        Y = Y.data
        Y[mask] = 0.

        # Calculate terms for the update
        ZW =  gpu_utils.array(Z).dot(gpu_utils.array(W.T))
        # ZW =  fast_dot(Z,W.T)
        ZW[mask] = 0.

        term1 = gpu_utils.square(gpu_utils.array(Y)) #.sum(axis=0)

        term2 = gpu_utils.array(ZZ).dot(gpu_utils.array(WW.T))
        # term2 = fast_dot(ZZ, WW.T)
        term2[mask] = 0
        # term2 = term2.sum(axis=0)

        term3 = - gpu_utils.dot(gpu_utils.square(gpu_utils.array(Z)),gpu_utils.square(gpu_utils.array(W)).T)
        term3[mask] = 0.
        # term3 = term3.sum(axis=0)
        term3 += gpu_utils.square(ZW)  #.sum(axis=0)

        ZW *= gpu_utils.array(Y)  # WARNING ZW becomes ZWY
        term4 = 2.*ZW #.sum(axis=0)

        tmp = gpu_utils.asnumpy(term1 + term2 + term3 - term4)

        # Perform updates of the Q distribution
        for g in range(self.n_groups):
            g_mask = (self.groups == g)
            Qb[g,:] = Pb[g,:] + .5 * tmp[g_mask,:].sum(axis=0)

        if s.sum(Qb<0)> 0:
            import pdb; pdb.set_trace()
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

    def sample(self, distrib='P'):
        #instead of overwriting sample, we should maybe change the dimensions of this node !
        P = Gamma(dim=(self.dim[1],1), a=self.P.params["a"][0,:], b=self.P.params["b"][0,:])
        self.samp = P.sample()
        return self.samp


class TauD_Node_bk(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)

    def precompute(self, options):
        self.N = self.dim[0]
        self.lbconst = s.sum(self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']))
        gpu_utils.gpu_mode = options['gpu_mode']

        # update of Qa
        Y = self.markov_blanket["Y"].getExpectation()
        mask = ma.getmask(Y)
        Y = Y.data

        self.Qa_pre = self.P.getParameters()['a'] + (Y.shape[0] - mask.sum(axis=0))/2.

    def getExpectations(self, expand=True):
        QExp = self.Q.getExpectations()
        if expand:
            N = self.markov_blanket['Z'].N
            expanded_E = s.repeat(QExp['E'][None, :], N, axis=0)
            expanded_lnE = s.repeat(QExp['lnE'][None, :], N, axis=0)
            return {'E': expanded_E, 'lnE': expanded_lnE}
        else:
            return QExp

    def getExpectation(self, expand=True):
        QExp = self.getExpectations(expand)
        return QExp['E']

    # @profile
    def updateParameters(self):
        # Collect expectations from other nodes
        Y = self.markov_blanket["Y"].getExpectation()
        mask = ma.getmask(Y)

        Wtmp = self.markov_blanket["W"].getExpectations()
        Ztmp = self.markov_blanket["Z"].getExpectations()

        W, WW = Wtmp["E"], Wtmp["E2"]
        Z, ZZ = Ztmp["E"], Ztmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P, Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']

        # Mask matrices
        Y = Y.data
        Y[mask] = 0.

        # Calculate terms for the update
        ZW =  gpu_utils.array(Z).dot(gpu_utils.array(W.T))
        # ZW =  fast_dot(Z,W.T)
        ZW[mask] = 0.

        import pdb; pdb.set_trace()
        term1 = gpu_utils.square(gpu_utils.array(Y)).sum(axis=0)

        term2 = gpu_utils.array(ZZ).dot(gpu_utils.array(WW.T))
        # term2 = fast_dot(ZZ, WW.T)
        term2[mask] = 0
        term2 = term2.sum(axis=0)

        term3 = gpu_utils.dot(gpu_utils.square(gpu_utils.array(Z)),gpu_utils.square(gpu_utils.array(W)).T)
        term3[mask] = 0.
        term3 = -term3.sum(axis=0)
        term3 += (gpu_utils.square(ZW)).sum(axis=0)

        ZW *= gpu_utils.array(Y)  # WARNING ZW becomes ZWY
        term4 = 2.*(ZW.sum(axis=0))

        tmp = gpu_utils.asnumpy(term1 + term2 + term3 - term4)

        # Perform updates of the Q distribution
        Qb = Pb + tmp/2.

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

    def sample(self, distrib='P'):
        #instead of overwriting sample, we should maybe change the dimensions of this node !
        P = Gamma(dim=(self.dim[1],1), a=self.P.params["a"][0,:], b=self.P.params["b"][0,:])
        self.samp = P.sample()
        return self.samp


class TauN_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)

    def precompute(self, options=None):
        self.D = self.dim[0]
        self.lbconst = s.sum(self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']))

    def getExpectations(self, expand=True):
        QExp = self.Q.getExpectations()
        if expand:
            D = self.markov_blanket['W'].D
            expanded_E = s.repeat(QExp['E'][:, None], D, axis=1)
            expanded_lnE = s.repeat(QExp['lnE'][:, None], D, axis=1)
            return {'E': expanded_E, 'lnE': expanded_lnE}
        else:
            return QExp

    def getExpectation(self, expand=True):
        QExp = self.getExpectations(expand)
        return QExp['E']

    def updateParameters(self):
        # Collect expectations from other nodes
        Y = self.markov_blanket["Y"].getExpectation()
        mask = ma.getmask(Y)

        Wtmp = self.markov_blanket["W"].getExpectations()
        Ztmp = self.markov_blanket["Z"].getExpectations()

        W, WW = Wtmp["E"], Wtmp["E2"]
        Z, ZZ = Ztmp["E"], Ztmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P, Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']

        # Mask matrices
        Y = Y.data
        Y[mask] = 0.

        # Calculate terms for the update

        ZW =  ma.array(Z.dot(W.T), mask=mask)

        term1 = s.square(Y.astype("float64")).sum(axis=1)

        term2 = 2.*(Y*ZW).sum(axis=1) # save to modify

        term3 = ma.array(ZZ.dot(WW.T), mask=mask).sum(axis=1)

        term4 = dotd(ZW, ZW.T) - ma.array(s.dot(s.square(Z), s.square(W).T), mask=mask).sum(axis=1)

        tmp = term1 - term2 + term3 + term4

        # Perform updates of the Q distribution
        Qa = Pa + (Y.shape[1] - mask.sum(axis=1))/2.
        Qb = Pb + tmp.copy()/2.

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=Qa, b=Qb)

    def calculateELBO(self):
        # Collect parameters and expectations from current node
        P, Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.expectations['E'], self.Q.expectations['lnE']

        # Do the calculations
        lb_p = self.lbconst + s.sum((Pa-1.)*QlnE) - s.sum(Pb*QE)
        lb_q = s.sum(Qa*s.log(Qb)) + s.sum((Qa-1.)*QlnE) - s.sum(Qb*QE) - s.sum(special.gammaln(Qa))

        return lb_p - lb_q

    def sample(self, distrib='P'):
        #instead of overwriting sample, we should maybe change the dimensions of this node
        P = Gamma(dim=(self.dim[0], 1), a=self.P.params["a"][:,0], b=self.P.params["b"][:,0])
        self.samp = P.sample()
        return self.samp
