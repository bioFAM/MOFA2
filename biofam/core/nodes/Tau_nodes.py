from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
import scipy.special as special

from biofam.core.utils import dotd

# Import manually defined functions
from .variational_nodes import Gamma_Unobserved_Variational_Node

from biofam.core.distributions import *


class TauD_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        self.N = self.dim[0]
        self.lbconst = s.sum(self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']))

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

    def updateParameters(self, ix=None, ro=None):
        """
        Public function to update the nodes parameters
        Optional arguments for stochastic updates are:
            - ix: list of indices of the minibatch
            - ro: step size of the natural gradient ascent
        """
        ########################################################################
        # get Expectations which are necessarry for the update
        ########################################################################
        Y = self.markov_blanket["Y"].getExpectation()
        mask = ma.getmask(Y)
        Wtmp = self.markov_blanket["W"].getExpectations()
        Ztmp = self.markov_blanket["Z"].getExpectations()
        N = Y.shape[0]

        W, WW = Wtmp["E"], Wtmp["E2"]
        Z, ZZ = Ztmp["E"], Ztmp["E2"]

        # Collect parameters from the P distributions of this node
        P = self.P.getParameters()
        Pa, Pb = P['a'], P['b']

        ########################################################################
        # subset matrices for stochastic inference
        ########################################################################
        Y = Y.data[ix,:].copy()
        mask = mask[ix,:].copy()
        Z = Z[ix,:].copy()
        ZZ = ZZ[ix,:].copy()

        ########################################################################
        # Masking
        ########################################################################
        Y[mask] = 0.

        ########################################################################
        # compute stochastic "anti-bias" coefficient
        ########################################################################
        coeff = float(N) / float(len(ix))

        ########################################################################
        # compute the update
        ########################################################################
        par_up = self._updateParameters(Y, W, WW, Z, ZZ, Pa, Pb, mask, coeff)

        ########################################################################
        # Do the asignment
        ########################################################################
        if ro is not None: # TODO have a default ro of 1 instead ? whats the overhead cost ?
            par_up['Qa'] = ro * par_up['Qa'] + (1-ro) * self.Q.getParameters()['a']
            par_up['Qb'] = ro * par_up['Qb'] + (1-ro) * self.Q.getParameters()['b']
        self.Q.setParameters(a=par_up['Qa'], b=par_up['Qb'])

    def _updateParameters(self, Y, W, WW, Z, ZZ, Pa, Pb, mask, coeff):
        ZW =  ma.array(Z.dot(W.T), mask=mask)

        # term1 = s.square(Y).sum(axis=0).data # not checked numerically with or without mask
        term1 = s.square(Y.astype("float64")).sum(axis=0)

        # term2 = 2.*(Y*s.dot(Z,W.T)).sum(axis=0).data
        term2 = 2.*(Y*ZW).sum(axis=0) # save to modify

        term3 = ma.array(ZZ.dot(WW.T), mask=mask).sum(axis=0)

        # dotd(WZ, WZ.T) should compute a sum for all k and k' of W_{d,k} Z_{d,k} W_{d,k'} Z_{d,k'} including for k=k', and then summed over n
        # s.dot(s.square(Z),s.square(SW).T) computes the sum over all k of wdk^2 * znk2 summed over all n too
        term4 = dotd(ZW.T, ZW) - ma.array(s.dot(s.square(Z),s.square(W).T), mask=mask).sum(axis=0)

        tmp = term1 - term2 + term3 + term4

        # Perform updates of the Q distribution
        Qa = Pa + coeff * .5 * (Y.shape[0] - mask.sum(axis=0))
        Qb = Pb + coeff * tmp.copy()/2.  # TODO why copy here ?

        return {'Qa': Qa, 'Qb': Qb}


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
        self.precompute()

    def precompute(self):
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

    def updateParameters(self, ix=None, ro=None):
        if ix is not None:
            print('stochastic inference not implemented for taud')
            exit(1)
            
        # Collect expectations from other nodes
        Y = self.markov_blanket["Y"].getExpectation().copy()
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
