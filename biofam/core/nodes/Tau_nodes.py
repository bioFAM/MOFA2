from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
import scipy.special as special

from biofam.core.utils import *

# Import manually defined functions
from .variational_nodes import Gamma_Unobserved_Variational_Node

from biofam.core.distributions import *


class TauD_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        # self.precompute()

    def precompute(self):
        self.N = self.dim[0]
        self.lbconst = s.sum(self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']))

        # update of Qa
        Y = self.markov_blanket["Y"].getExpectation()
        mask = ma.getmask(Y)
        Y = Y.data

        self.mini_batch = None

        # TODO is it ok to do that with stochastic ? here we kind of assume that missing values are evenly distributed
        self.Qa_pre = self.P.getParameters()['a'] + (Y.shape[0] - mask.sum(axis=0))/2.

    def getExpectations(self, expand=True):
        QExp = self.Q.getExpectations()
        if expand:
            N = self.markov_blanket['Z'].dim[0]
            expanded_E = s.repeat(QExp['E'][None, :], N, axis=0)
            expanded_lnE = s.repeat(QExp['lnE'][None, :], N, axis=0)
            return {'E': expanded_E, 'lnE': expanded_lnE}
        else:
            return QExp

    def getExpectation(self, expand=True):
        QExp = self.getExpectations(expand)
        return QExp['E']

    def define_mini_batch(self, ix):
        # define minibatch of data for all nodes to use
        QExp = self.Q.getExpectations()

        # expand only the size of the minibatch
        expanded_E = s.repeat(QExp['E'][None, :], len(ix), axis=0)
        # expanded_lnE = s.repeat(QExp['lnE'][None, :], len(ix), axis=0)
        self.mini_batch = expanded_E

    def get_mini_batch(self):
        if self.mini_batch is None:
            return self.getExpectation()
        return self.mini_batch

    def updateParameters(self, ix=None, ro=None):
        """
        Public function to update the nodes parameters
        Optional arguments for stochastic updates are:
            - ix: list of indices of the minibatch
            - ro: step size of the natural gradient ascent
        """
        ########################################################################
        # get full Expectations or minibatches
        ########################################################################
        Y = self.markov_blanket["Y"].get_mini_batch()
        Wtmp = self.markov_blanket["W"].getExpectations()
        Ztmp = self.markov_blanket["Z"].get_mini_batch()
        N = self.markov_blanket["Y"].dim[0]

        W, WW = Wtmp["E"], Wtmp["E2"]
        Z, ZZ = Ztmp["E"], Ztmp["E2"]

        # Collect parameters from the P distributions of this node
        P = self.P.getParameters()
        Pa, Pb = P['a'], P['b']

        ########################################################################
        # Masking
        ########################################################################
        mask = ma.getmask(Y)
        Y = Y.data
        Y[mask] = 0.

        ########################################################################
        # compute stochastic "anti-bias" coefficient
        ########################################################################
        coeff = float(N) / float(Y.shape[0])
        # TODO fix that in BayesNet instead
        if ro is None:
            ro =1.

        ########################################################################
        # compute the update
        ########################################################################
        par_up = self._updateParameters(Y, W, WW, Z, ZZ, Pa, Pb, mask, coeff, ro)

        ########################################################################
        # Do the asignment
        # ########################################################################
        # if ro is not None: # TODO have a default ro of 1 instead ? whats the overhead cost ?
        #     par_up['Qa'] = ro * par_up['Qa'] + (1-ro) * self.Q.getParameters()['a']
        #     par_up['Qb'] = ro * par_up['Qb'] + (1-ro) * self.Q.getParameters()['b']
        # self.Q.setParameters(a=par_up['Qa'], b=par_up['Qb'])

    def _updateParameters(self, Y, W, WW, Z, ZZ, Pa, Pb, mask, coeff, ro):

        # Calculate terms for the update
        ZW =  Z.dot(W.T)
        ZW[mask] = 0.

        term1 = s.square(Y).sum(axis=0)

        term2 = ZZ.dot(WW.T)
        term2[mask] = 0
        term2 = term2.sum(axis=0)

        term3 = np.dot(np.square(Z),np.square(W).T)
        term3[mask] = 0.
        term3 = -term3.sum(axis=0)
        term3 += (np.square(ZW)).sum(axis=0)

        ZW *= Y  # WARNING ZW becomes ZWY
        term4 = 2.*(ZW.sum(axis=0))

        tmp = term1 + term2 + term3 - term4
        tmp *= coeff

        Qb = self.Q.getParameters()['b']
        Qb *= (1-ro)
        Qb += ro * (Pb + tmp/2.)


        Qa = self.Q.getParameters()['a']
        Qa *= (1 - ro)
        Qa += ro * self.Qa_pre

        # Qb = Pb + tmp/2.
        #
        # # return updated parameters of the Q distribution
        # return {'Qa': self.Q.params['a'], 'Qb': Qb}

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
