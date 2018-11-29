
from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
import scipy.special as special

# Import manually defined functions
from .variational_nodes import Gamma_Unobserved_Variational_Node

class AlphaW_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None, qlnE=None):
        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE, qlnE=qlnE)

    def precompute(self, options=None):
        self.factors_axis = 0

    def getExpectations(self, expand=False):
        QExp = self.Q.getExpectations()
        QExp['E'] = QExp['E']
        QExp['lnE'] = QExp['lnE']
        if expand:
            D = self.markov_blanket['W'].dim[0]
            expanded_E = s.repeat(QExp['E'][None, :], D, axis=0)
            expanded_lnE = s.repeat(QExp['lnE'][None, :], D, axis=0)
            return {'E': expanded_E, 'lnE': expanded_lnE}
        else:
            return QExp

    def getExpectation(self, expand=False):
        QExp = self.getExpectations(expand)
        return QExp['E']

    def updateParameters(self, ix=None, ro=None):
        # NOTE Here we use a step of 1 because higher in the hierarchy means useless to decay the step size as W would converge anyway
        self._updateParameters()

    def _updateParameters(self):
        # Collect expectations from other nodes
        tmp = self.markov_blanket["W"].getExpectations()
        E  = tmp["E"]
        if "ENN" in tmp:
            EWW = tmp["ENN"]
        else:
            EWW = tmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P = self.P.getParameters()
        Pa, Pb = P['a'], P['b']

        # Perform updates
        Qa = Pa + 0.5*E.shape[0]
        Qb = Pb + 0.5*EWW.sum(axis=0)

        self.Q.setParameters(a=Qa, b=Qb)

    def calculateELBO(self):
        # Collect parameters and expectations
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.getExpectations()['E'], self.Q.getExpectations()['lnE']

        # Do the calculations
        lb_p = (Pa*s.log(Pb)).sum() - special.gammaln(Pa).sum() + ((Pa-1.)*QlnE).sum() - (Pb*QE).sum()
        lb_q = (Qa*s.log(Qb)).sum() - special.gammaln(Qa).sum() + ((Qa-1.)*QlnE).sum() - (Qb*QE).sum()

        return lb_p - lb_q

class AlphaZ_Node(Gamma_Unobserved_Variational_Node):
    """ """

    def __init__(self, dim, pa, pb, qa, qb, groups, qE=None, qlnE=None):
        self.groups = groups
        self.factors_axis = 1
        self.N = len(self.groups)
        self.n_groups = len(np.unique(groups))

        self.mini_batch = None

        assert self.n_groups == dim[0], "node dimension does not match number of groups"

        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE, qlnE=qlnE)

    def precompute(self, options):
        self.n_per_group = np.zeros(self.n_groups)
        for c in range(self.n_groups):
            self.n_per_group[c] = (self.groups == c).sum()

    def getExpectations(self, expand=False):
        QExp = self.Q.getExpectations()
        if expand:
            expanded_expectation = QExp['E'][self.groups, :]
            expanded_lnE = QExp['lnE'][self.groups, :]
            return {'E': expanded_expectation, 'lnE': expanded_lnE}
        else:
            return {'E': QExp['E'], 'lnE': QExp['lnE']}

    def getExpectation(self, expand=False):
        QExp = self.getExpectations(expand)
        return QExp['E']

    def define_mini_batch(self, ix):
        QExp = self.Q.getExpectations()
        tmp_group = self.groups[ix]
        expanded_expectation = QExp['E'][tmp_group, :]
        # expanded_lnE = QExp['lnE'][tmp_group, :]
        self.mini_batch = expanded_expectation

    def get_mini_batch(self):
        if self.mini_batch is None:
            return self.getExpectation(expand=True)
        return self.mini_batch

    def updateParameters(self, ix=None, ro=1.):
        """
        Public method to update the nodes parameters
        Optional arguments for stochastic updates are:
            - ix: list of indices of the minibatch
            - ro: step size of the natural gradient ascent
        """
        tmp = self.markov_blanket["Z"].get_mini_batch()
        if 'ENN' in tmp:
            EZZ = tmp["ENN"]
        else:
            EZZ = tmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P = self.P.getParameters()
        Pa, Pb = P['a'], P['b']

        #-----------------------------------------------------------------------
        # subset matrices for stochastic inference
        #-----------------------------------------------------------------------
        if ix is None:
            groups = self.groups
        else:
            groups = self.groups[ix]

        #-----------------------------------------------------------------------
        # compute the update
        #-----------------------------------------------------------------------
        self._updateParameters(Pa, Pb, EZZ, groups, ro)


    def _updateParameters(self, Pa, Pb, EZZ, groups, ro):
        """ Hidden method to compute parameter updates """

        Q = self.Q.getParameters()
        Q['a'] *= (1-ro)
        Q['b'] *= (1-ro)

        for c in range(self.n_groups):
            mask = (groups == c)

            # coeff for stochastic inference
            n_batch = mask.sum()
            if n_batch == 0: continue
            n_total = self.n_per_group[c]
            coeff = n_total/n_batch

            Q['a'][c,:] += ro * (Pa[c,:] + 0.5 * n_total)  # TODO should be precomputed
            Q['b'][c,:] += ro * (Pb[c,:] + 0.5 * coeff * EZZ[mask, :].sum(axis=0))

    def calculateELBO(self):
        # Collect parameters and expectations
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.getExpectations()['E'], self.Q.getExpectations()['lnE']

        # Do the calculations
        lb_p = (Pa*s.log(Pb)).sum() - special.gammaln(Pa).sum() + ((Pa-1.)*QlnE).sum() - (Pb*QE).sum()
        lb_q = (Qa*s.log(Qb)).sum() - special.gammaln(Qa).sum() + ((Qa-1.)*QlnE).sum() - (Qb*QE).sum()

        return lb_p - lb_q
