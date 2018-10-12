
from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
import scipy.special as special

# Import manually defined functions
from .variational_nodes import Constant_Variational_Node, Beta_Unobserved_Variational_Node

class ThetaW_Node(Beta_Unobserved_Variational_Node):
    """
    This class contain a Theta node associate to factors for which
    we dont have annotations.

    The inference is done per view and factor, so the dimension of the node is the
    number of non-annotated factors

    the updateParameters function needs to know what factors are non-annotated in
    order to choose from the S matrix
    """

    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        # Beta_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)

    def precompute(self, options=None):
        self.factors_axis = 0
        self.Ppar = self.P.getParameters()

    def getExpectations(self, expand=False):
        QExp = self.Q.getExpectations()
        if expand:
            D = self.markov_blanket['W'].D
            expanded_E = s.repeat(QExp['E'][None, :], D, axis=0)
            expanded_lnE = s.repeat(QExp['lnE'][None, :], D, axis=0)
            expanded_lnEInv = s.repeat(QExp['lnEInv'][None, :], D, axis=0)
            return {'E': expanded_E, 'lnE': expanded_lnE, 'lnEInv': expanded_lnEInv}
        else:
            return QExp

    def getExpectation(self, expand=False):
        QExp = self.getExpectations(expand)
        return QExp['E']

    def updateParameters(self, ix=None, factor=None):
        # NOTE Here we use a step of 1 because higher in the hierarchy means useless to decay the step size as W would converge anyway
        self._updateParameters()

    def _updateParameters(self, factors_selection=None):
        # factors_selection (np array or list): indices of factors that are non-annotated

        # Collect expectations from other nodes
        S = self.markov_blanket['W'].getExpectations()["EB"]

        # Precompute terms
        if factors_selection is not None:
            tmp1 = S[:,factors_selection].sum(axis=0)
        else:
            tmp1 = S.sum(axis=0)

        # Perform updates
        Qa = self.Ppar['a'] + tmp1
        Qb = self.Ppar['b'] + S.shape[0] - tmp1

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=Qa, b=Qb)

    def calculateELBO(self):

        # Collect parameters and expectations
        Qpar, Qexp = self.Q.getParameters(), self.Q.getExpectations()
        Pa, Pb, Qa, Qb = self.Ppar['a'], self.Ppar['b'], Qpar['a'], Qpar['b']
        QE, QlnE, QlnEInv = Qexp['E'], Qexp['lnE'], Qexp['lnEInv']

        # minus cross entropy of Q and P
        # lb_p = ma.masked_invalid( (Pa-1.)*QlnE + (Pb-1.)*QlnEInv - special.betaln(Pa,Pb) ).sum()
        lb_p = (Pa-1.)*QlnE + (Pb-1.)*QlnEInv - special.betaln(Pa,Pb)
        lb_p[np.isnan(lb_p)] = 0

        # minus entropy of Q
        # lb_q = ma.masked_invalid( (Qa-1.)*QlnE + (Qb-1.)*QlnEInv - special.betaln(Qa,Qb) ).sum()
        lb_q = (Qa-1.)*QlnE + (Qb-1.)*QlnEInv - special.betaln(Qa,Qb)
        lb_q[np.isnan(lb_q)] = 0

        return lb_p.sum() - lb_q.sum()


class ThetaZ_Node(Beta_Unobserved_Variational_Node):
    """
    Theta node on Z per group.
    Dimensions of the node are number of groups * number of factors
    Implementation is similar to the one of AlphaZ_Node_groups
    """

    def __init__(self, dim, pa, pb, qa, qb, groups, groups_dic, qE=None):

        self.groups = groups
        self.group_names = groups_dic
        self.factors_axis = 1
        self.N = len(self.groups)
        self.n_groups = len(np.unique(groups))

        self.mini_batch = None

        assert self.n_groups == dim[0], "node dimension does not match number of groups"

        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)

    def precompute(self, options=None):
        self.Ppar = self.P.getParameters()
        self.n_per_group = np.zeros(self.n_groups)
        for c in range(self.n_groups):
            self.n_per_group[c] = (self.groups == c).sum()

    def getExpectations(self, expand=False):
        QExp = self.Q.getExpectations()
        if expand:
            expanded_E = QExp['E'][self.groups, :]
            expanded_lnE = QExp['lnE'][self.groups, :]
            expanded_lnEInv = QExp['lnEInv'][self.groups, :]
            return {'E': expanded_E, 'lnE': expanded_lnE, 'lnEInv': expanded_lnEInv}
        else:
            return QExp

    def getExpectation(self, expand=False):
        QExp = self.getExpectations(expand)
        return QExp['E']

    def define_mini_batch(self, ix):
        QExp = self.Q.getExpectations()
        tmp_group = self.groups[ix]
        expanded_expectation = QExp['E'][tmp_group, :]
        expanded_lnE = QExp['lnE'][tmp_group, :]
        expanded_lnEInv = QExp['lnEInv'][tmp_group, :]
        self.mini_batch = {'E': expanded_expectation,
                           'lnE': expanded_lnE,
                           'lnEInv': expanded_lnEInv}

    def get_mini_batch(self):
        if self.mini_batch is None:
            return self.getExpectations(expand=True)
        return self.mini_batch

    def updateParameters(self, ix=None, ro=None, factors_selection=None):
        # Collect expectations from other nodes
        S = self.markov_blanket['Z'].get_mini_batch()["EB"]

        #-----------------------------------------------------------------------
        # subset matrices for stochastic inference
        #-----------------------------------------------------------------------
        if ix is None:
            groups = self.groups
        else:
            groups = self.groups[ix]

        #-----------------------------------------------------------------------
        # 
        #-----------------------------------------------------------------------
        if ro is None:
            ro = 1.

        #-----------------------------------------------------------------------
        # compute the update
        #-----------------------------------------------------------------------
        self._updateParameters(S, groups, ro)

    def _updateParameters(self, S, groups, ro):

        Q = self.Q.getParameters()
        Q['a'] *= (1-ro)
        Q['b'] *= (1-ro)

        # Perform update
        for c in range(self.n_groups):
            mask = (groups == c)

            # coeff for stochastic inference
            n_batch = mask.sum()
            if n_batch == 0: continue
            n_total = self.n_per_group[c]
            coeff = n_total/n_batch

            tmp1 = S[mask, :].sum(axis=0)

            Q['a'][c,:] += ro * (self.Ppar['a'][c,:] + coeff * tmp1)
            Q['b'][c,:] += ro * (self.Ppar['b'][c,:] + coeff * (S[mask, :].shape[0] - tmp1))

    def calculateELBO(self):

        # Collect parameters and expectations
        Qpar, Qexp = self.Q.getParameters(), self.Q.getExpectations()
        Pa, Pb, Qa, Qb = self.Ppar['a'], self.Ppar['b'], Qpar['a'], Qpar['b']
        QE, QlnE, QlnEInv = Qexp['E'], Qexp['lnE'], Qexp['lnEInv']

        # minus cross entropy of Q and P
        # lb_p = ma.masked_invalid( (Pa-1.)*QlnE + (Pb-1.)*QlnEInv - special.betaln(Pa,Pb) ).sum()
        lb_p = (Pa - 1.) * QlnE + (Pb - 1.) * QlnEInv - special.betaln(Pa, Pb)
        lb_p[np.isnan(lb_p)] = 0

        # minus entropy of Q
        # lb_q = ma.masked_invalid( (Qa-1.)*QlnE + (Qb-1.)*QlnEInv - special.betaln(Qa,Qb) ).sum()
        lb_q = (Qa - 1.) * QlnE + (Qb - 1.) * QlnEInv - special.betaln(Qa, Qb)
        lb_q[np.isnan(lb_q)] = 0

        return lb_p.sum() - lb_q.sum()
