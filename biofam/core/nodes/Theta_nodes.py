
from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
import scipy.special as special

# Import manually defined functions
from .variational_nodes import Constant_Variational_Node, Beta_Unobserved_Variational_Node

class ThetaW_Node_mk(Beta_Unobserved_Variational_Node):
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

    def updateParameters(self, factors_selection=None):
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


class ThetaW_Constant_Node_mk(Constant_Variational_Node):
    """
    Dimensions of Theta_Constant_Node should be (D[m], K)
    """
    def __init__(self, dim, value, N_cells=1):
        super().__init__(dim, value)
        self.N_cells = N_cells

    def precompute(self, options=None):
        self.E = self.value
        # TODO this is wrong with missing values -> need to correct N_cells to account for the cells in which a given gene is missing
        self.lnE = self.N_cells * s.log(self.value)
        self.lnEInv = self.N_cells * s.log(1.-self.value)

    def getExpectations(self):
        return { 'E':self.E, 'lnE':self.lnE, 'lnEInv':self.lnEInv }

    def removeFactors(self, idx, axis=1):
        # Ideally we want this node to use the removeFactors defined in Node()
        # but the problem is that we also need to update the "expectations", so i need
        # to call precompute()
        self.value = s.delete(self.value, idx, axis)
        self.precompute()
        self.updateDim(axis=axis, new_dim=self.dim[axis]-len(idx))

class ThetaZ_Node_k(Beta_Unobserved_Variational_Node):
    """
    This class contain a Theta node associate to factors for which
    we dont have annotations.
    The inference is done per factor, so the dimension of the node is the
    number of non-annotated factors
    the updateParameters function needs to know what factors are non-annotated in
    order to choose from the T matrix
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
            N = self.markov_blanket['Z'].N
            expanded_E = s.repeat(QExp['E'][None, :], N, axis=0)
            expanded_lnE = s.repeat(QExp['lnE'][None, :], N, axis=0)
            expanded_lnEInv = s.repeat(QExp['lnEInv'][None, :], N, axis=0)
            return {'E': expanded_E, 'lnE': expanded_lnE, 'lnEInv': expanded_lnEInv}
        else:
            return QExp

    def getExpectation(self, expand=False):
        QExp = self.getExpectations(expand)
        return QExp['E']

    def updateParameters(self, factors_selection=None):
        # factors_selection (np array or list): indices of factors that are non-annotated

        # Collect expectations from other nodes
        S = self.markov_blanket['Z'].getExpectations()["EB"]

        # Precompute terms
        if factors_selection is not None:
            tmp1 = S[:, factors_selection].sum(axis=0)
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
        lb_p = (Pa - 1.) * QlnE + (Pb - 1.) * QlnEInv - special.betaln(Pa, Pb)
        lb_p[np.isnan(lb_p)] = 0

        # minus entropy of Q
        # lb_q = ma.masked_invalid( (Qa-1.)*QlnE + (Qb-1.)*QlnEInv - special.betaln(Qa,Qb) ).sum()
        lb_q = (Qa - 1.) * QlnE + (Qb - 1.) * QlnEInv - special.betaln(Qa, Qb)
        lb_q[np.isnan(lb_q)] = 0

        return lb_p.sum() - lb_q.sum()

class ThetaZ_Node_groups(Beta_Unobserved_Variational_Node):
    """
    Theta node on Z per group.
    Dimensions of the node are number of groups * number of factors
    Implementation is similar to the one of AlphaZ_Node_groups
    """

    def __init__(self, dim, pa, pb, qa, qb, groups, groups_dic, qE=None):

        self.groups = groups
        self.group_names = groups_dic
        self.N = len(self.groups)
        self.n_groups = len(np.unique(groups))

        assert self.n_groups == dim[0], "node dimension does not match number of groups"

        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)

    def precompute(self, options=None):
        self.factors_axis = 1
        self.Ppar = self.P.getParameters()

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

    def updateParameters(self, factors_selection=None):
        # print(self.getExpectation())
        # factors_selection (np array or list): indices of factors that are non-annotated
        # collect local parameters
        Q = self.Q.getParameters()
        Qa, Qb = Q['a'], Q['b']

        # Collect expectations from other nodes
        S = self.markov_blanket['Z'].getExpectations()["EB"]

        # Precompute terms
        if factors_selection is not None:
            tmpS = S[:, factors_selection]
        else:
            tmpS = S

        # Perform update
        for c in range(self.n_groups):
            mask = (self.groups == c)
            tmp1 = tmpS[mask, :].sum(axis=0)
            Qa[c,:] = self.Ppar['a'][c,:] + tmp1
            Qb[c,:] = self.Ppar['b'][c,:] + S[mask, :].shape[0] - tmp1

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=Qa, b=Qb)

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


class ThetaZ_Constant_Node_k(Constant_Variational_Node):
    """
    Dimensions of Theta_Constant_Node should be (N, K)
    """

    def __init__(self, dim, value, N_cells=1):
        super().__init__(dim, value)
        self.N_cells = N_cells

    def precompute(self, options=None):
        self.E = self.value
        # TODO this is wrong with missing values -> need to correct N_cells to account for the cells in which a given gene is missing
        self.lnE = self.N_cells * s.log(self.value)
        self.lnEInv = self.N_cells * s.log(1. - self.value)

    def getExpectations(self):
        return {'E': self.E, 'lnE': self.lnE, 'lnEInv': self.lnEInv}

    def removeFactors(self, idx, axis=1):
        # Ideally we want this node to use the removeFactors defined in Node()
        # but the problem is that we also need to update the "expectations", so i need
        # to call precompute()
        self.value = s.delete(self.value, idx, axis)
        self.precompute()
        self.updateDim(axis=axis, new_dim=self.dim[axis] - len(idx))
