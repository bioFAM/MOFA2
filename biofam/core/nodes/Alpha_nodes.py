
from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
import scipy.special as special

# Import manually defined functions
from .variational_nodes import Gamma_Unobserved_Variational_Node

# TODO add sample functions everywhere
# TODO calculateELBO is the same and could be moved to the parent node ?
class AlphaW_Node_mk(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None, qlnE=None):
        # Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE, qlnE=qlnE)
        self.precompute()

    def precompute(self):
        # self.lbconst = self.K * ( self.P.a*s.log(self.P.b) - special.gammaln(self.P.a) )
        # self.lbconst = s.sum( self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']) )
        self.factors_axis = 0  # TODO check that !

    def updateParameters(self):

        # Collect expectations from other nodes
        if "SW" in self.markov_blanket:
            tmp = self.markov_blanket["SW"].getExpectations()
            ES  = tmp["EB"]
            # TODO what is ENN and is it really what we want and not E2 ? (eternal question)
            EWW = tmp["ENN"]
        else:
            tmp = self.markov_blanket["W"].getExpectations()
            ES  = tmp["E"]
            EWW = tmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']

        # Perform updates
        Qa = Pa + 0.5*ES.shape[0]
        Qb = Pb + 0.5*EWW.sum(axis=0)

        # Save updated parameters of the Q distribution
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


class AlphaZ_Node_k(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None, qlnE=None):
        # Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE, qlnE=qlnE)
        self.precompute()

    def precompute(self):
        # self.lbconst = self.K * ( self.P.a*s.log(self.P.b) - special.gammaln(self.P.a) )
        # self.lbconst = s.sum( self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']) )
        self.factors_axis = 0   # TODO check that !!!!!

    def updateParameters(self):

        # TODO: muZ node not accounted for here
        # Collect expectations from other node
        if "SZ" in self.markov_blanket:
            tmp = self.markov_blanket["SZ"].getExpectations()
            EZZ = tmp["ENN"]
        else:
            tmp = self.markov_blanket["Z"].getExpectations()
            EZZ = tmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']

        # Perform update
        Qa = Pa + 0.5*EZZ.shape[0]
        Qb = Pb + 0.5*EZZ.sum(axis=0)

        # Save updated parameters of the Q distribution
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


class AlphaZ_Node_groups(Gamma_Unobserved_Variational_Node):
    """ """

    def __init__(self, dim, pa, pb, qa, qb, groups, groups_dic=None, qE=None, qlnE=None):
        self.groups = groups
        self.factors_axis = 1
        self.N = len(self.groups)
        self.n_groups = len(np.unique(groups))

        assert self.n_groups == dim[0], "node dimension does not match number of groups"

        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE, qlnE=qlnE)

    def getExpectations(self):
        # reshape the values to N_samples * N_factors and return
        QExp = self.Q.getExpectations()
        expanded_expectation = QExp['E'][self.groups, :]
        expanded_lnE = QExp['lnE'][self.groups, :]
        # do we need to expand the variance as well -> not used I think
        return {'E': expanded_expectation, 'lnE': expanded_lnE}

    def getExpectation(self):
        # reshape the values to N_samples * N_factors and return
        QExp = self.Q.getExpectations()
        expanded_expectation = QExp['E'][self.groups, :]
        return expanded_expectation

    def updateParameters(self):
        # TODO: add an if MuZ is in markov blanket ?
        if "SZ" in self.markov_blanket:
            tmp = self.markov_blanket["SZ"].getExpectations()
            EZZ = tmp["ENN"]
        else:
            tmp = self.markov_blanket["Z"].getExpectations()
            EZZ = tmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']
        Qa, Qb = Q['a'], Q['b']

        # Perform update
        for c in range(self.n_groups):
            mask = (self.groups == c)

            Qa[c,:] = Pa[c,:] + 0.5*EZZ[mask, :].shape[0]
            Qb[c,:] = Pb[c,:] + 0.5*EZZ[mask, :].sum(axis=0)

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
