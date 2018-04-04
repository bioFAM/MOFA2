
from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
import scipy.special as special

# Import manually defined functions
from .variational_nodes import Gamma_Unobserved_Variational_Node

class AlphaW_Node_mk(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None, qlnE=None):
        # Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        super(AlphaW_Node_mk,self).__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE, qlnE=qlnE)
        self.precompute()

    def precompute(self):
        # self.lbconst = self.K * ( self.P.a*s.log(self.P.b) - special.gammaln(self.P.a) )
        # self.lbconst = s.sum( self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']) )
        self.factors_axis = 0

    def updateParameters(self):

        # Collect expectations from other nodes
        if "SW" in self.markov_blanket:
            tmp = self.markov_blanket["SW"].getExpectations()
            EWW = tmp["ENN"]
        else:
            tmp = self.markov_blanket["W"].getExpectations()
            EWW = tmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']

        # Perform updates
        Qa = Pa + 0.5*EWW.shape[0]
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
        super(AlphaZ_Node_k,self).__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE, qlnE=qlnE)
        self.precompute()

    def precompute(self):
        # self.lbconst = self.K * ( self.P.a*s.log(self.P.b) - special.gammaln(self.P.a) )
        # self.lbconst = s.sum( self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']) )
        self.factors_axis = 0

    def updateParameters(self):

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
