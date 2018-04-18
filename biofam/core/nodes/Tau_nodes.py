
from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
import scipy.special as special

from biofam.core.utils import dotd

# Import manually defined functions
from .variational_nodes import Gamma_Unobserved_Variational_Node


class Tau_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None):
        super(Tau_Node,self).__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.precompute()

    def precompute(self):
        self.D = self.dim[0]
        self.lbconst = s.sum(self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']))

    def updateParameters(self):
        # Collect expectations from other nodes
        Y = self.markov_blanket["Y"].getExpectation().copy()
        tmp = self.markov_blanket["SW"].getExpectations()
        SW,SWW = tmp["E"], tmp["ESWW"]
        Ztmp = self.markov_blanket["Z"].getExpectations()
        Z,ZZ = Ztmp["E"],Ztmp["E2"]
        mask = ma.getmask(Y)

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']

        # Mask matrices
        Y = Y.data
        Y[mask] = 0.

        # Calculate terms for the update
        # term1 = s.square(Y).sum(axis=0).data # not checked numerically with or without mask
        term1 = s.square(Y).sum(axis=0)

        # term2 = 2.*(Y*s.dot(Z,SW.T)).sum(axis=0).data
        term2 = 2.*(Y*s.dot(Z,SW.T)).sum(axis=0) # save to modify

        term3 = ma.array(ZZ.dot(SWW.T), mask=mask).sum(axis=0)

        SWZ = ma.array(SW.dot(Z.T), mask=mask.T)

        # dotd(SWZ, SWZ.T) should compute a sum for all k and k' of W_{d,k} Z_{d,k} W_{d,k'} Z_{d,k'} including for k=k', and then summed over n
        # s.dot(s.square(Z),s.square(SW).T) computes the sum over all k of wdk^2 * znk2 summed over all n too
        term4 = dotd(SWZ, SWZ.T) - ma.array(s.dot(s.square(Z),s.square(SW).T), mask=mask).sum(axis=0)

        tmp = term1 - term2 + term3 + term4

        # Perform updates of the Q distribution
        Qa = Pa + (Y.shape[0] - mask.sum(axis=0))/2.
        Qb = Pb + tmp/2.


        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=Qa, b=Qb)

    def calculateELBO(self):
        # Collect parameters and expectations from current node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.expectations['E'], self.Q.expectations['lnE'] 

        # Do the calculations
        lb_p = self.lbconst + s.sum((Pa-1.)*QlnE) - s.sum(Pb*QE)
        lb_q = s.sum(Qa*s.log(Qb)) + s.sum((Qa-1.)*QlnE) - s.sum(Qb*QE) - s.sum(special.gammaln(Qa))

        return lb_p - lb_q
