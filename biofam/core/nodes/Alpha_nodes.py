
from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
import scipy.special as special

# Import manually defined functions
from .variational_nodes import Gamma_Unobserved_Variational_Node

# TODO add sample functions everywhere 
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

    # TODO: for now it is assumed that clusters are a list of digits from 0 to n_clust. we should handle other types
    # and convert it (strings etc)
    def __init__(self, dim, pa, pb, qa, qb, clusters, cluster_dic=None, qE=None, qE2=None):
        self.clusters = clusters
        self.factors_axis = 1
        self.N = len(self.clusters)
        self.n_clusters = len(np.unique(clusters))

        assert self.n_clusters == dim[0], "node dimension does not match number of clusters"

        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE, qlnE=qlnE)

    def getExpectations(self):
        # reshape the values to N_samples * N_factors and return
        QExp = self.Q.getExpectations()
        expanded_expectation = QExp['E'][self.clusters, :]
        expanded_lnE = QExp['lnE'][self.clusters, :]
        # do we need to expand the variance as well -> not used I think
        return {'E': expanded_expectation, 'lnE': expanded_lnE}

    def updateParameters(self):
        # TODO: add an if MuZ is in markov blanket ?
        if "SZ" in self.markov_blanket:
            tmp = self.markov_blanket["SZ"].getExpectations()
            EZZ = tmp["ENN"]
        else:
            tmp = self.markov_blanket["Z"].getExpectations()
            EZZ = tmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters()
        Pa, Pb = P['a'], P['b']
        Qa, Qb = Q['a'], Q['b']

        # Perform update
        for c in range(self.n_clusters):
            mask = (self.clusters == c)

            Qa[c,:] = Pa[c,:] + 0.5*EZZ[mask, :].shape[0]
            Qb[c,:] = Pb[c,:] + 0.5*EZZ[mask, :].sum(axis=0)

        # TODO not even necessary as we already have a reference above ?
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
