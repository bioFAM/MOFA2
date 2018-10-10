
from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
import scipy.special as special

# Import manually defined functions
from .variational_nodes import Gamma_Unobserved_Variational_Node

# TODO add sample functions everywhere
# TODO calculateELBO is the same and could be moved to the parent node ?
# TODO actually all the nodes are exactly the same apart from labeling. Could be one single node
class AlphaW_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None, qlnE=None):
        # Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE, qlnE=qlnE)

    def precompute(self, options=None):
        # self.lbconst = self.K * ( self.P.a*s.log(self.P.b) - special.gammaln(self.P.a) )
        # self.lbconst = s.sum( self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']) )
        self.factors_axis = 0

    def getExpectations(self, expand=False):
        QExp = self.Q.getExpectations()
        QExp['E'] = QExp['E']
        QExp['lnE'] = QExp['lnE']
        if expand:
            D = self.markov_blanket['W'].D
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
        # TODO should we change that for stochastic ? considering that the forgetting rate already acts on W which is
        # TODO only node in the markov blanket
        # Collect expectations from other nodes
        tmp = self.markov_blanket["W"].getExpectations()
        E  = tmp["E"]
        if "ENN" in tmp:
            EWW = tmp["ENN"]
        else:
            EWW = tmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']
        Qa, Qb = Q['a'], Q['b']

        # Perform updates
        Qa = Pa + 0.5*E.shape[0]
        Qb = Pb + 0.5*EWW.sum(axis=0)

        # return {'Qa': Qa, 'Qb': Qb}


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

    def __init__(self, dim, pa, pb, qa, qb, groups, groups_dic, qE=None, qlnE=None):
        self.groups = groups
        self.group_names = groups_dic
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
            # reshape the values to N_samples * N_factors and return
            expanded_expectation = QExp['E'][self.groups, :]
            expanded_lnE = QExp['lnE'][self.groups, :]
            # do we need to expand the variance as well -> not used I think
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

    def updateParameters(self, ix=None, ro=None):
        # TODO: add an if MuZ is in markov blanket ?
        tmp = self.markov_blanket["Z"].getExpectations()
        if 'ENN' in tmp:
            EZZ = tmp["ENN"]
        else:
            EZZ = tmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters().copy()
        Pa, Pb = P['a'], P['b']
        Qa, Qb = Q['a'], Q['b']

        ########################################################################
        # subset matrices for stochastic inference
        ########################################################################
        # TODO could this not be replaced by get_mini_batch ?
        if ix is None:
            ix = range(EZZ.shape[0])
        EZZ = EZZ[ix,:].copy()
        groups = self.groups[ix].copy()

        ########################################################################
        # compute the update
        ########################################################################
        par_up = self._updateParameters(Qa, Qb, Pa, Pb, EZZ, groups)

        ########################################################################
        # Do the asignment
        ########################################################################
        if ro is not None: # TODO have a default ro of 1 instead ? whats the overhead cost ?
        # TODO also change. do no deep copy but instead the same as in the other nodes
            par_up['Qa'] = ro * par_up['Qa'] + (1-ro) * self.Q.getParameters()['a']
            par_up['Qb'] = ro * par_up['Qb'] + (1-ro) * self.Q.getParameters()['b']
        self.Q.setParameters(a=par_up['Qa'], b=par_up['Qb'])

    def _updateParameters(self, Qa, Qb, Pa, Pb, EZZ, groups):
        for c in range(self.n_groups):
            mask = (groups == c)

            # coeff for stochastic inference
            n_batch = mask.sum()
            if n_batch == 0: continue
            n_total = self.n_per_group[c]
            coeff = n_total/n_batch


            Qa[c,:] = Pa[c,:] + 0.5 * n_total  # TODO should be precomputed
            Qb[c,:] = Pb[c,:] + 0.5 * coeff * EZZ[mask, :].sum(axis=0)

        return {'Qa': Qa, 'Qb': Qb}

    def calculateELBO(self):
        # Collect parameters and expectations
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.getExpectations()['E'], self.Q.getExpectations()['lnE']

        # Do the calculations
        lb_p = (Pa*s.log(Pb)).sum() - special.gammaln(Pa).sum() + ((Pa-1.)*QlnE).sum() - (Pb*QE).sum()
        lb_q = (Qa*s.log(Qb)).sum() - special.gammaln(Qa).sum() + ((Qa-1.)*QlnE).sum() - (Qb*QE).sum()

        return lb_p - lb_q
