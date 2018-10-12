from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
import scipy.special as special

from biofam.core.utils import *
from biofam.core import gpu_utils

# Import manually defined functions
from .variational_nodes import Gamma_Unobserved_Variational_Node

from biofam.core.distributions import *


class TauD_Node(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, groups, groups_dic, qE=None):

        self.groups = groups
        self.group_names = groups_dic
        self.N = len(self.groups)
        self.n_groups = len(np.unique(groups))

        assert self.n_groups == dim[0], "node dimension does not match number of groups"

        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)

    def precompute(self, options):
        self.lbconst = s.sum(self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']))
        gpu_utils.gpu_mode = options['gpu_mode']

        # compute number of sample per group
        self.n_per_group = np.zeros(self.n_groups)
        for c in range(self.n_groups):
            self.n_per_group[c] = (self.groups == c).sum()

        # update of Qa
        Y = self.markov_blanket["Y"].getExpectation()
        mask = ma.getmask(Y)
        Y = Y.data

        self.mini_batch = None

        # TODO is it ok to do that with stochastic ? here we kind of assume that missing values are evenly distributed
        self.Qa_pre = self.P.getParameters()['a'].copy()

        for g in range(self.n_groups):
            g_mask = (self.groups == g)
            Y_tmp = Y[g_mask, :]
            mask_tmp = mask[g_mask, :]
            self.Qa_pre[g,:] += (Y_tmp.shape[0] - mask_tmp.sum(axis=0))/2.

    def getExpectations(self, expand=True):
        QExp = self.Q.getExpectations()
        if expand:
            expanded_E = QExp['E'][self.groups, :]
            expanded_lnE = QExp['lnE'][self.groups, :]
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
        self.groups_batch = self.groups[ix]
        expanded_E = QExp['E'][self.groups_batch, :]
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
        #-----------------------------------------------------------------------
        # get full Expectations or minibatches
        #-----------------------------------------------------------------------
        Y = self.markov_blanket["Y"].get_mini_batch()
        Wtmp = self.markov_blanket["W"].getExpectations()
        Ztmp = self.markov_blanket["Z"].get_mini_batch()
        N = self.markov_blanket["Y"].dim[0]

        W, WW = Wtmp["E"], Wtmp["E2"]
        Z, ZZ = Ztmp["E"], Ztmp["E2"]

        # Collect parameters from the P distributions of this node
        P = self.P.getParameters()
        Pa, Pb = P['a'], P['b']

        #-----------------------------------------------------------------------
        # Masking
        #-----------------------------------------------------------------------
        mask = ma.getmask(Y)
        Y = Y.data
        Y[mask] = 0.

        #-----------------------------------------------------------------------
        # subsetting
        #-----------------------------------------------------------------------
        # TODO check that a deep copy is made for the subsetting
        groups = self.groups
        if ix is not None:
            groups = groups[ix]

        #-----------------------------------------------------------------------
        # make sure ro is not None
        #-----------------------------------------------------------------------
        # TODO fix that in BayesNet instead
        if ro is None:
            ro =1.

        #-----------------------------------------------------------------------
        # compute the update
        #-----------------------------------------------------------------------
        self._updateParameters(Y, W, WW, Z, ZZ, Pa, Pb, mask, ro, groups)

        ########################################################################
        # Do the asignment
        # ########################################################################
        # if ro is not None: # TODO have a default ro of 1 instead ? whats the overhead cost ?
        #     par_up['Qa'] = ro * par_up['Qa'] + (1-ro) * self.Q.getParameters()['a']
        #     par_up['Qb'] = ro * par_up['Qb'] + (1-ro) * self.Q.getParameters()['b']
        # self.Q.setParameters(a=par_up['Qa'], b=par_up['Qb'])

    def _updateParameters(self, Y, W, WW, Z, ZZ, Pa, Pb, mask, ro, groups):

        # Calculate terms for the update
        ZW =  gpu_utils.array(Z).dot(gpu_utils.array(W.T))
        ZW[mask] = 0.

        term1 = gpu_utils.square(gpu_utils.array(Y))

        term2 = gpu_utils.array(ZZ).dot(gpu_utils.array(WW.T))
        term2[mask] = 0

        term3 = - gpu_utils.dot(gpu_utils.square(gpu_utils.array(Z)),gpu_utils.square(gpu_utils.array(W)).T)
        term3[mask] = 0.
        term3 += gpu_utils.square(ZW)

        ZW *= gpu_utils.array(Y)  # WARNING ZW becomes ZWY
        term4 = 2.*ZW

        tmp = gpu_utils.asnumpy(term1 + term2 + term3 - term4)

        # weighted update (TODO check that ro default is 1)
        Qb = self.Q.getParameters()['b']
        Qb *= (1-ro)

        # TODO check that this is fine
        for g in range(self.n_groups):
            g_mask = (groups == g)

            n_batch = g_mask.sum()
            if n_batch == 0: continue
            n_total = self.n_per_group[g]
            coeff = n_total/n_batch

            # TODO check what happens if tmp[g_mask,:] is empty for a group
            Qb[g,:] += ro * (Pb[g,:] + .5 * coeff * tmp[g_mask,:].sum(axis=0))

        # TODO this should be completly useless
        Qa = self.Q.getParameters()['a']
        Qa *= (1 - ro)
        Qa += ro * self.Qa_pre


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
