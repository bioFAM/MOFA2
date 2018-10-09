from __future__ import division
from copy import deepcopy
import numpy.ma as ma
import numpy as np
import scipy as s
from copy import deepcopy
import math

from biofam.core.utils import *
from biofam.core import gpu_utils


# Import manually defined functions
from .variational_nodes import BernoulliGaussian_Unobserved_Variational_Node
from .variational_nodes import UnivariateGaussian_Unobserved_Variational_Node
from .variational_nodes import UnivariateGaussian_Unobserved_Variational_Node_with_MultivariateGaussian_Prior

# TODO make more general to handle both cases with and without SL in the Markov Blanket
class W_Node(UnivariateGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, pmean, pvar, qmean, qvar, qE=None, qE2=None, idx_covariates=None):
        super().__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)

        # Define indices for covariates
        if idx_covariates is not None:
            self.covariates[idx_covariates] = True

    def precompute(self, options):
        # Precompute terms to speed up computation
        self.D = self.dim[0]
        self.K = self.dim[1]
        self.covariates = np.zeros(self.dim[1], dtype=bool)
        self.factors_axis = 1
        gpu_utils.gpu_mode = options['gpu_mode']

    def getLvIndex(self):
        # Method to return the index of the latent variables (without covariates)
        latent_variables = np.array(range(self.dim[1]))
        if any(self.covariates):
            # latent_variables = np.delete(latent_variables, latent_variables[self.covariates])
            latent_variables = latent_variables[~self.covariates]
        return latent_variables

    def removeFactors(self, idx):
        super().removeFactors(idx, axis=1)

    def updateParameters(self):
        # print(self.getExpectations())
        # Collect expectations from the markov blanket
        Y = self.markov_blanket["Y"].getExpectation()
        SZtmp = self.markov_blanket["Z"].getExpectations()
        tau = self.markov_blanket["Tau"].getExpectation()
        latent_variables = self.getLvIndex() # excluding covariates from the list of latent variables
        mask = ma.getmask(Y)

        # Collect parameters from the prior or expectations from the markov blanket
        if "MuW" in self.markov_blanket:
            Mu = self.markov_blanket['MuW'].getExpectation()
        else:
            Mu = self.P.getParameters()["mean"]

        if "AlphaW" in self.markov_blanket:
            Alpha = self.markov_blanket['AlphaW'].getExpectation(expand=True)
        else:
            Alpha = 1./self.P.params['var']

        # Mask tau
        tau[mask] = 0.
        # Mask Y
        Y = Y.data
        Y[mask] = 0.

        # Collect parameters from the P and Q distributions of this node
        Q = self.Q.getParameters()
        Qmean, Qvar = Q['mean'], Q['var']

        for k in latent_variables:

            foo = s.dot(SZtmp["E2"][:,k], tau)

            # GPU bit --------------------------------------------------------
            bar_tmp1 = gpu_utils.array(SZtmp["E"][:,k])

            bar_tmp2 = - gpu_utils.dot(gpu_utils.array(SZtmp["E"][:,s.arange(self.dim[1])!=k]),
                               gpu_utils.array(Qmean[:,s.arange(self.dim[1])!=k].T))
            bar_tmp2 += gpu_utils.array(Y)
            bar_tmp2 *= gpu_utils.array(tau)

            bar = gpu_utils.asnumpy(gpu_utils.dot(bar_tmp1, bar_tmp2))
            # ----------------------------------------------------------------

            Qvar[:,k] = 1./(Alpha[:,k]+foo)
            Qmean[:,k] = Qvar[:,k] * (bar + Alpha[:,k]*Mu[:,k])

        # Save updated parameters of the Q distribution
        self.Q.setParameters(mean=Qmean, var=Qvar)

    def calculateELBO(self):

        # Collect parameters and expectations of current node
        Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()
        Qmean, Qvar = Qpar['mean'], Qpar['var']
        QE, QE2 = Qexp['E'],Qexp['E2']

        if "MuW" in self.markov_blanket:
            PE, PE2 = self.markov_blanket['MuW'].getExpectations()['E'], self.markov_blanket['MuW'].getExpectations()['E2']
        else:
            PE, PE2 = self.P.getParameters()["mean"], s.zeros((self.D,self.dim[1]))

        # if "SigmaAlphaW" in self.markov_blanket:
        #     name_alpha = "SigmaAlphaW"
        # else:
        #     name_alpha = "AlphaW"

        if 'AlphaW' in self.markov_blanket:
            Alpha = self.markov_blanket["AlphaW"].getExpectations(expand=True).copy() # Notice that this Alpha is the ARD prior on Z, not on W.
        else:
            Alpha = dict()
            Alpha['E'] = 1./self.P.params['var']
            Alpha['lnE'] = s.log(1./self.P.params['var'])

        # This ELBO term contains only cross entropy between Q and P,and entropy of Q. So the covariates should not intervene at all
        latent_variables = self.getLvIndex()
        Alpha["E"], Alpha["lnE"] = Alpha["E"][:, latent_variables], Alpha["lnE"][:, latent_variables]
        Qmean, Qvar = Qmean[:, latent_variables], Qvar[:, latent_variables]
        PE, PE2 = PE[:, latent_variables], PE2[:, latent_variables]
        QE, QE2 = QE[:, latent_variables], QE2[:, latent_variables]

        # compute term from the exponential in the Gaussian
        tmp1 = 0.5 * QE2 - PE * QE + 0.5 * PE2
        tmp1 = -(tmp1 * Alpha['E']).sum()

        # compute term from the precision factor in front of the Gaussian
        tmp2 = 0.5 * Alpha["lnE"].sum()

        lb_p = tmp1 + tmp2
        # lb_q = -(s.log(Qvar).sum() + self.D*self.dim[1])/2. # I THINK THIS IS WRONG BECAUSE SELF.DIM[1] ICNLUDES COVARIATES
        lb_q = -(s.log(Qvar).sum() + self.D * len(latent_variables)) / 2.
        return lb_p-lb_q


    def sample(self, dist='P'):
        if "MuW" in self.markov_blanket:
            p_mean = self.markov_blanket['MuW'].sample()
        else:
            p_mean = self.P.params['mean']

        if "SigmaAlphaW" in self.markov_blanket:
            if self.markov_blanket["SigmaAlphaW"].__class__.__name__ == "AlphaW_Node_mk":
                alpha = self.markov_blanket['SigmaAlphaW'].sample()
                p_var = s.square(1./alpha)
                #p_cov = s.diag(p_var)
                p_cov = [p_var[k]*np.eye(self.D) for k in range(self.K)]
            else:
                p_cov = self.markov_blanket['SigmaAlphaW'].sample()
        elif "AlphaW" in self.markov_blanket:
            alpha = self.markov_blanket['AlphaW'].sample()
            p_var = s.square(1. / alpha)
            #p_cov = s.diag(p_var)
            p_cov = [p_var[k] * np.eye(self.D) for k in range(self.K)]
        else:
            p_cov = self.P.params['cov']

        # simulating

        samp_tmp = []
        for i in range(self.dim[1]):
            #samp_tmp.append(s.random.multivariate_normal(p_mean[:, i], p_cov[i])) #does it yield the correct result for sparse input matrices ?
            if p_cov[i].__class__.__name__ == 'dia_matrix':
                samp_tmp.append(s.random.multivariate_normal(p_mean[:, i], p_cov[i].toarray())) #inefficient
            elif p_cov[i].__class__.__name__ == 'ndarray':
                samp_tmp.append(s.random.multivariate_normal(p_mean[:, i], p_cov[i]))
            else:
                print("Not implemented yet")
                exit()

        # self.samp = s.array([tmp/tmp.std() for tmp in samp_tmp]).transpose()
        self.samp = s.array([tmp - tmp.mean() for tmp in samp_tmp]).transpose()

        return self.samp

class SW_Node(BernoulliGaussian_Unobserved_Variational_Node):
    # TOO MANY ARGUMENTS, SHOULD WE USE **KWARGS AND *KARGS ONLY?
    # def __init__(self, dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0=None, qEW_S1=None, qES=None):
    def __init__(self, dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0=None, qEW_S1=None, qES=None):
        super().__init__(dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0, qEW_S1, qES)

    def precompute(self, options):
        self.D = self.dim[0]
        self.factors_axis = 1
        gpu_utils.gpu_mode = options['gpu_mode']

    def updateParameters(self):
        # Collect expectations from other nodes
        Ztmp = self.markov_blanket["Z"].getExpectations()
        Z,ZZ = Ztmp["E"],Ztmp["E2"]
        tau = self.markov_blanket["Tau"].getExpectation(expand=True)  # TODO might be worth expanding only once when updating expectation
        Y = self.markov_blanket["Y"].getExpectation()
        if "AlphaW" in self.markov_blanket:
            alpha = self.markov_blanket["AlphaW"].getExpectation(expand=True)
        else:
            alpha = 1./self.P.params['var_B1']
        thetatmp = self.markov_blanket["ThetaW"].getExpectations(expand=True)
        theta_lnE, theta_lnEInv  = thetatmp['lnE'], thetatmp['lnEInv']


        mask = ma.getmask(Y)

        # Collect parameters and expectations from P and Q distributions of this node
        SW = self.Q.getExpectations()["E"]
        Q = self.Q.getParameters()
        Qmean_S1, Qvar_S1, Qtheta = Q['mean_B1'], Q['var_B1'], Q['theta']

        # Mask matrices
        Y = Y.data
        Y[mask] = 0.
        tau[mask] = 0.

        # precompute terms usful for all k
        tauYT = (gpu_utils.array(tau)*gpu_utils.array(Y)).T

        # Update each latent variable in turn
        for k in range(self.dim[1]):
            # Calculate intermediate steps
            term1 = (theta_lnE-theta_lnEInv)[:,k]

            # GPU --------------------------------------------------------------
            # Variables used in multiple operations snhould be loaded on GPU only once
            Zk_cp = gpu_utils.array(Z[:,k])
            tau_cp = gpu_utils.array(tau)
            ZZk_cp = gpu_utils.array(ZZ[:,k])
            alphak_cp = gpu_utils.array(alpha[:,k])

            term2 = gpu_utils.asnumpy(0.5*gpu_utils.log(alphak_cp))
            # term3 = 0.5*s.log(fast_dot(ZZ[:,k], tau) + alpha[:,k])
            term3 = gpu_utils.asnumpy(0.5*gpu_utils.log(gpu_utils.dot(ZZk_cp, tau_cp) + alphak_cp))

            term4_tmp1 = gpu_utils.dot(tauYT, Zk_cp)
            # term4_tmp1 = fast_dot(tauYT,Z[:,k])

            term4_tmp2_1 = gpu_utils.array(SW[:,s.arange(self.dim[1])!=k].T)
            term4_tmp2_2 = (Zk_cp * gpu_utils.array(Z[:,s.arange(self.dim[1])!=k]).T).T
            term4_tmp2 = gpu_utils.dot(term4_tmp2_2, term4_tmp2_1)
            # term4_tmp2 = fast_dot(term4_tmp2_2, term4_tmp2_1)
            term4_tmp2 *= tau_cp  # most expensive bit
            term4_tmp2 = term4_tmp2.sum(axis=0)

            term4_tmp3 = gpu_utils.dot(ZZk_cp.T,tau_cp) + alphak_cp # good to modify (I REPLACE MA.DOT FOR S.DOT, IT SHOULD BE SAFE )
            # term4_tmp3 = fast_dot(ZZ[:,k].T,tau) + alpha[:,k]


            term4 = gpu_utils.asnumpy(0.5*gpu_utils.divide(gpu_utils.square(term4_tmp1-term4_tmp2),term4_tmp3)) # good to modify, awsnt checked numerically

            # ------------------------------------------------------------------

            # Update S
            # NOTE there could be some precision issues in S --> loads of 1s in result
            Qtheta[:,k] = 1./(1.+s.exp(-(term1+term2-term3+term4)))

            # Update W
            Qvar_S1[:,k] = gpu_utils.asnumpy(1./term4_tmp3)
            Qmean_S1[:,k] = Qvar_S1[:,k]*gpu_utils.asnumpy(term4_tmp1-term4_tmp2)

            # Update Expectations for the next iteration
            SW[:,k] = Qtheta[:,k] * Qmean_S1[:,k]

        # Save updated parameters of the Q distribution
        self.Q.setParameters(mean_B0=s.zeros((self.D,self.dim[1])), var_B0=1./alpha, mean_B1=Qmean_S1, var_B1=Qvar_S1, theta=Qtheta )

    def calculateELBO(self):
        # Collect parameters and expectations
        Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()
        S,WW = Qexp["EB"], Qexp["ENN"]
        Qvar = Qpar['var_B1']

        theta = self.markov_blanket["ThetaW"].getExpectations(expand=True)


        # Get ARD sparsity or prior variance
        if "AlphaW" in self.markov_blanket:
            alpha = self.markov_blanket["AlphaW"].getExpectations(expand=True).copy()  # TODO is this copy necessarry ??
        else:
            alpha = dict()
            alpha['E'] = 1./self.P.params['var_B1']
            alpha['lnE'] = s.log(1./self.P.params['var_B1'])


        # Calculate ELBO for W
        lb_pw = (alpha["lnE"].sum() - s.sum(alpha["E"]*WW))/2.
        lb_qw = -0.5*self.dim[1]*self.D - 0.5*(S*s.log(Qvar) + (1.-S)*s.log(1./alpha["E"])).sum() # IS THE FIRST CONSTANT TERM CORRECT???
        # #NOT SURE ABOUT THE FORMULA for lb_qw (brackets of expectation propagating inside the log ?)
        lb_w = lb_pw - lb_qw

        # Calculate ELBO for S
        lb_ps = S*theta['lnE'] + (1.-S)*theta['lnEInv']
        lb_qs = S*s.log(S) + (1.-S)*s.log(1.-S)
        lb_ps[s.isnan(lb_ps)] = 0.
        lb_qs[s.isnan(lb_qs)] = 0.
        lb_s = s.sum(lb_ps) - s.sum(lb_qs)

        return lb_w + lb_s

    def sample(self, dist='P'):
        # get necessary parameters
        mu_w_hat = self.P.getParameters()['mean_B1']
        mu_w_hat = s.ones(self.dim) * mu_w_hat

        theta = self.markov_blanket["ThetaW"].sample()
        #if theta.shape != mu_w_hat.shape:  #if theta had already mu_z_hat shape, then the sampling above would be wrong
        theta = s.repeat(theta[None,:],mu_w_hat.shape[0],0)

        #if "AlphaW" in self.markov_blanket:
        alpha = self.markov_blanket["AlphaW"].sample()
            #if alpha.shape[0] == 1:
            #    alpha = s.repeat(alpha[:], self.dim[1], axis=0)
            #if alpha.shape != mu_w_hat.shape:
            #    alpha = s.repeat(alpha[None,:], self.dim[0], axis=0)
        #else:
        #    print("Not implemented")
        #    exit()

        # simulate
        bernoulli_s = s.random.binomial(1, theta)

        w_hat = np.array([s.random.normal(mu_w_hat[:, i], math.sqrt(1./alpha[i])) for i in range(mu_w_hat.shape[1])]).T
        #w_hat = np.random.normal(mu_w_hat, np.sqrt(1./alpha))

        self.samp = bernoulli_s * w_hat
        return self.samp
