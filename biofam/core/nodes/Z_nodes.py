from __future__ import division
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

class Z_Node(UnivariateGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, pmean, pvar, qmean, qvar, qE=None, qE2=None, idx_covariates=None):
        super().__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)

        self.N = self.dim[0]
        self.K = self.dim[1]

        # Define indices for covariates
        if idx_covariates is not None:
            self.covariates[idx_covariates] = True

    def precompute(self, options):
        # Precompute terms to speed up computation
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

    def removeFactors(self, idx, axis=1):
        super(Z_Node, self).removeFactors(idx, axis)

    def updateParameters(self):
        # Collect expectations from the markov blanket
        Y = self.markov_blanket["Y"].getExpectation()
        SWtmp = self.markov_blanket["W"].getExpectations()
        tau = self.markov_blanket["Tau"].getExpectation()
        latent_variables = self.getLvIndex()  # excluding covariates from the list of latent variables
        mask = [ma.getmask(Y[m]) for m in range(len(Y))]

        # Collect parameters from the prior or expectations from the markov blanket
        if "MuZ" in self.markov_blanket:
            Mu = self.markov_blanket['MuZ'].getExpectation()
        else:
            Mu = self.P.getParameters()["mean"]

        if "AlphaZ" in self.markov_blanket:
            Alpha = self.markov_blanket['AlphaZ'].getExpectation(expand=True)
        else:
            Alpha = 1./self.P.params['var']

        # Check dimensionality of Tau and expand if necessary (for Jaakola's bound only)
        for m in range(len(Y)):
            # Mask tau
            # tau[m] = ma.masked_where(ma.getmask(Y[m]), tau[m]) # important to keep this out of the loop to mask non-gaussian tau
            tau[m][mask[m]] = 0.
            # Mask Y
            Y[m] = Y[m].data
            Y[m][mask[m]] = 0.

        # Collect parameters from the P and Q distributions of this node
        Q = self.Q.getParameters()
        Qmean, Qvar = Q['mean'], Q['var']

        # import pdb; pdb.set_trace()
        M = len(Y)
        for k in latent_variables:
            foo = s.zeros((self.N,))
            bar = s.zeros((self.N,))
            for m in range(M):
                tau_cp = gpu_utils.array(tau[m])
                foo += gpu_utils.asnumpy(gpu_utils.dot(tau_cp, gpu_utils.array(SWtmp[m]["E2"][:, k])))

                bar_tmp1 = gpu_utils.array(SWtmp[m]["E"][:,k])

                # NOTE slow bit but hard to optimise
                # bar_tmp2 = - fast_dot(Qmean[:, s.arange(self.dim[1]) != k], SWtmp[m]["E"][:, s.arange(self.dim[1]) != k].T)
                tmp_cp1 = gpu_utils.array(Qmean[:, s.arange(self.dim[1]) != k])
                tmp_cp2 = gpu_utils.array(SWtmp[m]["E"][:, s.arange(self.dim[1]) != k].T)
                bar_tmp2 = - gpu_utils.dot(tmp_cp1, tmp_cp2)
                bar_tmp2 += gpu_utils.array(Y[m])
                bar_tmp2 *= tau_cp
                ##############################

                bar += gpu_utils.asnumpy(gpu_utils.dot(bar_tmp2, bar_tmp1))


            Qvar[:, k] = 1. / (Alpha[:, k] + foo)
            Qmean[:, k] = Qvar[:, k] * (bar + Alpha[:, k] * Mu[:, k])
        # Save updated parameters of the Q distribution
        self.Q.setParameters(mean=Qmean, var=Qvar)

    def calculateELBO(self):
        # Collect parameters and expectations of current node
        Qpar, Qexp = self.Q.getParameters(), self.Q.getExpectations()
        Qmean, Qvar = Qpar['mean'], Qpar['var']
        QE, QE2 = Qexp['E'], Qexp['E2']

        if "MuZ" in self.markov_blanket:
            PE, PE2 = self.markov_blanket['MuZ'].getExpectations()['E'], \
                      self.markov_blanket['MuZ'].getExpectations()['E2']
        else:
            PE, PE2 = self.P.getParameters()["mean"], s.zeros((self.N, self.dim[1]))

        if 'AlphaZ' in self.markov_blanket:
            Alpha = self.markov_blanket[
                'AlphaZ'].getExpectations(expand=True).copy()  # Notice that this Alpha is the ARD prior on Z, not on W.
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
        # lb_q = -(s.log(Qvar).sum() + self.N*self.dim[1])/2. # I THINK THIS IS WRONG BECAUSE SELF.DIM[1] ICNLUDES COVARIATES
        lb_q = -(s.log(Qvar).sum() + self.N * len(latent_variables)) / 2.

        return lb_p - lb_q

    def sample(self, dist='P'):
        if "MuZ" in self.markov_blanket:
            p_mean = self.markov_blanket['MuZ'].sample()
        else:
            p_mean = self.P.params['mean']

        if "AlphaZ" in self.markov_blanket:
            alpha = self.markov_blanket['AlphaZ'].sample()
            p_var = s.square(1. / alpha)
            #p_cov = s.diag(p_var)
            p_cov = [p_var[k] * np.eye(self.N) for k in range(self.K)]
        else:
            if "SigmaZ" in self.markov_blanket:
                p_cov = self.markov_blanket['SigmaZ'].sample()
            else:
                p_cov = self.P.params['cov']

        # simulating

        samp_tmp = []
        for i in range(self.dim[1]):
            # samp_tmp.append(s.random.multivariate_normal(p_mean[:, i], p_cov[i])) #does it yield the correct result for sparse input matrices ?
            if p_cov[i].__class__.__name__ == 'dia_matrix':
                samp_tmp.append(s.random.multivariate_normal(p_mean[:, i], p_cov[i].toarray()))  # inefficient
            elif p_cov[i].__class__.__name__ == 'ndarray':
                samp_tmp.append(s.random.multivariate_normal(p_mean[:, i], p_cov[i]))
            else:
                print("Not implemented yet")
                exit()

        # self.samp = s.array([tmp/tmp.std() for tmp in samp_tmp]).transpose()

        self.samp = s.array([tmp - tmp.mean() for tmp in samp_tmp]).transpose()
        #self.samp = np.array(samp_tmp).T

        return self.samp


class SZ_Node(BernoulliGaussian_Unobserved_Variational_Node):
    # TOO MANY ARGUMENTS, SHOULD WE USE **KWARGS AND *KARGS ONLY?
    def __init__(self, dim, pmean_T0, pmean_T1, pvar_T0, pvar_T1, ptheta, qmean_T0, qmean_T1, qvar_T0, qvar_T1, qtheta,
                 qEZ_T0=None, qEZ_T1=None, qET=None, idx_covariates=None):
        super().__init__(dim, pmean_T0, pmean_T1, pvar_T0, pvar_T1, ptheta, qmean_T0, qmean_T1, qvar_T0,
                                      qvar_T1, qtheta, qEZ_T0, qEZ_T1, qET)

        # Define indices for covariates
        if idx_covariates is not None:
            self.covariates[idx_covariates] = True

    def precompute(self, options):
        self.N = self.dim[0]
        self.K = self.dim[1]
        self.covariates = np.zeros(self.dim[1], dtype=bool)

        gpu_utils.gpu_mode = options['gpu_mode']

        self.factors_axis = 1

    def getLvIndex(self):
        # Method to return the index of the latent variables (without covariates)
        latent_variables = np.array(range(self.dim[1]))
        if any(self.covariates):
            # latent_variables = np.delete(latent_variables, latent_variables[self.covariates])
            latent_variables = latent_variables[~self.covariates]

        return latent_variables

    def updateParameters(self):
        # Collect expectations from other nodes
        # why .copy() ?

        Wtmp = [Wtmp_m for Wtmp_m in self.markov_blanket["W"].getExpectations()]
        W = [Wtmp_m["E"] for Wtmp_m in Wtmp]
        WW = [Wtmp_m["E2"] for Wtmp_m in Wtmp]

        tau = self.markov_blanket["Tau"].getExpectation()
        Y = self.markov_blanket["Y"].getExpectation()
        # TODO sort that out
        if "AlphaZ" in self.markov_blanket:
            alpha = self.markov_blanket["AlphaZ"].getExpectation(expand=True)
        else:
            alpha = 1./self.P.params['var_B1']
        thetatmp = self.markov_blanket['ThetaZ'].getExpectations(expand=True)
        theta_lnE, theta_lnEInv = thetatmp['lnE'], thetatmp['lnEInv']

        # print(theta_lnEInv)

        # Collect parameters and expectations from P and Q distributions of this node
        SZ = self.Q.getExpectations()["E"]
        Q = self.Q.getParameters()
        Qmean_T1, Qvar_T1, Qtheta = Q['mean_B1'], Q['var_B1'], Q['theta']

        M = len(Y)  # number of views

        # Mask matrices (excluding covariates from the list of latent variables)
        mask = [ma.getmask(Y[m]) for m in range(len(Y))]
        for m in range(M):
            Y[m] = Y[m].data
            Y[m][mask[m]] = 0.
            tau[m][mask[m]] = 0.

        # TODO : remove loops working with tensors ?
        # Update each latent variable in turn
        for k in range(self.dim[1]):
            # Calculate intermediate stept
            term1 = (theta_lnE - theta_lnEInv)[:, k]
            # TODO this is not right: alpha should be expended to full matrix before and used as full matrix
            # term2 = 0.5 * s.log(alpha[:, k]) should work
            # TODO modify everywhere else
            # GPU --------------------------------------------------------------
            # load variables on GPUs
            alphak_cp = gpu_utils.array(alpha[:,k])

            term2 = gpu_utils.asnumpy(0.5 * gpu_utils.log(alphak_cp))

            # term3 = 0.5*s.log(ma.dot(WW[:,k],tau) + alpha[k])
            # term3 = gpu_utils.zeros((self.N,))
            term4_tmp1 = gpu_utils.zeros((self.N,))
            term4_tmp2 = gpu_utils.zeros((self.N,))
            term4_tmp3 = gpu_utils.zeros((self.N,))
            for m in range(M):
                tau_cp = gpu_utils.array(tau[m])
                WWk_cp = gpu_utils.array(WW[m][:, k])
                Wk_cp = gpu_utils.array(W[m][:, k])

                term4_tmp1 += gpu_utils.dot(tau_cp * gpu_utils.array(Y[m]), Wk_cp) # good to modify # TODO critical time here  22 %
                term4_tmp2 += (tau_cp * gpu_utils.dot(gpu_utils.array(SZ[:, s.arange(self.dim[1]) != k]),
                                             (Wk_cp * gpu_utils.array(W[m][:, s.arange(self.dim[1]) != k].T)))).sum(axis=1)  # good to modify  # TODO critical time here  37 %
                term4_tmp3 += gpu_utils.dot(tau_cp, WWk_cp) # TODO critical time here  3 %

            # term4 = 0.5*s.divide((term4_tmp1-term4_tmp2)**2,term4_tmp3)
            term3 =gpu_utils.asnumpy(0.5 * gpu_utils.log(term4_tmp3 + alphak_cp))

            term4_tmp3 += alphak_cp
            term4 = 0.5 * gpu_utils.divide(gpu_utils.square(term4_tmp1 - term4_tmp2), term4_tmp3)  # good to modify, awsnt checked numerically
            term4 = gpu_utils.asnumpy(term4)

            # Update S
            # NOTE there could be some precision issues in T --> loads of 1s in result
            Qtheta[:, k] = 1. / (1. + s.exp(-(term1 + term2 - term3 + term4)))
            Qtheta[:,k] = np.nan_to_num(Qtheta[:,k])

            # Update Z
            Qvar_T1[:, k] = 1. / gpu_utils.asnumpy(term4_tmp3)
            Qmean_T1[:, k] = Qvar_T1[:, k] * gpu_utils.asnumpy(term4_tmp1 - term4_tmp2)

            # Update Expectations for the next iteration
            SZ[:, k] = Qtheta[:, k] * Qmean_T1[:, k]

        # Save updated parameters of the Q distribution
        self.Q.setParameters(mean_B0=s.zeros((self.N, self.dim[1])), var_B0=1. / alpha,
                             mean_B1=Qmean_T1, var_B1=Qvar_T1, theta=Qtheta)

    def calculateELBO(self):

        # Collect parameters and expectations
        Qpar, Qexp = self.Q.getParameters(), self.Q.getExpectations()
        T, ZZ = Qexp["EB"], Qexp["ENN"]
        Qvar = Qpar['var_B1']
        theta = self.markov_blanket['ThetaZ'].getExpectations(expand=True)

        # Get ARD sparsity or prior variance
        if "AlphaZ" in self.markov_blanket:
            alpha = self.markov_blanket['AlphaZ'].getExpectations(expand=True).copy()
        else:
            alpha = dict()
            alpha['E'] = 1./self.P.params['var_B1']
            alpha['lnE'] = s.log(1./self.P.params['var_B1'])

        # This ELBO term contains only cross entropy between Q and P, and entropy of Q. So the covariates should not intervene at all
        latent_variables = self.getLvIndex()
        alpha["E"], alpha["lnE"] = alpha["E"][:,latent_variables], alpha["lnE"][:,latent_variables]
        T, ZZ = T[:, latent_variables], ZZ[:, latent_variables]
        Qvar = Qvar[:, latent_variables]
        theta['lnE'] = theta['lnE'][:,latent_variables]
        theta['lnEInv'] = theta['lnEInv'][:,latent_variables]

        # Calculate ELBO for Z
        lb_pz = (alpha["lnE"].sum() - s.sum(alpha["E"] * ZZ)) / 2.
        lb_qz = -0.5 * self.dim[1] * self.N - 0.5 * (
                    T * s.log(Qvar) + (1. - T) * s.log(1. / alpha["E"])).sum()  # IS THE FIRST CONSTANT TERM CORRECT???
        lb_z = lb_pz - lb_qz

        # Calculate ELBO for T
        lb_pt = T * theta['lnE'] + (1. - T) * theta['lnEInv']
        lb_qt = T * s.log(T) + (1. - T) * s.log(1. - T)
        lb_pt[s.isnan(lb_pt)] = 0.
        lb_qt[s.isnan(lb_qt)] = 0.
        lb_t = s.sum(lb_pt) - s.sum(lb_qt)

        return lb_z + lb_t

    def sample(self, dist='P'):
        # get necessary parameters
        mu_z_hat = self.P.getParameters()['mean_B1']
        mu_z_hat = s.ones(self.dim) * mu_z_hat

        theta = self.markov_blanket['ThetaZ'].sample()
        #if theta.shape != mu_z_hat.shape: #if theta had already mu_z_hat shape, then the sampling above would be wrong
        theta = s.repeat(theta[None, :], mu_z_hat.shape[0], 0)

        alpha = self.markov_blanket['AlphaZ'].sample()
        #if alpha.shape[0] == 1:
        #    alpha = s.repeat(alpha[:], self.dim[1], axis=0)
        #if alpha.shape != mu_z_hat.shape:
        #    alpha = s.repeat(alpha[None, :], self.dim[0], axis=0)

        # simulate
        bernoulli_t = s.random.binomial(1, theta)

        z_hat = np.array([s.random.normal(mu_z_hat[:, i], math.sqrt(1./alpha[i])) for i in range(mu_z_hat.shape[1])]).T
        #z_hat = s.random.normal(mu_z_hat, np.sqrt(1. / alpha))

        self.samp = bernoulli_t * z_hat
        self.samp[:, self.covariates] = self.getExpectation()[:, self.covariates]

        return self.samp


