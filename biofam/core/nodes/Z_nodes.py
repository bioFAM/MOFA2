from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
from copy import deepcopy
import math
from biofam.core.utils import *

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

    def precompute(self):
        # Precompute terms to speed up computation
        self.covariates = np.zeros(self.dim[1], dtype=bool)
        self.factors_axis = 1

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

        M = len(Y)
        for k in latent_variables:
            foo = s.zeros((self.N,))
            bar = s.zeros((self.N,))
            for m in range(M):
                foo += np.dot(tau[m], SWtmp[m]["E2"][:, k])

                bar_tmp1 = SWtmp[m]["E"][:,k]

                # NOTE slow bit but hard to optimise
                # bar_tmp2 = - fast_dot(Qmean[:, s.arange(self.dim[1]) != k], SWtmp[m]["E"][:, s.arange(self.dim[1]) != k].T)
                bar_tmp2 = - s.dot(Qmean[:, s.arange(self.dim[1]) != k], SWtmp[m]["E"][:, s.arange(self.dim[1]) != k].T)
                bar_tmp2 += Y[m]
                bar_tmp2 *= tau[m]
                ##############################

                bar += np.dot(bar_tmp2, bar_tmp1)

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
        self.precompute()

        # Define indices for covariates
        if idx_covariates is not None:
            self.covariates[idx_covariates] = True

    def precompute(self):
        self.N = self.dim[0]
        self.K = self.dim[1]
        self.covariates = np.zeros(self.dim[1], dtype=bool)

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
        # import pdb; pdb.set_trace()
        for k in range(self.dim[1]):
            # Calculate intermediate stept
            term1 = (theta_lnE - theta_lnEInv)[:, k]
            # TODO this is not right: alpha should be expended to full matrix before and used as full matrix
            # term2 = 0.5 * s.log(alpha[:, k]) should work
            # TODO modify everywhere else
            term2 = 0.5 * s.log(alpha[:,k])

            # term3 = 0.5*s.log(ma.dot(WW[:,k],tau) + alpha[k])
            term3 = 0.5 * s.log(np.sum(np.array([s.dot(tau[m], WW[m][:, k]) for m in range(M)]), axis=0) + alpha[:,k])  # good to modify # TODO critical time here  2 %

            # term4_tmp1 = ma.dot((tau*Y).T,W[:,k]).data
            term4_tmp1 = np.sum(np.array([s.dot(tau[m] * Y[m], W[m][:, k]) for m in range(M)]), axis=0)  # good to modify # TODO critical time here  22 %

            # term4_tmp2 = ( tau * s.dot((W[:,k]*W[:,s.arange(self.dim[1])!=k].T).T, SZ[:,s.arange(self.dim[1])!=k].T) ).sum(axis=0)
            term4_tmp2 = np.sum(np.array([(tau[m] * s.dot(SZ[:, s.arange(self.dim[1]) != k], (W[m][:, k] * W[m][:, s.arange(self.dim[1]) != k].T))).sum(axis=1) for m in range(M)]), axis=0)  # good to modify  # TODO critical time here  37 %

            # term4_tmp3 = s.dot(WW[:,k].T,tau) + alpha[k] # good to modify (I REPLACE MA.DOT FOR S.DOT, IT SHOULD BE SAFE )
            term4_tmp3 = np.sum(np.array([s.dot(tau[m], WW[m][:, k]) for m in range(M)]), axis=0) + alpha[:,k]  # TODO critical time here  3 %

            # term4 = 0.5*s.divide((term4_tmp1-term4_tmp2)**2,term4_tmp3)
            term4 = 0.5 * s.divide(s.square(term4_tmp1 - term4_tmp2), term4_tmp3)  # good to modify, awsnt checked numerically

            # Update S
            # NOTE there could be some precision issues in T --> loads of 1s in result
            Qtheta[:, k] = 1. / (1. + s.exp(-(term1 + term2 - term3 + term4)))
            # import pdb; pdb.set_trace()
            Qtheta[:,k] = np.nan_to_num(Qtheta[:,k])

            # Update Z
            Qvar_T1[:, k] = 1. / term4_tmp3
            Qmean_T1[:, k] = Qvar_T1[:, k] * (term4_tmp1 - term4_tmp2)

            # Update Expectations for the next iteration
            SZ[:, k] = Qtheta[:, k] * Qmean_T1[:, k]

        # Save updated parameters of the Q distribution
        # import pdb; pdb.set_trace()
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



class MuZ_Node(UnivariateGaussian_Unobserved_Variational_Node):
    """ """

    def __init__(self, pmean, pvar, qmean, qvar, clusters, n_Z, cluster_dic=None, qE=None, qE2=None):
        # For now clusters have to be integers from 0 to n_clusters
        # compute dim from numbers of clusters (n_clusters * Z)
        self.clusters = clusters
        self.N = len(self.clusters)
        self.n_clusters = len(np.unique(clusters))
        dim = (self.n_clusters, n_Z)
        self.factors_axis = 1
        super(Cluster_Node_Gaussian, self).__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE,
                                                    qE2=qE2)

    def getExpectations(self):
        # reshape the values to N_samples * N_factors and return
        QExp = self.Q.getExpectations()
        expanded_expectation = QExp['E'][self.clusters, :]
        expanded_E2 = QExp['E2'][self.clusters, :]
        # do we need to expand the variance as well -> not used I think
        return {'E': expanded_expectation, 'E2': expanded_E2}

    def updateParameters(self):
        Ppar = self.P.getParameters()
        Z = self.markov_blanket['Z'].Q.getExpectation()

        if "AlphaZ" in self.markov_blanket:
            Alpha = self.markov_blanket[
                'AlphaZ'].getExpectation().copy()  # Notice that this Alpha is the ARD prior on Z, not on W.
        elif "SigmaZ" in self.markov_blanket:
            Sigma = self.markov_blanket['SigmaZ'].getExpectation().copy()
        else:
            Sigma = self.markov_blanket['Z'].P.getParameters()["cov"]

        Qmean, Qvar = self.Q.getParameters()['mean'], self.Q.getParameters()['var']

        # update of the variance

        for c in range(self.n_clusters):
            mask = (self.clusters == c)
            if "AlphaZ" in self.markov_blanket:
                tmp = (Alpha[mask, :]).sum(axis=0)
            else:
                tmp = np.matrix.trace(Sigma[:, mask, mask])
            Qvar[c, :] = tmp
        Qvar += 1. / Ppar['var']
        Qvar = 1. / Qvar

        # update of the mean

        if "AlphaZ" in self.markov_blanket:
            tmp = Z * Alpha
        else:
            # TODO : check if we should ask the mask l. 462
            tmp = np.zeros(self.dim)
            for k in range(self.dim[1]):
                tmp[:, k] = np.dot(Sigma[k, :, :] - np.diag(Sigma[k, :, :]), Z[:, k])

        for c in range(self.n_clusters):
            mask = (self.clusters == c)
            tmp = (tmp[mask, :]).sum(axis=0)
            Qmean[c, :] = tmp
        Qmean = Qmean + Ppar['mean'] / Ppar['var']
        Qmean *= Qvar

        self.Q.setParameters(mean=Qmean, var=Qvar)

    def calculateELBO(self):
        PParam = self.P.getParameters()
        PVar, Pmean = PParam['var'], PParam['mean']

        QExp = self.Q.getExpectations()
        QE2, QE = QExp['E2'], QExp['E']

        Qvar = self.Q.getParameters()['var']

        # Cluster terms corresponding to covariates should not intervene
        # filtering the covariates out
        latent_variables = self.markov_blanket['Z'].getLvIndex()
        PVar, Pmean = PVar[:, latent_variables], Pmean[:, latent_variables]
        QE2, QE = QE2[:, latent_variables], QE[:, latent_variables]
        Qvar = Qvar[:, latent_variables]

        # minus cross entropy
        tmp = -(0.5 * s.log(PVar)).sum()
        tmp2 = - ((0.5 / PVar) * (QE2 - 2. * QE * Pmean + Pmean ** 2.)).sum()

        # entropy of Q
        tmp3 = 0.5 * (s.log(Qvar)).sum()
        tmp3 += 0.5 * self.dim[0] * len(latent_variables)

        return tmp + tmp2 + tmp3


# Z node with multivariate prior
# class Z_Node(UnivariateGaussian_Unobserved_Variational_Node_with_MultivariateGaussian_Prior):
#     def __init__(self, dim, pmean, pcov, qmean, qvar, qE=None, qE2=None, idx_covariates=None, precompute_pcovinv=True):
#         super().__init__(dim=dim, pmean=pmean, pcov=pcov, qmean=qmean, qvar=qvar, axis_cov=0, qE=qE, qE2=qE2)
#
#         self.precompute_pcovinv = precompute_pcovinv
#
#         self.N = self.dim[0]
#         self.K = self.dim[1]
#
#         # Define indices for covariates
#         if idx_covariates is not None:
#             self.covariates[idx_covariates] = True
#
#     def precompute(self):
#         # Precompute terms to speed up computation
#         self.covariates = np.zeros(self.dim[1], dtype=bool)
#         self.factors_axis = 1
#
#         if self.precompute_pcovinv:
#             p_cov = self.P.params["cov"]
#
#             self.p_cov_inv = []
#             self.p_cov_inv_diag = []
#
#             for k in range(self.K):
#                 if p_cov[k].__class__.__name__ == 'dia_matrix':
#                     diag_inv = 1 / p_cov[k].data
#                     diag_inv = diag_inv.flatten()
#                     inv = np.diag(diag_inv)
#                 elif p_cov[k].__class__.__name__ == 'ndarray':
#                     diagonal = np.diagonal(p_cov[k])
#                     if np.all(np.diag(diagonal) == p_cov[k]):
#                         diag_inv = 1. / diagonal
#                         inv = np.diag(diag_inv)
#                     else:
#                         inv = np.linalg.inv(p_cov[k])
#                         diag_inv = np.diagonal(inv)
#                 else:
#                     #TODO : deal with sparse non diagonal input matrices as pcov
#                     print("Not implemented yet")
#                     exit()
#
#                 self.p_cov_inv.append(inv)
#                 self.p_cov_inv_diag.append(diag_inv)
#
#         else:
#             self.p_cov_inv = None
#             self.p_cov_inv_diag = None
#
#     def getLvIndex(self):
#         # Method to return the index of the latent variables (without covariates)
#         latent_variables = np.array(range(self.dim[1]))
#         if any(self.covariates):
#             # latent_variables = np.delete(latent_variables, latent_variables[self.covariates])
#             latent_variables = latent_variables[~self.covariates]
#         return latent_variables
#
#     def removeFactors(self, idx, axis=1):
#         super(Z_Node, self).removeFactors(idx, axis)
#
#         if self.p_cov_inv is not None:
#             for i in idx :
#                 del self.p_cov_inv[i]
#                 del self.p_cov_inv_diag[i]
#
#     def updateParameters(self):
#         # Collect expectations from the markov blanket
#         Y = self.markov_blanket["Y"].getExpectation()
#         SWtmp = self.markov_blanket["W"].getExpectations()
#         tau = self.markov_blanket["Tau"].getExpectation()
#         latent_variables = self.getLvIndex()  # excluding covariates from the list of latent variables
#         mask = [ma.getmask(Y[m]) for m in range(len(Y))]
#
#         # Collect parameters from the prior or expectations from the markov blanket
#         if "MuZ" in self.markov_blanket:
#             Mu = self.markov_blanket['MuZ'].getExpectation()
#         else:
#             Mu = self.P.getParameters()["mean"]
#
#         if "AlphaZ" in self.markov_blanket:
#             Alpha = self.markov_blanket['AlphaZ'].getExpectation(expand=True)
#
#         else:
#             if "SigmaZ" in self.markov_blanket:
#                 Sigma = self.markov_blanket['SigmaZ'].getExpectations()
#                 p_cov_inv = Sigma['inv']
#                 p_cov_inv_diag = Sigma['inv_diag']
#             else:
#                 p_cov_inv = self.p_cov_inv
#                 p_cov_inv_diag = self.p_cov_inv_diag
#
#         # Check dimensionality of Tau and expand if necessary (for Jaakola's bound only)
#         for m in range(len(Y)):
#             # Mask tau
#             # tau[m] = ma.masked_where(ma.getmask(Y[m]), tau[m]) # important to keep this out of the loop to mask non-gaussian tau
#             tau[m][mask[m]] = 0.
#             # Mask Y
#             Y[m] = Y[m].data
#             Y[m][mask[m]] = 0.
#
#         # Collect parameters from the P and Q distributions of this node
#         Q = self.Q.getParameters()
#         Qmean, Qvar = Q['mean'], Q['var']
#
#         M = len(Y)
#         for k in latent_variables:
#             foo = s.zeros((self.N,))
#             bar = s.zeros((self.N,))
#             for m in range(M):
#                 foo += np.dot(tau[m], SWtmp[m]["E2"][:, k])
#
#                 bar_tmp1 = SWtmp[m]["E"][:,k]
#
#                 # NOTE slow bit but hard to optimise
#                 # bar_tmp2 = - fast_dot(Qmean[:, s.arange(self.dim[1]) != k], SWtmp[m]["E"][:, s.arange(self.dim[1]) != k].T)
#                 bar_tmp2 = - s.dot(Qmean[:, s.arange(self.dim[1]) != k], SWtmp[m]["E"][:, s.arange(self.dim[1]) != k].T)
#                 bar_tmp2 += Y[m]
#                 bar_tmp2 *= tau[m]
#                 ##############################
#
#                 bar += np.dot(bar_tmp2, bar_tmp1)
#
#             if "AlphaZ" in self.markov_blanket:
#                 Qvar[:, k] = 1. / (Alpha[:, k] + foo)
#                 Qmean[:, k] = Qvar[:, k] * (bar + Alpha[:, k] * Mu[:, k])
#
#             else:
#                 Qvar[:, k] = 1. / (foo + p_cov_inv_diag[k])
#
#                 if self.P.params["cov"][k].__class__.__name__ == 'dia_matrix':
#                     Qmean[:, k] = Qvar[:, k] * bar
#                 else:
#                     tmp = p_cov_inv[k] - p_cov_inv_diag[k] * s.eye(self.N)
#                     for n in range(self.N):
#                         Qmean[n, k] = Qvar[n, k] * (bar[n] + np.dot(tmp[n, :], Mu[:, k] - Qmean[:, k]))
#
#         # Save updated parameters of the Q distribution
#         self.Q.setParameters(mean=Qmean, var=Qvar)
#
#     # TODO, problem here is that we need to make sure k is in the latent variables first
#     def calculateELBO_k(self, k):
#         '''Compute the ELBO for factor k in absence of Alpha node in the markov blanket of Z'''
#         Qpar, Qexp = self.Q.getParameters(), self.Q.getExpectations()
#         Qmean, Qvar = Qpar['mean'], Qpar['var']
#         QE, QE2 = Qexp['E'], Qexp['E2']
#
#         if "SigmaZ" in self.markov_blanket:
#             Sigma = self.markov_blanket['SigmaZ'].getExpectations()
#             p_cov = Sigma['cov']
#             p_cov_inv = Sigma['inv']
#             p_cov_inv_diag = Sigma['inv_diag']
#         else:
#             p_cov = self.P.params['cov']
#             p_cov_inv = self.p_cov_inv
#             p_cov_inv_diag = self.p_cov_inv_diag
#
#         # compute cross entropy term
#         tmp1 = 0
#         if p_cov[k].__class__.__name__ == 'ndarray':
#             mat_tmp = p_cov_inv[k] - p_cov_inv_diag[k] * s.eye(self.N)
#             tmp1 += QE[:, k].transpose().dot(mat_tmp).dot(QE[:, k])
#         tmp1 += p_cov_inv_diag[k].dot(QE2[:, k])
#         tmp1 = -.5 * tmp1
#         # tmp1 = 0.5*QE2 - PE*QE + 0.5*PE2
#         # tmp1 = -(tmp1 * Alpha['E']).sum()
#
#         # compute term from the precision factor in front of the Gaussian
#         tmp2 = 0  # constant here
#         # if self.n_iter> 4:
#         if p_cov[k].__class__.__name__ == 'dia_matrix':
#             tmp2 += np.sum(np.log(p_cov[k].data.flatten()))
#         elif p_cov[k].__class__.__name__ == 'ndarray':
#             tmp2 += np.linalg.slogdet(p_cov[k])[1]
#         else:
#             print("Not implemented yet")
#             exit()
#         tmp2 *= (-.5)
#         # tmp2 = 0.5*Alpha["lnE"].sum()
#
#         lb_p = tmp1 + tmp2
#         # lb_q = -(s.log(Qvar).sum() + self.N*self.dim[1])/2. # I THINK THIS IS WRONG BECAUSE SELF.DIM[1] ICNLUDES COVARIATES
#         lb_q = -.5 * s.log(Qvar[:, k]).sum()
#
#         return lb_p - lb_q
#
#     def calculateELBO(self):
#         if not ("AlphaZ" in self.markov_blanket):
#             latent_variables = self.getLvIndex()
#             elbo = 0
#             for k in latent_variables:
#                 elbo += self.calculateELBO_k(k)
#
#             elbo += .5 * self.N * len(latent_variables)
#
#             return elbo
#
#         else:
#             # Collect parameters and expectations of current node
#             Qpar, Qexp = self.Q.getParameters(), self.Q.getExpectations()
#             Qmean, Qvar = Qpar['mean'], Qpar['var']
#             QE, QE2 = Qexp['E'], Qexp['E2']
#
#             if "MuZ" in self.markov_blanket:
#                 PE, PE2 = self.markov_blanket['MuZ'].getExpectations()['E'], \
#                           self.markov_blanket['MuZ'].getExpectations()['E2']
#             else:
#                 PE, PE2 = self.P.getParameters()["mean"], s.zeros((self.N, self.dim[1]))
#
#             Alpha = self.markov_blanket[
#                 'AlphaZ'].getExpectations(expand=True).copy()  # Notice that this Alpha is the ARD prior on Z, not on W.
#
#             # This ELBO term contains only cross entropy between Q and P,and entropy of Q. So the covariates should not intervene at all
#             latent_variables = self.getLvIndex()
#             Alpha["E"], Alpha["lnE"] = Alpha["E"][:, latent_variables], Alpha["lnE"][:, latent_variables]
#             Qmean, Qvar = Qmean[:, latent_variables], Qvar[:, latent_variables]
#             PE, PE2 = PE[:, latent_variables], PE2[:, latent_variables]
#             QE, QE2 = QE[:, latent_variables], QE2[:, latent_variables]
#
#             # compute term from the exponential in the Gaussian
#             tmp1 = 0.5 * QE2 - PE * QE + 0.5 * PE2
#             tmp1 = -(tmp1 * Alpha['E']).sum()
#
#             # compute term from the precision factor in front of the Gaussian
#             tmp2 = 0.5 * Alpha["lnE"].sum()
#
#             lb_p = tmp1 + tmp2
#             # lb_q = -(s.log(Qvar).sum() + self.N*self.dim[1])/2. # I THINK THIS IS WRONG BECAUSE SELF.DIM[1] ICNLUDES COVARIATES
#             lb_q = -(s.log(Qvar).sum() + self.N * len(latent_variables)) / 2.
#
#             return lb_p - lb_q
#
#     def sample(self, dist='P'):
#         if "MuZ" in self.markov_blanket:
#             p_mean = self.markov_blanket['MuZ'].sample()
#         else:
#             p_mean = self.P.params['mean']
#
#         if "AlphaZ" in self.markov_blanket:
#             alpha = self.markov_blanket['AlphaZ'].sample()
#             p_var = s.square(1. / alpha)
#             #p_cov = s.diag(p_var)
#             p_cov = [p_var[k] * np.eye(self.N) for k in range(self.K)]
#         else:
#             if "SigmaZ" in self.markov_blanket:
#                 p_cov = self.markov_blanket['SigmaZ'].sample()
#             else:
#                 p_cov = self.P.params['cov']
#
#         # simulating
#
#         samp_tmp = []
#         for i in range(self.dim[1]):
#             # samp_tmp.append(s.random.multivariate_normal(p_mean[:, i], p_cov[i])) #does it yield the correct result for sparse input matrices ?
#             if p_cov[i].__class__.__name__ == 'dia_matrix':
#                 samp_tmp.append(s.random.multivariate_normal(p_mean[:, i], p_cov[i].toarray()))  # inefficient
#             elif p_cov[i].__class__.__name__ == 'ndarray':
#                 samp_tmp.append(s.random.multivariate_normal(p_mean[:, i], p_cov[i]))
#             else:
#                 print("Not implemented yet")
#                 exit()
#
#         # self.samp = s.array([tmp/tmp.std() for tmp in samp_tmp]).transpose()
#
#         self.samp = s.array([tmp - tmp.mean() for tmp in samp_tmp]).transpose()
#         #self.samp = np.array(samp_tmp).T
#
#         return self.samp
