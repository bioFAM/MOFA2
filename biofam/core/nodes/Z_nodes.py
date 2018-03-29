from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
from copy import deepcopy

# Import manually defined functions
from .variational_nodes import UnivariateGaussian_Unobserved_Variational_Node
from .variational_nodes import BernoulliGaussian_Unobserved_Variational_Node


class Z_Node(UnivariateGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, pmean, pvar, qmean, qvar, qE=None, qE2=None, idx_covariates=None):
        super().__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)
        self.precompute()

        # Define indices for covariates
        if idx_covariates is not None:
            self.covariates[idx_covariates] = True

    def precompute(self):
        # Precompute terms to speed up computation
        self.N = self.dim[0]
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

        # Collect expectations from the markov blanket
        Y = deepcopy(self.markov_blanket["Y"].getExpectation())
        SWtmp = self.markov_blanket["SW"].getExpectations()
        tau = deepcopy(self.markov_blanket["Tau"].getExpectation())
        latent_variables = self.getLvIndex() # excluding covariates from the list of latent variables
        mask = [ma.getmask(Y[m]) for m in range(len(Y))]

        # Collect parameters from the prior or expectations from the markov blanket
        if "Mu" in self.markov_blanket:
            Mu = self.markov_blanket['Mu'].getExpectation()
        else:
            Mu = self.P.getParameters()["mean"]

        if "Alpha" in self.markov_blanket:
            Alpha = self.markov_blanket['Alpha'].getExpectation()
            Alpha = s.repeat(Alpha[None,:], self.N, axis=0)
        else:
            Alpha = 1./self.P.getParameters()["var"]

        # Check dimensionality of Tau and expand if necessary (for Jaakola's bound only)
        for m in range(len(Y)):
            # DEPRECATED: tau is expanded inside the node
            # if tau[m].shape != Y[m].shape:
            #     tau[m] = s.repeat(tau[m].copy()[None,:], self.N, axis=0)

            # Mask tau
            # tau[m] = ma.masked_where(ma.getmask(Y[m]), tau[m]) # important to keep this out of the loop to mask non-gaussian tau
            tau[m][mask[m]] = 0.
            # Mask Y
            Y[m] = Y[m].data
            Y[m][mask[m]] = 0.

        # Collect parameters from the P and Q distributions of this node
        Q = self.Q.getParameters().copy()
        Qmean, Qvar = Q['mean'], Q['var']

        M = len(Y)
        for k in latent_variables:
            foo = s.zeros((self.N,))
            bar = s.zeros((self.N,))
            for m in range(M):
                foo += np.dot(tau[m],SWtmp[m]["E2"][:,k])
                bar += np.dot(tau[m]*(Y[m] - s.dot( Qmean[:,s.arange(self.dim[1])!=k] , SWtmp[m]["E"][:,s.arange(self.dim[1])!=k].T )), SWtmp[m]["E"][:,k])
            Qvar[:,k] = 1./(Alpha[:,k]+foo)
            Qmean[:,k] = Qvar[:,k] * (  Alpha[:,k]*Mu[:,k] + bar )

        # Save updated parameters of the Q distribution
        self.Q.setParameters(mean=Qmean, var=Qvar)

    def calculateELBO(self):
        # Collect parameters and expectations of current node
        Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()
        Qmean, Qvar = Qpar['mean'], Qpar['var']
        QE, QE2 = Qexp['E'],Qexp['E2']

        if "Mu" in self.markov_blanket:
            PE, PE2 = self.markov_blanket['Mu'].getExpectations()['E'], self.markov_blanket['Mu'].getExpectations()['E2']
        else:
            PE, PE2 = self.P.getParameters()["mean"], s.zeros((self.N,self.dim[1]))

        if "Alpha" in self.markov_blanket:
            Alpha = self.markov_blanket['Alpha'].getExpectations().copy() # Notice that this Alpha is the ARD prior on Z, not on W.
            Alpha["E"] = s.repeat(Alpha["E"][None,:], self.N, axis=0)
            Alpha["lnE"] = s.repeat(Alpha["lnE"][None,:], self.N, axis=0)
        else:
            Alpha = { 'E':1./self.P.getParameters()["var"], 'lnE':s.log(1./self.P.getParameters()["var"]) }

        # This ELBO term contains only cross entropy between Q and P,and entropy of Q. So the covariates should not intervene at all
        latent_variables = self.getLvIndex()
        Alpha["E"], Alpha["lnE"] = Alpha["E"][:,latent_variables], Alpha["lnE"][:,latent_variables]
        Qmean, Qvar = Qmean[:, latent_variables], Qvar[:, latent_variables]
        PE, PE2 = PE[:, latent_variables], PE2[:, latent_variables]
        QE, QE2 = QE[:, latent_variables], QE2[:, latent_variables]

        # compute term from the exponential in the Gaussian
        tmp1 = 0.5*QE2 - PE*QE + 0.5*PE2
        tmp1 = -(tmp1 * Alpha['E']).sum()

        # compute term from the precision factor in front of the Gaussian
        tmp2 = 0.5*Alpha["lnE"].sum()

        lb_p = tmp1 + tmp2
        # lb_q = -(s.log(Qvar).sum() + self.N*self.dim[1])/2. # I THINK THIS IS WRONG BECAUSE SELF.DIM[1] ICNLUDES COVARIATES
        lb_q = -(s.log(Qvar).sum() + self.N*len(latent_variables))/2.

        return lb_p-lb_q

    def sample(self, dist='P'):
        if "Mu" in self.markov_blanket:
            print('Not implemented')
            exit(1)
        if "Alpha" in self.markov_blanket:
            print("Not implemented")  # TODO where is the AlphaZ node ?
            exit(1)

        p_mean = self.P.params['mean']
        p_var = self.P.params['var']

        # simulating and handling covariates
        self.samp = s.random.normal(p_mean, np.sqrt(p_var))
        self.samp[:, self.covariates] = self.getExpectation()[:, self.covariates]

        return self.samp


class SZ_Node(BernoulliGaussian_Unobserved_Variational_Node):
    # TOO MANY ARGUMENTS, SHOULD WE USE **KWARGS AND *KARGS ONLY?
    def __init__(self, dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEZ_S0=None, qEZ_S1=None, qES=None, idx_covariates=None):
        super().__init__(dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEZ_S0, qEZ_S1, qES)
        self.precompute()

        # Define indices for covariates
        if idx_covariates is not None:
            self.covariates[idx_covariates] = True

    def precompute(self):
        self.N = self.dim[0]
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
        Wtmp = [Wtmp_m.copy() for Wtmp_m in self.markov_blanket["W"].getExpectations()]
        W = [Wtmp_m["E"] for Wtmp_m in Wtmp]
        WW = [Wtmp_m["E2"] for Wtmp_m in Wtmp]

        tau = [tau_m.copy() for tau_m in self.markov_blanket["Tau"].getExpectation()]
        Y = [Y_m.copy() for Y_m in self.markov_blanket["Y"].getExpectation()]
        alpha = self.markov_blanket["AlphaZ"].getExpectation().copy()
        thetatmp = self.markov_blanket['ThetaZ'].getExpectations().copy()
        theta_lnE, theta_lnEInv  = thetatmp['lnE'], thetatmp['lnEInv']
        mask = ma.getmask(Y)

        M = len(Y)  #number of views

        # Collect parameters and expectations from P and Q distributions of this node
        SZ = self.Q.getExpectations()["E"]
        Q = self.Q.getParameters()
        Qmean_S1, Qvar_S1, Qtheta = Q['mean_B1'], Q['var_B1'], Q['theta']

        # Check dimensions of theta and and expand if necessary
        if theta_lnE.shape != Qmean_S1.shape:
            theta_lnE = s.repeat(theta_lnE[None,:], Qmean_S1.shape[0],0)
        if theta_lnEInv.shape != Qmean_S1.shape:
            theta_lnEInv = s.repeat(theta_lnEInv[None,:], Qmean_S1.shape[0],0)

        # DEPRECATED: tau is expanded inside the node
        # Check dimensions of Tau and and expand if necessary
        # M = len(Y) #number of views
        # for m in range(M):
        #     if tau[m].shape != Y[m].shape:
        #         tau[m] = s.repeat(tau[m][None,:], Y[m].shape[0], axis=0)
        # tau = ma.masked_where(ma.getmask(Y), tau)

        # Check dimensions of alpha and and expand if necessary
        if alpha.shape[0] == 1:
            alpha = s.repeat(alpha[:], self.dim[1], axis=0)

        # Mask matrices (excluding covariates from the list of latent variables)
        latent_variables = self.getLvIndex()
        mask = [ma.getmask(Y[m]) for m in range(len(Y))]
        for m in range(M):
            Y[m] = Y[m].data
            Y[m][mask[m]] = 0.
            tau[m][mask[m]] = 0.

        # TODO : remove loops working with tensors ?
        # Update each latent variable in turn
        for k in range(self.dim[1]):

            # Calculate intermediate stept
            term1 = (theta_lnE - theta_lnEInv)[:,k]
            term2 = 0.5*s.log(alpha[k])
            # term3 = 0.5*s.log(ma.dot(WW[:,k],tau) + alpha[k])
            term3 = 0.5*s.log(s.sum([s.dot(tau[m], WW[m][:,k]) for m in range(M)], axis=0) + alpha[k]) # good to modify

            # term4_tmp1 = ma.dot((tau*Y).T,W[:,k]).data
            term4_tmp1 = s.sum([s.dot(tau[m]*Y[m], W[m][:,k]) for m in range(M)], axis=0) # good to modify
            # term4_tmp2 = ( tau * s.dot((W[:,k]*W[:,s.arange(self.dim[1])!=k].T).T, SZ[:,s.arange(self.dim[1])!=k].T) ).sum(axis=0)

            # term4_tmp2 = s.sum([(tau[m] * s.dot(SZ[:,s.arange(self.dim[1])!=k], W[m][:,k]*W[m][:,s.arange(self.dim[1])!=k].T) ).sum(axis=1) for m in range(M)], axis=0)
            term4_tmp2 = s.sum([(tau[m] * s.dot(SZ[:, s.arange(self.dim[1]) != k], (W[m][:, k] * W[m][:, s.arange(self.dim[1]) != k].T))).sum(axis=1) for m in range(M)], axis=0)  # good to modify

            term4_tmp3 = s.sum([s.dot(tau[m], WW[m][:,k]) for m in range(M)], axis=0) + alpha[k]

            term4 = 0.5 * s.divide(s.square(term4_tmp1-term4_tmp2), term4_tmp3) # good to modify, awsnt checked numerically

            # Update S
            # NOTE there could be some precision issues in S --> loads of 1s in result
            Qtheta[:,k] = 1./(1.+s.exp(-(term1+term2-term3+term4)))

            # Update Z
            Qvar_S1[:,k] = 1./term4_tmp3
            Qmean_S1[:,k] = Qvar_S1[:,k] * (term4_tmp1-term4_tmp2)

            # Update Expectations for the next iteration
            SZ[:,k] = Qtheta[:,k] * Qmean_S1[:,k]

        # Save updated parameters of the Q distribution
        self.Q.setParameters(mean_B0=s.zeros((self.N,self.dim[1])), var_B0=s.repeat(1./alpha[None,:],self.N,0), mean_B1=Qmean_S1, var_B1=Qvar_S1, theta=Qtheta )

    def calculateELBO(self):

        # Collect parameters and expectations
        Qpar, Qexp = self.Q.getParameters(), self.Q.getExpectations()
        S, ZZ = Qexp["EB"], Qexp["ENN"]
        Qvar = Qpar['var_B1']
        theta = self.markov_blanket['ThetaZ'].getExpectations()

        # Get ARD sparsity or prior variance
        if "AlphaZ" in self.markov_blanket:
            alpha = self.markov_blanket['AlphaZ'].getExpectations().copy()
            if alpha["E"].shape[0] == 1:
                alpha["E"] = s.repeat(alpha["E"][:], self.dim[1], axis=0)
                alpha["lnE"] = s.repeat(alpha["lnE"][:], self.dim[1], axis=0)
        else:
            print("Not implemented")
            exit()

        # This ELBO term contains only cross entropy between Q and P, and entropy of Q. So the covariates should not intervene at all
        latent_variables = self.getLvIndex()
        alpha["E"], alpha["lnE"] = alpha["E"][latent_variables], alpha["lnE"][latent_variables]
        S, ZZ = S[:,latent_variables], ZZ[:,latent_variables]
        Qvar = Qvar[:,latent_variables]
        theta['lnE'] = theta['lnE'][latent_variables]
        theta['lnEInv'] = theta['lnEInv'][latent_variables]

        # Calculate ELBO for Z
        lb_pz = (self.N*alpha["lnE"].sum() - s.sum(alpha["E"]*ZZ))/2.
        lb_qz = -0.5*self.dim[1]*self.N - 0.5*(S*s.log(Qvar) + (1.-S)*s.log(1./alpha["E"])).sum() # IS THE FIRST CONSTANT TERM CORRECT???
        lb_z = lb_pz - lb_qz

        # Calculate ELBO for S
        lb_ps = S*theta['lnE'] + (1.-S)*theta['lnEInv']
        lb_qs = S*s.log(S) + (1.-S)*s.log(1.-S)
        lb_ps[s.isnan(lb_ps)] = 0.
        lb_qs[s.isnan(lb_qs)] = 0.
        lb_s = s.sum(lb_ps) - s.sum(lb_qs)

        return lb_z + lb_s

    def sample(self, dist='P'):
        # get necessary parameters
        mu_z_hat = self.P.getParameters()['mean_S1']
        mu_z_hat = s.ones(self.dim) * mu_z_hat

        theta = self.markov_blanket['ThetaZ'].sample()
        if theta.shape != mu_z_hat.shape:
            theta = s.repeat(theta[None,:],mu_z_hat.shape[0],0)

        alpha = self.markov_blanket['AlphaZ'].sample()
        if alpha.shape[0] == 1:
            alpha = s.repeat(alpha[:], self.dim[1], axis=0)
        if alpha.shape != mu_z_hat.shape:
            alpha = s.repeat(alpha[None,:], self.N, axis=0)

        # simulate
        bernoulli_t = s.random.binomial(1, theta)
        z_hat = s.random.normal(mu_z_hat, np.sqrt(1./alpha))

        self.samp = bernoulli_t * z_hat
        self.samp[:, self.covariates] = self.getExpectation()[:, self.covariates]

        return self.samp

# class MuZ_Node(UnivariateGaussian_Unobserved_Variational_Node):
#     """ """
#     def __init__(self, pmean, pvar, qmean, qvar, clusters, n_Z, cluster_dic=None, qE=None, qE2=None):
#         # compute dim from numbers of clusters (n_clusters * Z)
#         self.clusters = clusters
#         self.N = len(self.clusters)
#         self.n_clusters = len(np.unique(clusters))
#         dim = (self.n_clusters, n_Z)
#         self.factors_axis = 1
#         super().__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)

#     def getExpectations(self):
#         # reshape the values to N_samples * N_factors and return
#         QExp = self.Q.getExpectations()
#         expanded_expectation = QExp['E'][self.clusters, :]
#         expanded_E2 = QExp['E2'][self.clusters, :]
#         # do we need to expand the variance as well -> not used I think
#         return {'E': expanded_expectation , 'E2': expanded_E2}

#     def updateParameters(self):
#         Ppar = self.P.getParameters()
#         Z = self.markov_blanket['Z'].Q.getExpectation()

#         if "Alpha" in self.markov_blanket:
#             Alpha = self.markov_blanket['Alpha'].getExpectation().copy() # Notice that this Alpha is the ARD prior on Z, not on W.
#             Alpha = s.repeat(Alpha[None,:], self.N, axis=0)
#         else:
#             Alpha = 1./self.markov_blanket['Z'].P.getParameters()["var"]

#         Qmean, Qvar = self.Q.getParameters()['mean'], self.Q.getParameters()['var']
#         ZTauMean = Z * Alpha

#         # TODO merge two loops when sure it's clean
#         # update of the variance
#         for c in range(self.n_clusters):
#             mask = (self.clusters == c)
#             tmp = (Alpha[mask, :]).sum(axis=0)
#             Qvar[c,:] = tmp
#         Qvar += 1./Ppar['var']
#         Qvar = 1./Qvar

#         # update of the mean
#         for c in range(self.n_clusters):
#             mask = (self.clusters == c)
#             tmp = (ZTauMean[mask, :]).sum(axis=0)
#             Qmean[c,:] = tmp
#         Qmean = Qmean + Ppar['mean']/Ppar['var']
#         Qmean *= Qvar

#         self.Q.setParameters(mean=Qmean, var=Qvar)

#     def calculateELBO(self):
#         PParam = self.P.getParameters()
#         PVar, Pmean = PParam['var'], PParam['mean']

#         QExp = self.Q.getExpectations()
#         QE2, QE = QExp['E2'], QExp['E']

#         Qvar = self.Q.getParameters()['var']

#         # Cluster terms corersponding to covariates should not intervene
#         # filtering the covariates out
#         latent_variables = self.markov_blanket['Z'].getLvIndex()
#         PVar, Pmean = PVar[:, latent_variables], Pmean[:, latent_variables]
#         QE2, QE = QE2[:, latent_variables], QE[:, latent_variables]
#         Qvar = Qvar[:, latent_variables]

#         # minus cross entropy
#         tmp = -(0.5 * s.log(PVar)).sum()
#         tmp2 = - ((0.5/PVar) * (QE2 - 2.*QE*Pmean + Pmean**2.)).sum()

#         # entropy of Q
#         tmp3 = 0.5 * (s.log(Qvar)).sum()
#         tmp3 += 0.5 * self.dim[0] * len(latent_variables)

#         return tmp + tmp2 + tmp3
