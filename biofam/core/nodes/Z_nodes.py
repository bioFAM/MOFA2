from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s
from copy import deepcopy

# Import manually defined functions
from .variational_nodes import UnivariateGaussian_Unobserved_Variational_Node


class Z_Node(UnivariateGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, pmean, pvar, qmean, qvar, qE=None, qE2=None, idx_covariates=None):
        super(Z_Node,self).__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)
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
            if tau[m].shape != Y[m].shape:
                tau[m] = s.repeat(tau[m].copy()[None,:], self.N, axis=0)
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
                foo += np.dot(tau[m],SWtmp[m]["ESWW"][:,k])
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


class MuZ_Node(UnivariateGaussian_Unobserved_Variational_Node):
    """ """
    def __init__(self, pmean, pvar, qmean, qvar, clusters, n_Z, cluster_dic=None, qE=None, qE2=None):
        # compute dim from numbers of clusters (n_clusters * Z)
        self.clusters = clusters
        self.N = len(self.clusters)
        self.n_clusters = len(np.unique(clusters))
        dim = (self.n_clusters, n_Z)
        self.factors_axis = 1
        super(Cluster_Node_Gaussian, self).__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)

    def getExpectations(self):
        # reshape the values to N_samples * N_factors and return
        QExp = self.Q.getExpectations()
        expanded_expectation = QExp['E'][self.clusters, :]
        expanded_E2 = QExp['E2'][self.clusters, :]
        # do we need to expand the variance as well -> not used I think
        return {'E': expanded_expectation , 'E2': expanded_E2}

    def updateParameters(self):
        Ppar = self.P.getParameters()
        Z = self.markov_blanket['Z'].Q.getExpectation()

        if "Alpha" in self.markov_blanket:
            Alpha = self.markov_blanket['Alpha'].getExpectation().copy() # Notice that this Alpha is the ARD prior on Z, not on W.
            Alpha = s.repeat(Alpha[None,:], self.N, axis=0)
        else:
            Alpha = 1./self.markov_blanket['Z'].P.getParameters()["var"]

        Qmean, Qvar = self.Q.getParameters()['mean'], self.Q.getParameters()['var']
        ZTauMean = Z * Alpha

        # TODO merge two loops when sure it's clean
        # update of the variance
        for c in range(self.n_clusters):
            mask = (self.clusters == c)
            tmp = (Alpha[mask, :]).sum(axis=0)
            Qvar[c,:] = tmp
        Qvar += 1./Ppar['var']
        Qvar = 1./Qvar

        # update of the mean
        for c in range(self.n_clusters):
            mask = (self.clusters == c)
            tmp = (ZTauMean[mask, :]).sum(axis=0)
            Qmean[c,:] = tmp
        Qmean = Qmean + Ppar['mean']/Ppar['var']
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
        tmp2 = - ((0.5/PVar) * (QE2 - 2.*QE*Pmean + Pmean**2.)).sum()

        # entropy of Q
        tmp3 = 0.5 * (s.log(Qvar)).sum()
        tmp3 += 0.5 * self.dim[0] * len(latent_variables)

        return tmp + tmp2 + tmp3
