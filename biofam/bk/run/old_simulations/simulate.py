"""
Module to simulate test data

To-do:
- Currently non-gaussian likelihoods have no noise
- Maybe integrate this into the corresponding Node classes with a Simulate() method?
- Fix binomial
"""

from __future__ import division
import scipy as s
import pandas as pd
import warnings
from scipy.stats import bernoulli, norm, gamma, uniform, poisson, binom
from random import sample

from biofam.core.utils import sigmoid


class Simulate(object):
    def __init__(self, M, N, D, K):
        """General method to Simulate from the generative model

        PARAMETERS
        ----------
        M (int): number of views
        N (int): number of samples
        D (list/tuple of length M): dimensionality of each view
        K (int): number of latent variables
        """

        # Sanity checks
        assert len(D) == M
        assert K < min(D)
        assert K < N

        self.M = M
        self.N = N
        self.K = K
        self.D = D

    def initAlpha(self):
        """ Initialisation of ARD on the weights"""
        alpha = [ s.zeros(self.K,) for m in range(self.M) ]
        for m in range(self.M):
            tmp = bernoulli.rvs(p=0.5, size=self.K)
            tmp[tmp==1] = 1.
            tmp[tmp==0] = 1E5
            alpha[m] = tmp
        return alpha

    def initW_ard(self, alpha=None):
        """ Initialisation of weights in automatic relevance determination prior"""
        if alpha is None:
            alpha = self.initAlpha()
        W = [ s.zeros((self.D[m],self.K)) for m in range(self.M) ]
        for m in range(self.M):
            for k in range(self.K):
                W[m][:,k] = norm.rvs(loc=0, scale=1/s.sqrt(alpha[m][k]), size=self.D[m])
        return W,alpha

    def initW_spikeslab(self, theta, alpha=None):
        """ Initialisation of weights in spike and slab prior"""

        # Simualte ARD precision
        if alpha is None:
            alpha = self.initAlpha()
        else:
            assert not any([0 in a for a in alpha]), 'alpha cannot be zero'

        # Simulate bernoulli variable S
        S = [ s.zeros((self.D[m],self.K)) for m in range(self.M) ]
        for m in range(self.M):

            # Completely vectorised, not sure if it works
            # S[m] = bernoulli.rvs(p=theta[m].flatten(), size=self.D[m]*self.K).reshape((self.D[m],self.K))

            # Partially vectorised
            for k in range(self.K):
                S[m][:,k] = bernoulli.rvs(p=theta[m][:,k], size=self.D[m])

            # Unvectorised
            # for d in range(self.D[m]):
            #     for k in range(self.K):
            #         S[m][d,k] = bernoulli.rvs(p=theta[m][d,k], size=1)


        # Simulate gaussian weights W
        W_hat = [ s.empty((self.D[m],self.K)) for m in range(self.M) ]
        W = [ s.empty((self.D[m],self.K)) for m in range(self.M) ]
        for m in range(self.M):
            for k in range(self.K):
                W_hat[m][:,k] = norm.rvs(loc=0, scale=s.sqrt(1./alpha[m][k]), size=self.D[m])
            W[m] = W_hat[m] * S[m]

        return S, W, W_hat, alpha

    def initZ(self):
        """ Initialisation of latent variables"""
        Z = s.empty((self.N,self.K))
        for n in range(self.N):
            for k in range(self.K):
                Z[n,k] = norm.rvs(loc=0, scale=1, size=1)
        return Z

    def initTau(self):
        """ Initialisation of noise precision"""
        return [ uniform.rvs(loc=1,scale=3,size=self.D[m]) for m in range(self.M) ]

    def initMu(self):
        """ Initialisation of hyperprior means of latent variables"""
        # Means are initialised to zero by default
        return [ s.zeros(self.D[m]) for m in range(self.M) ]

    def generateData(self, W, Z, Tau, Mu, likelihood, missingness=0.0, missing_view=False):
        """ Initialisation of observations

        PARAMETERS
        ----------
        W (list of length M where each element is a np array with shape (Dm,K)): weights
        Z (np array with shape (N,K): latent variables
        Tau (list of length M where each element is a np array with shape (Dm,)): precision of the normally-distributed noise
        Mu (list of length M where each element is a np array with shape (Dm,)): feature-wise means
        likelihood (str): type of likelihood
        missingness (float): percentage of missing values
        """

        Y = [ s.zeros((self.N,self.D[m])) for m in range(self.M) ]

        if likelihood == "gaussian":
            # Vectorised
            for m in range(self.M):
                Y[m] = s.dot(Z,W[m].T) + Mu[m] + norm.rvs(loc=0, scale=1/s.sqrt(Tau[m]), size=[self.N, self.D[m]])
            # Non-vectorised, slow
            # for m in range(self.M):
                # for n in range(self.N):
                    # for d in range(self.D[m]):
                        # Y[m][n,d] = s.dot(Z[n,:],W[m][d,:].T) + Mu[m][d] + norm.rvs(loc=0,scale=1/s.sqrt(Tau[m][d]))

        elif likelihood == "warp":
            for m in range(self.M):
                Y[m] = s.exp(s.dot(Z,W[m].T) + Mu[m] + norm.rvs(loc=0, scale=1/s.sqrt(Tau[m]), size=[self.N, self.D[m]]))


        # Sample observations using a poisson likelihood
        elif likelihood == "poisson":
            # Slow way
            # for m in range(self.M):
            #     for n in range(self.N):
            #         for d in range(self.D[m]):
            #             f = s.dot(Z[n,:],W[m][d,:].T)
            #             # f = s.dot(Z[n,:],W[m][d,:].T) + norm.rvs(loc=0,scale=s.sqrt(1/Tau[m][d]))
            #             rate = s.log(1+s.exp(f))
            #             # Sample from the Poisson distribution
            #             # Y[m][n,d] = poisson.rvs(rate)
            #             # Use the more likely values
            #             Y[m][n,d] = s.special.round(rate)

            # Fast way
            for m in range(self.M):
                F = s.dot(Z,W[m].T)
                # F = s.dot(Z,W[m].T) + norm.rvs(loc=0,scale=s.sqrt(1/Tau[m]))
                rate = s.log(1+s.exp(F))
                # Sample from the Poisson distribution
                # MAYBE THIS REQUIRES RESHAPING
                # Y[m] = poisson.rvs(rate)
                # Use the more likely values
                Y[m] = s.special.round(rate)

        # Sample observations using a bernoulli likelihood
        elif likelihood == "bernoulli":
            for m in range(self.M):
                # for n in range(self.N):
                    # for d in range(self.D[m]):
                        # Without noise
                        # f = sigmoid( s.dot(Z[n,:],W[m][d,:].T) )
                        # With noise, problem: it shifts the sigmoid...
                        # f = sigmoid( s.dot(Z[n,:],W[m][d,:].T) + norm.rvs(loc=0,scale=1/s.sqrt(Tau[m][d])) )

                        # Sample from the Bernoulli distributionn
                        # Y[m][n,d] = bernoulli.rvs(f)
                        # Use the more likely state
                        # Y[m][n,d] = s.special.round(f)
                f = sigmoid( s.dot(Z,W[m].T) )
                Y[m] = s.special.round(f)

        # Introduce missing values into the data
        if missingness > 0.0:
            for m in range(self.M):
                nas = s.random.choice(range(self.N*self.D[m]), size=int(missingness*self.N*self.D[m]), replace=False)
                tmp = Y[m].flatten()
                tmp[nas] = s.nan
                Y[m] = tmp.reshape((self.N,self.D[m]))
        if missing_view > 0.0:   # percentage of samples missing a view
            # select samples missing one view
            n_missing = s.random.choice(range(self.N), int(missing_view * self.N), replace=False)
            Y[0][n_missing,:] = s.nan

        # Convert data to pandas data frame
        for m in range(self.M):
            Y[m] = pd.DataFrame(data=Y[m])

        return Y
