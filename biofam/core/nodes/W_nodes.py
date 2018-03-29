from __future__ import division
from copy import deepcopy
import numpy.ma as ma
import numpy as np
import scipy as s

# Import manually defined functions
from .variational_nodes import BernoulliGaussian_Unobserved_Variational_Node
from. variational_nodes import UnivariateGaussian_Unobserved_Variational_Node

class W_Node(UnivariateGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, pmean, pvar, qmean, qvar, qE=None, qE2=None, idx_covariates=None):
        super().__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)
        self.precompute()

        # Define indices for covariates
        if idx_covariates is not None:
            self.covariates[idx_covariates] = True

    def precompute(self):
        # Precompute terms to speed up computation
        self.D = self.dim[0]
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
        SZtmp = self.markov_blanket["SZ"].getExpectations()
        tau = deepcopy(self.markov_blanket["Tau"].getExpectation())
        latent_variables = self.getLvIndex() # excluding covariates from the list of latent variables
        mask = ma.getmask(Y)

        # Collect parameters from the prior or expectations from the markov blanket
        if "MuZ" in self.markov_blanket:
            Mu = self.markov_blanket['MuW'].getExpectation()
        else:
            Mu = self.P.getParameters()["mean"]

        if "AlphaW" in self.markov_blanket:
            Alpha = self.markov_blanket['AlphaW'].getExpectation()
            Alpha = s.repeat(Alpha[None,:], self.D, axis=0)
        else:
            Alpha = 1./self.P.getParameters()["var"]

        # DEPRECATED: tau is expanded inside the node
        # Check dimensionality of Tau and expand if necessary (for Jaakola's bound only)
        # if tau.shape != Y.shape:
        #     tau = s.repeat(tau.copy()[None,:], np.shape(Y)[0], axis=0)
        # Mask tau
        tau[mask] = 0.
        # Mask Y
        Y = Y.data
        Y[mask] = 0.

        # Collect parameters from the P and Q distributions of this node
        Q = self.Q.getParameters().copy()
        Qmean, Qvar = Q['mean'], Q['var']

        for k in latent_variables:
            foo = s.zeros((self.D,))
            bar = s.zeros((self.D,))
            foo += np.dot(SZtmp["E2"][:,k], tau)
            bar += np.dot(SZtmp["E"][:,k], tau * (Y - s.dot(SZtmp["E"][:,s.arange(self.dim[1])!=k], Qmean[:,s.arange(self.dim[1])!=k].T )))
            Qvar[:,k] = 1./(Alpha[:,k] + foo)
            Qmean[:,k] = Qvar[:,k] * ( Alpha[:,k] * Mu[:,k] + bar )

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

        if "AlphaW" in self.markov_blanket:
            Alpha = self.markov_blanket['AlphaW'].getExpectations().copy() # Notice that this Alpha is the ARD prior on Z, not on W.
            Alpha["E"] = s.repeat(Alpha["E"][None,:], self.D, axis=0)
            Alpha["lnE"] = s.repeat(Alpha["lnE"][None,:], self.D, axis=0)
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
        # lb_q = -(s.log(Qvar).sum() + self.D*self.dim[1])/2. # I THINK THIS IS WRONG BECAUSE SELF.DIM[1] ICNLUDES COVARIATES
        lb_q = -(s.log(Qvar).sum() + self.D*len(latent_variables))/2.

        return lb_p - lb_q

    def sample(self, dist='P'):
        if "MuW" in self.markov_blanket:
            p_mean = self.markov_blanket['MuW'].sample()
        else:
            p_mean = self.P.params['mean']
        if "AlphaW" in self.markov_blanket:
            alpha = self.markov_blanket['AlphaW'].sample()
            p_var = s.square(1./alpha)
        else:
            p_var = self.P.params['var']

        # simulating and handling covariates
        self.samp = s.random.normal(p_mean, np.sqrt(p_var))
        self.samp[:, self.covariates] = self.getExpectation()[:, self.covariates]

        return self.samp

class SW_Node(BernoulliGaussian_Unobserved_Variational_Node):
    # TOO MANY ARGUMENTS, SHOULD WE USE **KWARGS AND *KARGS ONLY?
    # def __init__(self, dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0=None, qEW_S1=None, qES=None):
    def __init__(self, dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0=None, qEW_S1=None, qES=None):
        super().__init__(dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0, qEW_S1, qES)
        self.precompute()

    def precompute(self):
        self.D = self.dim[0]
        self.factors_axis = 1

    def updateParameters(self):
        # Collect expectations from other nodes
        Ztmp = self.markov_blanket["Z"].getExpectations()
        Z,ZZ = Ztmp["E"],Ztmp["E2"]
        tau = self.markov_blanket["Tau"].getExpectation().copy()
        Y = self.markov_blanket["Y"].getExpectation().copy()
        alpha = self.markov_blanket["Alpha"].getExpectation().copy()
        thetatmp = self.markov_blanket['Theta'].getExpectations()
        theta_lnE, theta_lnEInv  = thetatmp['lnE'], thetatmp['lnEInv']
        mask = ma.getmask(Y)

        # Collect parameters and expectations from P and Q distributions of this node
        SW = self.Q.getExpectations()["E"]
        Q = self.Q.getParameters()
        Qmean_S1, Qvar_S1, Qtheta = Q['mean_B1'], Q['var_B1'], Q['theta']

        # Check dimensions of Theta and and expand if necessary
        if theta_lnE.shape != Qmean_S1.shape:
            theta_lnE = s.repeat(theta_lnE[None,:],Qmean_S1.shape[0],0)
        if theta_lnEInv.shape != Qmean_S1.shape:
            theta_lnEInv = s.repeat(theta_lnEInv[None,:],Qmean_S1.shape[0],0)

        # DEPRECATED: tau is expanded inside the node
        # Check dimensions of Tau and and expand if necessary
        # if tau.shape != Y.shape:
        #     tau = s.repeat(tau[None,:], Y.shape[0], axis=0)
        
        # tau = ma.masked_where(ma.getmask(Y), tau)

        # Check dimensions of Alpha and and expand if necessary
        if alpha.shape[0] == 1:
            alpha = s.repeat(alpha[:], self.dim[1], axis=0)

        # Mask matrices
        Y = Y.data
        Y[mask] = 0.
        tau[mask] = 0.

        # Update each latent variable in turn
        for k in range(self.dim[1]):

            # Calculate intermediate steps
            term1 = (theta_lnE-theta_lnEInv)[:,k]
            term2 = 0.5*s.log(alpha[k])
            # term3 = 0.5*s.log(ma.dot(ZZ[:,k],tau) + alpha[k])
            term3 = 0.5*s.log(s.dot(ZZ[:,k],tau) + alpha[k]) # good to modify
            # term4_tmp1 = ma.dot((tau*Y).T,Z[:,k]).data
            term4_tmp1 = s.dot((tau*Y).T,Z[:,k]) # good to modify
            # term4_tmp2 = ( tau * s.dot((Z[:,k]*Z[:,s.arange(self.dim[1])!=k].T).T, SW[:,s.arange(self.dim[1])!=k].T) ).sum(axis=0)
            term4_tmp2 = ( tau * s.dot((Z[:,k]*Z[:,s.arange(self.dim[1])!=k].T).T, SW[:,s.arange(self.dim[1])!=k].T) ).sum(axis=0) # good to modify

            term4_tmp3 = ma.dot(ZZ[:,k].T,tau) + alpha[k]
            # term4_tmp3 = s.dot(ZZ[:,k].T,tau) + alpha[k] # good to modify (I REPLACE MA.DOT FOR S.DOT, IT SHOULD BE SAFE )

            # term4 = 0.5*s.divide((term4_tmp1-term4_tmp2)**2,term4_tmp3)
            term4 = 0.5*s.divide(s.square(term4_tmp1-term4_tmp2),term4_tmp3) # good to modify, awsnt checked numerically

            # Update S
            # NOTE there could be some precision issues in S --> loads of 1s in result
            Qtheta[:,k] = 1./(1.+s.exp(-(term1+term2-term3+term4)))

            # Update W
            Qvar_S1[:,k] = 1./term4_tmp3
            Qmean_S1[:,k] = Qvar_S1[:,k]*(term4_tmp1-term4_tmp2)

            # Update Expectations for the next iteration
            SW[:,k] = Qtheta[:,k] * Qmean_S1[:,k]

        # Save updated parameters of the Q distribution
        self.Q.setParameters(mean_B0=s.zeros((self.D,self.dim[1])), var_B0=s.repeat(1./alpha[None,:],self.D,0), mean_B1=Qmean_S1, var_B1=Qvar_S1, theta=Qtheta )

    def calculateELBO(self):

        # Collect parameters and expectations
        Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()
        S,WW = Qexp["EB"], Qexp["ENN"]
        Qvar = Qpar['var_B1']
        theta = self.markov_blanket['Theta'].getExpectations()

        # Get ARD sparsity or prior variance
        if "Alpha" in self.markov_blanket:
            alpha = self.markov_blanket['Alpha'].getExpectations().copy()
            if alpha["E"].shape[0] == 1:
                alpha["E"] = s.repeat(alpha["E"][:], self.dim[1], axis=0)
                alpha["lnE"] = s.repeat(alpha["lnE"][:], self.dim[1], axis=0)
        else:
            print("Not implemented")
            exit()

        # Calculate ELBO for W
        lb_pw = (self.D*alpha["lnE"].sum() - s.sum(alpha["E"]*WW))/2.
        lb_qw = -0.5*self.dim[1]*self.D - 0.5*(S*s.log(Qvar) + (1.-S)*s.log(1./alpha["E"])).sum() # IS THE FIRST CONSTANT TERM CORRECT???
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
        mu_w_hat = self.P.getParameters()['mean_S1']
        mu_w_hat = s.ones(self.dim) * mu_w_hat

        theta = self.markov_blanket['Theta'].sample()
        if theta.shape != mu_w_hat.shape:
            theta = s.repeat(theta[None,:],mu_w_hat.shape[0],0)

        alpha = self.markov_blanket['Alpha'].sample()
        if alpha.shape[0] == 1:
            alpha = s.repeat(alpha[:], self.dim[1], axis=0)
        if alpha.shape != mu_w_hat.shape:
            alpha = s.repeat(alpha[None,:], self.dim[0], axis=0)

        # simulate
        bernoulli_s = s.random.binomial(1, theta)
        w_hat = s.random.normal(mu_w_hat, np.sqrt(1./alpha))

        self.samp = bernoulli_s * w_hat
        return self.samp
