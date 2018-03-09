
from __future__ import division
import numpy.ma as ma
import numpy as np
import scipy as s

# Import manually defined functions
from .variational_nodes import BernoulliGaussian_Unobserved_Variational_Node

class SW_Node(BernoulliGaussian_Unobserved_Variational_Node):
    # TOO MANY ARGUMENTS, SHOULD WE USE **KWARGS AND *KARGS ONLY?
    # def __init__(self, dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0=None, qEW_S1=None, qES=None):
    def __init__(self, dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0=None, qEW_S1=None, qES=None):
        super(SW_Node,self).__init__(dim, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0, qEW_S1, qES)
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
        Qmean_S1, Qvar_S1, Qtheta = Q['mean_S1'], Q['var_S1'], Q['theta']

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
        self.Q.setParameters(mean_S0=s.zeros((self.D,self.dim[1])), var_S0=s.repeat(1./alpha[None,:],self.D,0), mean_S1=Qmean_S1, var_S1=Qvar_S1, theta=Qtheta )

    def calculateELBO(self):

        # Collect parameters and expectations
        Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()
        S,WW = Qexp["ES"], Qexp["EWW"]
        Qvar = Qpar['var_S1']
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
