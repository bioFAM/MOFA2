from __future__ import division
from copy import deepcopy
import numpy.ma as ma
import numpy as np
import scipy as s
from copy import deepcopy
import math

# Import manually defined functions
from .variational_nodes import BernoulliGaussian_Unobserved_Variational_Node
from .variational_nodes import UnivariateGaussian_Unobserved_Variational_Node
from .variational_nodes import UnivariateGaussian_Unobserved_Variational_Node_with_MultivariateGaussian_Prior

# TODO make more general to handle both cases with and without SL in the Markov Blanket
class W_Node(UnivariateGaussian_Unobserved_Variational_Node_with_MultivariateGaussian_Prior):
    def __init__(self, dim, pmean, pcov, qmean, qvar, qE=None, qE2=None, idx_covariates=None,precompute_pcovinv=True):
        super().__init__(dim=dim, pmean=pmean, pcov=pcov, qmean=qmean, qvar=qvar, axis_cov=0, qE=qE, qE2=qE2)

        self.precompute(precompute_pcovinv=precompute_pcovinv)

        # Define indices for covariates
        if idx_covariates is not None:
            self.covariates[idx_covariates] = True

    def precompute(self,precompute_pcovinv=True):
        # Precompute terms to speed up computation
        self.D = self.dim[0]
        self.K = self.dim[1]
        self.covariates = np.zeros(self.dim[1], dtype=bool)
        self.factors_axis = 1

        if precompute_pcovinv:
            p_cov = self.P.params["cov"]

            self.p_cov_inv = []
            self.p_cov_inv_diag = []

            for k in range(self.K):
                if p_cov[k].__class__.__name__ == 'dia_matrix':
                    diag_inv = 1 / p_cov[k].data
                    diag_inv = diag_inv.flatten()
                    inv = np.diag(diag_inv)
                elif p_cov[k].__class__.__name__ == 'ndarray':
                    diagonal = np.diagonal(p_cov[k])
                    if np.all(np.diag(diagonal) == p_cov[k]):
                        diag_inv = 1. / diagonal
                        inv = np.diag(diag_inv)
                    else:
                        inv = np.linalg.inv(p_cov[k])
                        diag_inv = np.diagonal(inv)
                else:
                    #TODO : deal with sparse non diagonal input matrices as pcov
                    print("Not implemented yet")
                    exit()

                self.p_cov_inv.append(inv)
                self.p_cov_inv_diag.append(diag_inv)

        else:
            self.p_cov_inv = None
            self.p_cov_inv_diag = None


    def getLvIndex(self):
        # Method to return the index of the latent variables (without covariates)
        latent_variables = np.array(range(self.dim[1]))
        if any(self.covariates):
            # latent_variables = np.delete(latent_variables, latent_variables[self.covariates])
            latent_variables = latent_variables[~self.covariates]
        return latent_variables

    def removeFactors(self, idx):
        super().removeFactors(idx, axis=1)
        if self.p_cov_inv is not None:
            for i in idx :
                del self.p_cov_inv[i]
                del self.p_cov_inv_diag[i]
            #self.p_cov_inv = s.delete(self.p_cov_inv, axis=0, obj=idx)
            #self.p_cov_inv_diag = s.delete(self.p_cov_inv_diag, axis=0, obj=idx)

    def updateParameters(self, ix=None, ro=None):
        """
        Public function to update the nodes parameters
        Optional arguments for stochastic updates are:
            - ix: list of indices of the minibatch
            - ro: step size of the natural gradient ascent
        """
        if ix is None:
            ix = range(self.dim[0])

        ########################################################################
        # get Expectations which are necessarry for the update
        ########################################################################
        Y = self.markov_blanket["Y"].getExpectation()
        ZE = self.markov_blanket["Z"].getExpectations()['E']
        ZE2 = self.markov_blanket["Z"].getExpectations()['E2']
        tau = self.markov_blanket["Tau"].getExpectation()
        mask = ma.getmask(Y)
        N = Y.shape[0]

        # Collect parameters from the prior or expectations from the markov blanket
        if "MuW" in self.markov_blanket:
            Mu = self.markov_blanket['MuW'].getExpectation()
        else:
            Mu = self.P.getParameters()["mean"]

        if "AlphaW" in self.markov_blanket:
            Alpha = self.markov_blanket['AlphaW'].getExpectation(expand=True)
            p_cov_inv = None
            p_cov_inv_diag = None
        elif "SigmaAlphaW" in self.markov_blanket:
            if self.markov_blanket["SigmaAlphaW"].__class__.__name__=="AlphaW_Node_mk":
                Alpha = self.markov_blanket['SigmaAlphaW'].getExpectation(expand=True)
                p_cov_inv = None
                p_cov_inv_diag = None
            else:
                Sigma = self.markov_blanket['SigmaAlphaW'].getExpectations()
                p_cov_inv = Sigma['inv']
                p_cov_inv_diag = Sigma['inv_diag']
                Alpha = None
        else:
            Alpha = None
            p_cov_inv = self.p_cov_inv
            p_cov_inv_diag = self.p_cov_inv_diag

        # Collect parameters from the Q distributions of this node
        Q = self.Q.getParameters().copy()
        Qmean, Qvar = Q['mean'], Q['var']

        ########################################################################
        # subset matrices for stochastic inference
        ########################################################################
        Y = Y.data[ix,:].copy()
        mask = mask[ix,:].copy()
        tau = tau[ix,:].copy()
        ZE = ZE[ix,:].copy()
        ZE2 = ZE2[ix,:].copy()
        Z = {'E': ZE, 'E2': ZE2}

        ########################################################################
        # Masking
        ########################################################################
        tau[mask] = 0.
        Y[mask] = 0.

        ########################################################################
        # compute stochastic "anti-bias" coefficient
        ########################################################################
        coeff = float(N) / float(len(ix))

        ########################################################################
        # compute the update
        ########################################################################
        par_up = self._updateParameters(Y, Z, tau, Mu, Alpha, p_cov_inv, p_cov_inv_diag,
                               Qmean, Qvar, coeff)

        ########################################################################
        # Do the asignment
        ########################################################################
        if ro is not None: # TODO have a default ro of 1 instead ? whats the overhead cost ?
            par_up['Qmean'] = ro * par_up['Qmean'] + (1-ro) * self.Q.getParameters()['mean']
            par_up['Qvar'] = ro * par_up['Qvar'] + (1-ro) * self.Q.getParameters()['var']
        self.Q.setParameters(mean=par_up['Qmean'], var=par_up['Qvar'])


    # TODO: use coef where appropriate
    def _updateParameters(self, Y, Z, tau, Mu, Alpha, p_cov_inv, p_cov_inv_diag,
                           Qmean, Qvar, coeff):

        latent_variables = self.getLvIndex() # excluding covariates from the list of latent variables

        for k in latent_variables:
            foo = s.zeros((self.D,))
            bar = s.zeros((self.D,))

            foo += coeff * np.dot(Z["E2"][:,k],tau)
            bar += coeff * np.dot(Z["E"][:,k],tau*(Y - s.dot(Z["E"][:,s.arange(self.dim[1])!=k], Qmean[:,s.arange(self.dim[1])!=k].T )))

            b = ("SigmaAlphaW" in self.markov_blanket) and (
                    self.markov_blanket["SigmaAlphaW"].__class__.__name__ == "AlphaW_Node_mk")
            b = b or ("AlphaW" in self.markov_blanket)
            if b:
                Qvar[:,k] = 1./(Alpha[:,k]+foo)
                Qmean[:,k] = Qvar[:,k] * (bar + Alpha[:,k]*Mu[:,k])

            else:
                Qvar[:, k] = 1. / (foo + p_cov_inv_diag[k])

                if self.P.params["cov"][k].__class__.__name__ == 'dia_matrix':
                    Qmean[:, k] = Qvar[:, k] * bar
                else:
                    tmp = p_cov_inv[k] - p_cov_inv_diag[k] * s.eye(self.D)
                    for d in range(self.D):
                        Qmean[d, k] = Qvar[d, k] * (bar[d] + np.dot(tmp[d, :],Mu[:,k]-Qmean[:, k])) #-Qmean[:, k]))

        # Save updated parameters of the Q distribution
        return {'Qmean': Qmean, 'Qvar': Qvar}

    # TODO, problem here is that we need to make sure k is in the latent variables first
    def calculateELBO_k(self, k):
        '''Compute the ELBO for factor k in absence of Alpha node in the markov blanket of W'''
        Qpar, Qexp = self.Q.getParameters(), self.Q.getExpectations()
        Qmean, Qvar = Qpar['mean'], Qpar['var']
        QE, QE2 = Qexp['E'], Qexp['E2']

        if "SigmaAlphaW" in self.markov_blanket:
            Sigma = self.markov_blanket['SigmaAlphaW'].getExpectations()
            p_cov = Sigma['cov']
            p_cov_inv = Sigma['inv']
            p_cov_inv_diag = Sigma['inv_diag']
        else:
            p_cov = self.P.params['cov']
            p_cov_inv = self.p_cov_inv
            p_cov_inv_diag = self.p_cov_inv_diag

        # compute cross entropy term
        tmp1 = 0
        if p_cov[k].__class__.__name__ == 'ndarray':
            mat_tmp = p_cov_inv[k] - p_cov_inv_diag[k] * s.eye(self.D)
            tmp1 += QE[:, k].transpose().dot(mat_tmp).dot(QE[:, k])
        tmp1 += p_cov_inv_diag[k].dot(QE2[:, k])
        tmp1 = -.5 * tmp1
        # tmp1 = 0.5*QE2 - PE*QE + 0.5*PE2
        # tmp1 = -(tmp1 * Alpha['E']).sum()

        # compute term from the precision factor in front of the Gaussian
        tmp2 = 0  # constant here
        # if self.n_iter> 4:
        #     import pdb; pdb.set_trace()
        if p_cov[k].__class__.__name__ == 'dia_matrix':
            tmp2 += np.sum(np.log(p_cov[k].data.flatten()))
        elif p_cov[k].__class__.__name__ == 'ndarray':
            tmp2 += np.linalg.slogdet(p_cov[k])[1]
        else:
            print("Not implemented yet")
            exit()
        tmp2 *= (-.5)
        # tmp2 = 0.5*Alpha["lnE"].sum()

        lb_p = tmp1 + tmp2
        # lb_q = -(s.log(Qvar).sum() + self.N*self.dim[1])/2. # I THINK THIS IS WRONG BECAUSE SELF.DIM[1] ICNLUDES COVARIATES
        lb_q = -.5 * s.log(Qvar[:, k]).sum()
        # import pdb; pdb.set_trace()

        #print("d")

        return lb_p - lb_q

    def calculateELBO(self):
        b = ("SigmaAlphaW" in self.markov_blanket) and (
                self.markov_blanket["SigmaAlphaW"].__class__.__name__ == "AlphaW_Node_mk")
        b = b or ("AlphaW" in self.markov_blanket)

        if not(b):
            latent_variables = self.getLvIndex()
            elbo = 0
            for k in latent_variables:
                elbo += self.calculateELBO_k(k)

            elbo += .5 * self.D * len(latent_variables)

            return elbo

        else:
            # Collect parameters and expectations of current node
            Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()
            Qmean, Qvar = Qpar['mean'], Qpar['var']
            QE, QE2 = Qexp['E'],Qexp['E2']

            if "MuW" in self.markov_blanket:
                PE, PE2 = self.markov_blanket['MuW'].getExpectations()['E'], self.markov_blanket['MuW'].getExpectations()['E2']
            else:
                PE, PE2 = self.P.getParameters()["mean"], s.zeros((self.D,self.dim[1]))

            if "SigmaAlphaW" in self.markov_blanket:
                name_alpha = "SigmaAlphaW"
            else:
                name_alpha = "AlphaW"
            Alpha = self.markov_blanket[name_alpha].getExpectations(expand=True).copy() # Notice that this Alpha is the ARD prior on Z, not on W.

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
        self.precompute()

    def precompute(self):
        self.D = self.dim[0]
        self.factors_axis = 1

    def updateParameters(self, ix=None, ro=None):
        if ix is not None:
            print('stochastic inference not implemented for SW node')
            exit(1)
        # Collect expectations from other nodes
        Ztmp = self.markov_blanket["Z"].getExpectations()
        Z,ZZ = Ztmp["E"],Ztmp["E2"]
        tau = self.markov_blanket["Tau"].getExpectation(expand=True).copy()
        Y = self.markov_blanket["Y"].getExpectation().copy()
        if "AlphaW" not in self.markov_blanket:
            print("SW node not implemented wihtout ARD")
            exit(1)
        alpha = self.markov_blanket["AlphaW"].getExpectation(expand=True).copy()
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

        # Update each latent variable in turn
        for k in range(self.dim[1]):

            # Calculate intermediate steps
            term1 = (theta_lnE-theta_lnEInv)[:,k]
            term2 = 0.5*s.log(alpha[:,k])

            term3 = 0.5*s.log(s.dot(ZZ[:,k], tau) + alpha[:,k]) # good to modify
            # term4_tmp1 = ma.dot((tau*Y).T,Z[:,k]).data
            term4_tmp1 = s.dot((tau*Y).T,Z[:,k]) # good to modify

            # term4_tmp2 = ( tau * s.dot((Z[:,k]*Z[:,s.arange(self.dim[1])!=k].T).T, SW[:,s.arange(self.dim[1])!=k].T) ).sum(axis=0)
            term4_tmp2 = ( tau * s.dot((Z[:,k]*Z[:,s.arange(self.dim[1])!=k].T).T, SW[:,s.arange(self.dim[1])!=k].T) ).sum(axis=0) # good to modify

            term4_tmp3 = s.dot(ZZ[:,k].T,tau) + alpha[:,k] # good to modify (I REPLACE MA.DOT FOR S.DOT, IT SHOULD BE SAFE )

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
        self.Q.setParameters(mean_B0=s.zeros((self.D,self.dim[1])), var_B0=1./alpha, mean_B1=Qmean_S1, var_B1=Qvar_S1, theta=Qtheta )

    def calculateELBO(self):
        # import pdb; pdb.set_trace()
        # Collect parameters and expectations
        Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()
        S,WW = Qexp["EB"], Qexp["ENN"]
        Qvar = Qpar['var_B1']

        theta = self.markov_blanket["ThetaW"].getExpectations(expand=True)


        # Get ARD sparsity or prior variance
        if "AlphaW" in self.markov_blanket:
            alpha = self.markov_blanket["AlphaW"].getExpectations(expand=True).copy()
        else:
            print("Not implemented")
            exit(1)

        # Calculate ELBO for W
        # import pdb; pdb.set_trace()
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


class MuW_Node(UnivariateGaussian_Unobserved_Variational_Node):
    """ """

    def __init__(self, pmean, pvar, qmean, qvar, clusters, n_W, cluster_dic=None, qE=None, qE2=None):
        # compute dim from numbers of clusters (n_clusters * W)
        self.clusters = clusters
        self.D = len(self.clusters)
        self.n_clusters = len(np.unique(clusters))
        dim = (self.n_clusters, n_W)
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
        W = self.markov_blanket['W'].Q.getExpectation()

        if "SigmaAlphaW" in self.markov_blanket:
            if self.markov_blanket["SigmaAlphaW"].__class__.__name__ == "AlphaW_Node_mk":
                Alpha = self.markov_blanket['SigmaAlphaW'].getExpectation().copy()  # Notice that this Alpha is the ARD prior on W, not on W.
            else:
                Sigma = self.markov_blanket['SigmaAlphaW'].getExpectation().copy()
        else:
            Sigma = self.markov_blanket['W'].P.getParameters()["cov"]

        Qmean, Qvar = self.Q.getParameters()['mean'], self.Q.getParameters()['var']


        b = ("SigmaAlphaW" in self.markov_blanket) and (
                self.markov_blanket["SigmaAlphaW"].__class__.__name__ == "AlphaW_Node_mk")

        # update of the variance

        for c in range(self.n_clusters):
            mask = (self.clusters == c)
            if b:
                tmp = (Alpha[mask, :]).sum(axis=0)
            else:
                tmp = np.matrix.trace(Sigma[:, mask, mask])
            Qvar[c, :] = tmp
        Qvar += 1. / Ppar['var']
        Qvar = 1. / Qvar

        # update of the mean

        if b:
            tmp = W * Alpha
        else:
            # TODO : check below if we should ask the mask
            tmp = np.zeros(self.dim)
            for k in range(self.dim[1]):
                tmp[:, k] = np.dot(Sigma[k, :, :] - np.diag(Sigma[k, :, :]), W[:, k])

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
        latent_variables = self.markov_blanket['W'].getLvIndex()
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
