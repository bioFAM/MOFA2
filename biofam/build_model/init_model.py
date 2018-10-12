"""
Module to initialise a bioFAM model
"""

import scipy as s
import scipy.stats as stats
from sys import path
import sklearn.decomposition

from biofam.core.nodes import *

class initModel(object):
    def __init__(self, dim, data, lik):
        """
        PARAMETERS
        dim: dictionary with keyworded dimensionalities:
            N for the number of samples
            M for the number of views
            K for the number of factors or latent variables,
            D for the number of features per feature group, or view (a list)
        data: list of length M with ndarrays of dimensionality (N,Dm):
            observed data
        lik: list of strings with length M
            likelihood for each view, choose from ('gaussian','poisson','bernoulli')
        """
        self.data = data
        self.lik = lik
        self.N = dim["N"]
        self.K = dim["K"]
        self.M = dim["M"]
        self.D = dim["D"]

        self.nodes = {}

    def initZ(self, pmean=0., pvar=1., qmean='random', qvar=1., qE=None, qE2=None, covariates=None,
        scale_covariates=None, precompute_pcovinv=True):
        """Method to initialise the latent variables

        PARAMETERS
        ----------
        pmean: mean of the prior distribution
        pcov: covariance of the prior distribution # NOTE was changed to pvar -> univariate
        qmean: initialisation of the mean of the variational distribution
            "random" for a random initialisation sampled from a standard normal distribution
            "pca" for an initialisation based on PCA
            "orthogonal" for a random initialisation with orthogonal factors
        qvar: initial value of the variance of the variational distribution
        qE: initial value of the expectation of the variational distribution
        qE2: initial value of the second moment of the variational distribution
        covariates: covariates to be included as non-updated factors
            None if no covariates are present, or a ndarray covariates with dimensions (N,Kcovariates)
        scale_covariates: scale covariates to zero-mean and unit variance to match the prior?
            None if no covariates are present, or a ndarray with dimensions (Kcov,) indicating which covariates to scale
        precompute_pcovinv: precompute the inverse of the covariance matrice of the prior of Z
        """

        ## Initialise prior distribution (P)

        # mean
        pmean = s.ones((self.N, self.K)) * pmean

        # TODO add sanity check if not float (dim of the matrices)
        # if isinstance(pcov, (int, float)):
        #     pcov = [s.sparse.eye(self.N) * pcov] * self.K

        ## Initialise variational distribution (Q)

        # variance
        qvar = s.ones((self.N, self.K)) * qvar

        # mean
        if qmean is not None:
            if isinstance(qmean, str):

                # Random initialisation
                if qmean == "random":
                    qmean = stats.norm.rvs(loc=0, scale=1, size=(self.N, self.K))

                # Random initialisation with orthogonal factors
                elif qmean == "orthogonal":
                    pca = sklearn.decomposition.PCA(n_components=self.K, copy=True, whiten=True)
                    pca.fit(stats.norm.rvs(loc=0, scale=1, size=(self.N, 9999)).T)
                    qmean = pca.components_.T

                # PCA initialisation
                elif qmean == "pca":
                    pca = sklearn.decomposition.PCA(n_components=self.K, copy=True, whiten=True)
                    pca.fit(s.concatenate(self.data, axis=0).T)
                    qmean = pca.components_.T

            elif isinstance(qmean, s.ndarray):
                assert qmean.shape == (self.N, self.K), "Wrong shape for the expectation of the Q distribution of Z"

            elif isinstance(qmean, (int, float)):
                qmean = s.ones((self.N, self.K)) * qmean

            else:
                print("Wrong initialisation for Z")
                exit()

        # Add covariates
        if covariates is not None:

            # Sanity checks
            assert scale_covariates != None, "If you use covariates you also need define model_opts['scale_covariates']"

            # Select indices for covariates
            idx_covariates = s.array(range(covariates.shape[1]))

            # Center and scale the covariates to match the prior distribution N(0,1)
            # Notice that if the spherical prior distribution is changed, this should also be modified
            scale_covariates = s.array(scale_covariates)
            covariates[:, scale_covariates] = (covariates[:, scale_covariates] - s.nanmean(
                covariates[:, scale_covariates], axis=0)) / s.nanstd(covariates[:, scale_covariates], axis=0)

            # Set to zero the missing values in the covariates
            covariates[s.isnan(covariates)] = 0.
            qmean[:, idx_covariates] = covariates

            # Remove prior and variational distributions from the covariates
            # pcov[idx_covariates] = s.nan
            #pvar[:, idx_covariates] = s.nan
            qvar[:, idx_covariates] = s.nan

        else:
            idx_covariates = None

        # Initialise the node
        self.nodes["Z"] = Z_Node(
            dim=(self.N, self.K),
            pmean=pmean, pvar=pvar,
            qmean=qmean, qvar=qvar,
            qE=qE, qE2=qE2,
            idx_covariates=idx_covariates
        )

        # self.nodes["Z"] = Z_Node(
        #     dim=(self.N, self.K),
        #     pmean=pmean, pcov=pcov,
        #     qmean=qmean, qvar=qvar,
        #     qE=qE, qE2=qE2,
        #     idx_covariates=idx_covariates,
        #     precompute_pcovinv=precompute_pcovinv
        # )

    def initSZ(self, pmean_T0=0., pmean_T1=0., pvar_T0=1., pvar_T1=1., ptheta=1., qmean_T0=0., qmean_T1='random', qvar_T0=1.,
        qvar_T1=1., qtheta=1., qEZ_T0=None, qEZ_T1=None, qET=None):
        """Method to initialise sparse factors with a (reparametrised) spike and slab prior
        PARAMETERS
        ----------
        (...)
        """

        # TODO:
        # - add different ways to initialise the expectations of the posterior: random, orthogonal, PCA
        # - enable covaraites, currently intercept does not work


        ## Initialise prior distribution (P)

        ## Initialise variational distribution (Q)
        if isinstance(qmean_T1, str):

            if qmean_T1 == "random":
                qmean_T1_tmp = stats.norm.rvs(loc=0, scale=1, size=(self.N, self.K))
            else:
                print("%s initialisation not implemented for Z" % qmean_T1)
                exit()

        elif isinstance(qmean_T1, s.ndarray):
            assert qmean_T1.shape == (self.N, self.K), "Wrong dimensionality"

        elif isinstance(qmean_T1, (int, float)):
            qmean_T1_tmp = s.ones((self.N, self.K)) * qmean_T1

        else:
            print("Wrong initialisation for Z")
            exit(1)

        self.nodes["Z"] = SZ_Node(
            dim=(self.N, self.K),

            ptheta=ptheta,
            pmean_T0=pmean_T0,
            pvar_T0=pvar_T0,
            pmean_T1=pmean_T1,
            pvar_T1=pvar_T1,

            qtheta=qtheta,
            qmean_T0=qmean_T0,
            qvar_T0=qvar_T0,
            qmean_T1=qmean_T1_tmp,
            qvar_T1=qvar_T1,

            qET=qET,
            qEZ_T0=qEZ_T0,
            qEZ_T1=qEZ_T1,
        )

    def initW(self, pmean=0., pvar=1., qmean=0, qvar=1.,
        qE=None, qE2=None, covariates=None,
        scale_covariates=None):
        """Method to initialise the weights
        PARAMETERS
        ----------
        pmean: mean of the prior distribution
        pvar: variance of the prior distribution
        qmean: initial value of the mean of the variational distribution
        qvar: initial value of the variance of the variational distribution
        qE: initial value of the expectation of the variational distribution
        qE2: initial value of the second moment of the variational distribution
        (NOT FUNCTIONAL) covariates: covariates to be included as non-updated weights
            None if no covariates are present, or a ndarray covariates with dimensions (N,Kcovariates)
        (NOT FUNCTIONAL) scale_covariates: scale covariates to zero-mean and unit variance to match the prior?
            None if no covariates are present, or a ndarray with dimensions (Kcov,) indicating which covariates to scale
        precompute_pcovinv: precompute the inverse of the covariance matrice of the prior of W
        """

        # Precompute the inverse of the covariance matrix of the gaussian prior of W
        # if precompute_pcovinv is None:
        #     precompute_pcovinv = [True] * self.M

        # Add covariates
        if covariates is not None:
            assert scale_covariates != None, "If you use covariates also define data_opts['scale_covariates']"

            # Select indices for covariates
            idx_covariates = s.array(range(covariates.shape[1]))

            # Center and scale the covariates to match the prior distribution N(0,1)
            scale_covariates = s.array(scale_covariates)
            covariates[:, scale_covariates] = (covariates[:, scale_covariates] - s.nanmean(
                covariates[:, scale_covariates], axis=0)) / s.nanstd(covariates[:, scale_covariates], axis=0)

            # Set to zero the missing values in the covariates
            covariates[s.isnan(covariates)] = 0.
        else:
            idx_covariates = None

        W_list = [None] * self.M
        for m in range(self.M):

            ## Initialise prior distribution (P) ##

            # mean
            pmean_m = s.ones((self.D[m], self.K)) * pmean

            # covariance
            # if isinstance(pcov, (int, float)):
            #     pcov_mk = s.sparse.eye(self.D[m]) * pcov
            #     pcov_m = [pcov_mk] * self.K

            ## Initialise variational distribution (Q) ##

            # variance
            qvar_m = s.ones((self.D[m], self.K)) * qvar

            # mean
            if qmean is not None:
                if isinstance(qmean, str):

                    # Random initialisation
                    if qmean == "random":
                        qmean_m = stats.norm.rvs(loc=0, scale=1., size=(self.D[m], self.K))

                    # Random and orthogonal initialisation
                    # elif qmean == "orthogonal":
                    #     pca = sklearn.decomposition.PCA(n_components=self.K, copy=True, whiten=True)
                    #     pca.fit(stats.norm.rvs(loc=0, scale=1, size=(self.D[m], 9999)).T)
                    #     qmean_m = pca.components_.T

                    # PCA initialisation
                    # elif qmean == "pca":
                    #     pca = sklearn.decomposition.PCA(n_components=self.K, copy=True, whiten=True)
                    #     pca.fit(s.concatenate(self.data, axis=0).T)
                    #     qmean_m = pca.components_.T

                elif isinstance(qmean, s.ndarray):
                    assert qmean.shape == (
                    self.D[m], self.K), "Wrong shape for the expectation of the Q distribution of W"
                    qmean_m = qmean

                elif isinstance(qmean, (int, float)):
                    qmean_m = s.ones((self.D[m], self.K)) * qmean

                else:
                    print("Wrong initialisation for W")
                    exit()

            else:
                qmean_m = None

            # Add covariates
            if covariates is not None:
                print("Covariates not implemented")
                exit()

                qmean_m[:, idx_covariates] = covariates

                # Remove prior and variational distributions from the covariates
                # pcov_m[idx_covariates] = s.nan
                qvar_m[:, idx_covariates] = s.nan  # MAYBE SET IT TO 0

            # Initialise the node
            W_list[m] = W_Node(dim=(self.D[m], self.K), pmean=pmean_m, pvar=pvar, qmean=qmean_m, qvar=qvar_m, qE=qE, qE2=qE2,
                               idx_covariates=idx_covariates)

        self.nodes["W"] = Multiview_Variational_Node(self.M, *W_list)

    def initSW(self, pmean_S0=0., pmean_S1=0., pvar_S0=1., pvar_S1=1., ptheta=1.,
        qmean_S0=0., qmean_S1=0, qvar_S0=1., qvar_S1=1., qtheta=1.,
        qEW_S0=None, qEW_S1=None, qES=None):
        """Method to initialise sparse weights with a (reparametrised) spike and slab prior

        PARAMETERS
        ----------
        (...)
        """

        W_list = [None]*self.M
        for m in range(self.M):

            ## Initialise prior distribution (P)

            ## Initialise variational distribution (Q)
            if isinstance(qmean_S1,str):

                if qmean_S1 == "random": # random
                    qmean_S1_tmp = stats.norm.rvs(loc=0, scale=1., size=(self.D[m],self.K))
                else:
                    print("%s initialisation not implemented for W" % qmean_S1)
                    exit()

            elif isinstance(qmean_S1,s.ndarray):
                assert qmean_S1.shape == (self.D[m],self.K), "Wrong dimensionality"

            elif isinstance(qmean_S1,(int,float)):
                qmean_S1_tmp = s.ones((self.D[m],self.K)) * qmean_S1

            else:
                print("Wrong initialisation for W")
                exit(1)

            W_list[m] = SW_Node(
                dim=(self.D[m],self.K),

                ptheta=ptheta,
                pmean_S0=pmean_S0,
                pvar_S0=pvar_S0,
                pmean_S1=pmean_S1,
                pvar_S1=pvar_S1,

                qtheta=qtheta,
                qmean_S0=qmean_S0,
                qvar_S0=qvar_S0,
                qmean_S1=qmean_S1_tmp,
                qvar_S1=qvar_S1,

                qES=qES,
                qEW_S0=qEW_S0,
                qEW_S1=qEW_S1,
            )

        self.nodes["W"] = Multiview_Variational_Node(self.M, *W_list)


    def initAlphaZ(self, groups, pa=1e-14, pb=1e-14, qa=1., qb=1., qE=None, qlnE=None):
        """Method to initialise the ARD prior on Z per sample group

        PARAMETERS
        ----------
        groups: dictionary
            a dictionariy with mapping between sample names (keys) and samples_groups (values)
        pa: float
            'a' parameter of the prior distribution
        pb :float
            'b' parameter of the prior distribution
        qb: float
            initialisation of the 'b' parameter of the variational distribution
        qE: float
            initial expectation of the variational distribution
        qlnE: float
            initial log expectation of the variational distribution
        """

        # Sanity checks
        assert len(groups) == self.N, 'sample groups labels do not match number of samples'

        # convert groups into integers from 0 to n_groups and keep the corresponding group names in groups_dic
        tmp = np.unique(groups, return_inverse=True)
        groups_dic = tmp[0]
        groups_ix = tmp[1]

        n_group = len(np.unique(groups_ix))
        assert len(groups_dic) == n_group, 'problem in np.unique'

        self.nodes["AlphaZ"] = AlphaZ_Node(
            dim=(n_group, self.K),
            pa=pa, pb=pb,
            qa=qa, qb=qb,
            groups=groups_ix,
            groups_dic=groups_dic,
            qE=qE, qlnE=qlnE
        )

    def initAlphaW(self, pa=1e-14, pb=1e-14, qa=1., qb=1., qE=None, qlnE=None):
        """Method to initialise the ARD prior on W

        PARAMETERS
        ----------
        pa: float
            'a' parameter of the prior distribution
        pb :float
            'b' parameter of the prior distribution
        qa: float
            initialisation of the 'b' parameter of the variational distribution
        qb: float
            initialisation of the 'b' parameter of the variational distribution
        qE: float
            initial expectation of the variational distribution
        qlnE: float
            initial log expectation of the variational distribution
        """

        alpha_list = [None]*self.M
        for m in range(self.M):
            alpha_list[m] = AlphaW_Node(dim=(self.K,), pa=pa, pb=pb, qa=qa, qb=qb, qE=qE, qlnE=qlnE)
        self.nodes["AlphaW"] = Multiview_Variational_Node(self.M, *alpha_list)

    def initTau(self, groups, pa=1., pb=1., qa=1., qb=1., qE=None, on='features'):
        """Method to initialise the precision of the noise

        PARAMETERS
        ----------
        pa: float
            'a' parameter of the prior distribution
        pb :float
            'b' parameter of the prior distribution
        qb: float
            initialisation of the 'b' parameter of the variational distribution
        qE: float
            initial expectation of the variational distribution
        on: string
            'features' to define a noise per feature
            'samples' to define a noise per sample
        """
        pa = 1.
        pb = 1.

        tau_list = [None]*self.M

        # Sanity checks
        assert len(groups) == self.N, 'sample groups labels do not match number of samples'

        # convert groups into integers from 0 to n_groups and keep the corresponding group names in groups_dic
        tmp = np.unique(groups, return_inverse=True)
        groups_dic = tmp[0]
        groups_ix = tmp[1]

        n_group = len(np.unique(groups_ix))
        assert len(groups_dic) == n_group, 'problem in np.unique'

        for m in range(self.M):

            # Poisson noise model for count data
            if self.lik[m] == "poisson":
                tmp = 0.25 + 0.17*s.amax(self.data[m],axis=0)
                tau_list[m] = Constant_Node(dim=(self.N, self.D[m]), value=tmp)

            # Bernoulli noise model for binary data
            elif self.lik[m] == "bernoulli":
                # tau_list[m] = Constant_Node(dim=(self.D[m],), value=s.ones(self.D[m])*0.25)
                # tau_list[m] = Tau_Jaakkola(dim=(self.D[m],), value=0.25)
                tau_list[m] = Tau_Jaakkola(dim=((self.N, self.D[m])), value=1.)

            # Binomial noise model for proportion data
            elif self.lik[m] == "binomial":
                print("Not implemented")
                exit()
                # tmp = 0.25*s.amax(self.data["tot"][m],axis=0)
                # tau_list[m] = Constant_Node(dim=(self.N, self.D[m]), value=tmp)

            # Gaussian noise model for continuous data
            elif self.lik[m] == "gaussian":

                # Noise per samples
                if on == 'samples':
                    tau_list[m] = TauN_Node(dim=(self.N,), pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
                    print("Using TauN noise!")
                # Noise per feature
                elif on == 'features':
                    tau_list[m] = TauD_Node(dim=(n_group, self.D[m]), pa=pa, pb=pb, qa=qa, qb=qb, groups=groups_ix, groups_dic=groups_dic, qE=qE)
                else:
                    print('did not understand noise option on =', on)
                    exit(1)

        self.nodes["Tau"] = Multiview_Mixed_Node(self.M, *tau_list)

    # TODO: make independent of noise but the problem for that is that precompute() needs to know and is called
    # before the markov_blanket is defined
    def initY(self, noise_on='features'):
        """Method to initialise the observations"""
        Y_list = [None]*self.M
        for m in range(self.M):
            if self.lik[m]=="gaussian":
                Y_list[m] = Y_Node(dim=(self.N,self.D[m]), value=self.data[m], noise_on=noise_on)
            elif self.lik[m]=="poisson":
                Y_list[m] = Poisson_PseudoY(dim=(self.N,self.D[m]), obs=self.data[m], E=None, noise_on=noise_on)
            elif self.lik[m]=="bernoulli":
                # Y_list[m] = Bernoulli_PseudoY(dim=(self.N,self.D[m]), obs=self.data[m], E=None)
                Y_list[m] =  Bernoulli_PseudoY_Jaakkola(dim=(self.N,self.D[m]), obs=self.data[m], E=None, noise_on=noise_on)
        self.Y = Multiview_Mixed_Node(self.M, *Y_list)
        self.nodes["Y"] = self.Y

    def initThetaZ(self, groups, pa=1, pb=1, qa=1., qb=1., qE=None, Klearn = None):
        """Method to initialise the ARD prior on Z per sample group

        PARAMETERS
        ----------
        groups: dictionary
            a dictionariy with mapping between sample names (keys) and samples_groups (values)
        pa: float
            'a' parameter of the prior distribution
        pb :float
            'b' parameter of the prior distribution
        qb: float
            initialisation of the 'b' parameter of the variational distribution
        qE: float
            initial expectation of the variational distribution
        qlnE: float
            initial log expectation of the variational distribution
        K: number of factors for which we learn theta. If no argument is given, we'll just use the
            total number of factors
        """

        # Sanity checks
        assert len(groups) == self.N, 'sample groups labels do not match number of samples'

        # convert groups into integers from 0 to n_groups and keep the corresponding group names in groups_dic
        tmp = np.unique(groups, return_inverse=True)
        groups_dic = tmp[0]
        groups_ix = tmp[1]

        n_group = len(np.unique(groups_ix))
        assert len(groups_dic) == n_group, 'problem in np.unique'

        # number of factors for which we learn theta
        if Klearn is None:
            Klearn = self.K

        self.nodes["ThetaZ"] = ThetaZ_Node(
            dim=(n_group, Klearn),
            pa=pa, pb=pb,
            qa=qa, qb=qb,
            groups=groups_ix,
            groups_dic=groups_dic,
            qE=qE)

    def initThetaW(self, pa=1., pb=1., qa=1., qb=1., qE=None):
        """Method to initialise the sparsity parameter of the spike and slab weights

        PARAMETERS
        ----------
         pa: float
            'a' parameter of the prior distribution
         pb :float
            'b' parameter of the prior distribution
         qb: float
            initialisation of the 'b' parameter of the variational distribution
         qE: float
            initial expectation of the variational distribution
        """
        Theta_list = [None] * self.M
        for m in range(self.M):
            Theta_list[m] = ThetaW_Node(dim=(self.K,), pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        self.nodes["ThetaW"] = Multiview_Variational_Node(self.M, *Theta_list)

    def initExpectations(self, *nodes):
        """ Method to initialise all expectations """
        for node in nodes:
            self.nodes[node].updateExpectations()

    def getNodes(self):
        """ Get method to return the nodes """
        return self.nodes
        #return { k:v for (k,v) in self.nodes.items()}
