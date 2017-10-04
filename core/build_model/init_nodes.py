
"""
Module to initalise the nodes

Z:
    MuZ:

SW:
    Alpha:

Tau:

Y:

Theta:
    ThetaConst
    ThetaLearn
    ThetaMixed
"""

import scipy as s
import scipy.stats as stats
from sys import path
import sklearn.decomposition


# Define which nodes to import
from biofam.nodes.basic_nodes import *
from biofam.nodes.multiview_nodes import *
from biofam.nodes.nongaussian_nodes import *
from biofam.nodes.Y_nodes import Y_Node
from biofam.nodes.Z_nodes import Z_Node, MuZ_Node
from biofam.nodes.W_nodes import SW_Node
from biofam.nodes.Alpha_nodes import AlphaW_Node_mk
from biofam.nodes.Tau_nodes import Tau_Node
from biofam.nodes.Theta_nodes import Theta_Node, Theta_Constant_Node


class initModel(object):
    def __init__(self, dim, data, lik, seed=None):
        """
        PARAMETERS
        ----------
         dim: dictionary
            keyworded dimensionalities: N for the number of samples, M for the number of views, K for the number of latent variables, D for the number of features per view (a list)
         data: list of ndarrays of length M: 
            observed data
         lik: list of strings 
            likelihood for each view
        """
        self.data = data
        self.lik = lik
        self.N = dim["N"]
        self.K = dim["K"]
        self.M = dim["M"]
        self.D = dim["D"]

        self.nodes = {}

        # Set the seed
        s.random.seed(seed)

    def initZ(self, pmean, pvar, qmean, qvar, qE=None, qE2=None, covariates=None, scale_covariates=None):
        """Method to initialise the latent variables

        PARAMETERS
        ----------
        pmean:
        pvar:
        qmean
        qvar
        qE
        qE2
        covariates: nd array
            matrix of covariates with dimensions (nsamples,ncovariates)
        scale_covariates: 
        """

        # Initialise mean of the Q distribution
        if qmean is not None:
            if isinstance(qmean,str):
                if qmean == "random": # Random initialisation of latent variables
                    qmean = stats.norm.rvs(loc=0, scale=1, size=(self.N,self.K))

                elif qmean == "orthogonal": # Latent variables are initialised randomly but ensuring orthogonality
                    pca = sklearn.decomposition.PCA(n_components=self.K, copy=True, whiten=True)
                    pca.fit(stats.norm.rvs(loc=0, scale=1, size=(self.N,9999)).T)
                    qmean = pca.components_.T

                elif qmean == "pca": # Latent variables are initialised from PCA in the concatenated matrix
                    pca = sklearn.decomposition.PCA(n_components=self.K, copy=True, whiten=True)
                    pca.fit(s.concatenate(self.data,axis=0).T)
                    qmean = pca.components_.T

            elif isinstance(qmean,s.ndarray):
                assert qmean.shape == (self.N,self.K)

            elif isinstance(qmean,(int,float)):
                qmean = s.ones((self.N,self.K)) * qmean

            else:
                print("Wrong initialisation for Z")
                exit()

        # Add covariates
        if covariates is not None:
            assert scale_covariates != None, "If you use covariates also define data_opts['scale_covariates']"

            # Select indices for covaraites
            idx_covariates = s.array(range(covariates.shape[1]))

            # Center and scale the covariates to match the prior distribution N(0,1)
            # to-do: this needs to be improved to take the particular mean and var into account
            # covariates[scale_covariates] = (covariates - covariates.mean(axis=0)) / covariates.std(axis=0)
            scale_covariates = s.array(scale_covariates)
            covariates[:,scale_covariates] = (covariates[:,scale_covariates] - s.nanmean(covariates[:,scale_covariates], axis=0)) / s.nanstd(covariates[:,scale_covariates], axis=0)

            # Set to zero the missing values in the covariates
            covariates[s.isnan(covariates)] = 0.
            qmean[:,idx_covariates] = covariates
        else:
            idx_covariates = None


        # Initialise the node
        # self.Z = Constant_Node(dim=(self.N,self.K), value=qmean)
        self.Z = Z_Node(dim=(self.N,self.K),
                        pmean=s.ones((self.N,self.K))*pmean,
                        pvar=s.ones((self.K,))*pvar,
                        qmean=s.ones((self.N,self.K))*qmean,
                        qvar=s.ones((self.N,self.K))*qvar,
                        qE=qE, qE2=qE2,
                        idx_covariates=idx_covariates)
        self.nodes["Z"] = self.Z

    def initSW(self, pmean_S0, pmean_S1, pvar_S0, pvar_S1, ptheta, qmean_S0, qmean_S1, qvar_S0, qvar_S1, qtheta, qEW_S0, qEW_S1, qES):
        """Method to initialise the spike-slab variable (product of bernoulli and gaussian variables)

        PARAMETERS
        ----------
        """
        SW_list = [None]*self.M
        for m in range(self.M):

           # Initialise first moment
            if isinstance(qmean_S1[m],str):
                if qmean_S1[m] == "random":
                    qmean_S1[m] = stats.norm.rvs(loc=0, scale=1, size=(self.D[m],self.K))
                else:
                    print("%s initialisation not implemented for SW" % qmean_S1[m])
                    exit()
            elif isinstance(qmean_S1[m],s.ndarray):
                assert qmean_S1[m].shape == (self.D[m],self.K), "Wrong dimensionality"
            elif isinstance(qmean_S1[m],(int,float)):
                qmean_S1[m] = s.ones((self.D[m],self.K)) * qmean_S1[m]
            else:
                print("Wrong initialisation for SW")
                exit()

            SW_list[m] = SW_Node(
                dim=(self.D[m],self.K),

                ptheta=ptheta[m],
                pmean_S0=pmean_S0[m],
                pvar_S0=pvar_S0[m],
                pmean_S1=pmean_S1[m],
                pvar_S1=pvar_S1[m],

                qtheta=qtheta[m],
                qmean_S0=qmean_S0[m],
                qvar_S0=qvar_S0[m],
                qmean_S1=qmean_S1[m],
                qvar_S1=qvar_S1[m],

                qES=qES[m],
                qEW_S0=qEW_S0[m],
                qEW_S1=qEW_S1[m],
                )

        self.SW = Multiview_Variational_Node(self.M, *SW_list)
        self.nodes["SW"] = self.SW

    def initAlphaW_mk(self, pa, pb, qa, qb, qE):

        """Method to initialise the precision of the group-wise ARD prior

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
        
        alpha_list = [None]*self.M
        for m in range(self.M):
            alpha_list[m] = AlphaW_Node_mk(dim=(self.K,), pa=pa[m], pb=pb[m], qa=qa[m], qb=qb[m], qE=qE[m])
            # alpha_list[m] = Constant_Node(dim=(self.K,), value=qE[m])
            # alpha_list[m].factors_axis = 0
        self.AlphaW = Multiview_Variational_Node(self.M, *alpha_list)
        # self.AlphaW = Multiview_Constant_Node(self.M, *alpha_list)
        self.nodes["AlphaW"] = self.AlphaW

    def initTau(self, pa, pb, qa, qb, qE):
        # Method to initialise the precision of the noise
        # Inputs:
        #  pa (float): 'a' parameter of the prior distribution
        #  pb (float): 'b' parameter of the prior distribution
        #  qb (float): initialisation of the 'b' parameter of the variational distribution
        #  qE (float): initial expectation of the variational distribution
        tau_list = [None]*self.M
        for m in range(self.M):
            if self.lik[m] == "poisson":
                tmp = 0.25 + 0.17*s.amax(self.data[m],axis=0)
                tau_list[m] = Constant_Node(dim=(self.D[m],), value=tmp)
            elif self.lik[m] == "bernoulli":
                # tmp = s.ones(self.D[m])*0.25
                # tau_list[m] = Constant_Node(dim=(self.D[m],), value=tmp)
                # tau_list[m] = Tau_Jaakkola(dim=(self.D[m],), value=0.25)
                tau_list[m] = Tau_Jaakkola(dim=((self.N,self.D[m])), value=1.)
            elif self.lik[m] == "binomial":
                tmp = 0.25*s.amax(self.data["tot"][m],axis=0)
                tau_list[m] = Constant_Node(dim=(self.D[m],), value=tmp)
            elif self.lik[m] == "gaussian":
                tau_list[m] = Tau_Node(dim=(self.D[m],), pa=pa[m], pb=pb[m], qa=qa[m], qb=qb[m], qE=qE[m])
            elif self.lik[m] == "warp":
                tau_list[m] = Tau_Node(dim=(self.D[m],), pa=pa[m], pb=pb[m], qa=qa[m], qb=qb[m], qE=qE[m])
        self.Tau = Multiview_Mixed_Node(self.M,*tau_list)
        self.nodes["Tau"] = self.Tau

    def initY(self):
        # Method to initialise the observed data
        Y_list = [None]*self.M
        for m in range(self.M):
            if self.lik[m]=="gaussian":
                Y_list[m] = Y_Node(dim=(self.N,self.D[m]), value=self.data[m])
            elif self.lik[m]=="poisson":
                # tmp = stats.norm.rvs(loc=0, scale=1, size=(self.N,self.D[m]))
                Y_list[m] = Poisson_PseudoY(dim=(self.N,self.D[m]), obs=self.data[m], E=None)
            elif self.lik[m]=="bernoulli":
                # Y_list[m] = Bernoulli_PseudoY(dim=(self.N,self.D[m]), obs=self.data[m], E=None)
                # tmp = stats.norm.rvs(loc=0, scale=1, size=(self.N,self.D[m]))
                Y_list[m] =  Bernoulli_PseudoY_Jaakkola(dim=(self.N,self.D[m]), obs=self.data[m], E=None)
                # Y_list[m] =  Bernoulli_PseudoY_Jaakkola(dim=(self.N,self.D[m]), obs=self.data[m], E=self.data[m])
            elif self.lik[m]=="warp":
                print "Not implemented"
                exit()
                # Y_list[m] = Warped_PseudoY_Node(dim=(self.N,self.D[m]), obs=self.data[m], func_type='tanh', I=3, E=None)
                # Y_list[m] = Warped_PseudoY_Node(dim=(self.N,self.D[m]), obs=self.data[m], func_type='logistic', I=1, E=None)
        self.Y = Multiview_Mixed_Node(self.M, *Y_list)
        self.nodes["Y"] = self.Y

    def initThetaMixed(self, pa, pb, qa, qb, qE, learnTheta):
        # Method to initialie a general theta node
        # Inputs:
        #  pa (float): 'a' parameter of the prior distribution
        #  pb (float): 'b' parameter of the prior distribution
        #  qb (float): initialisation of the 'b' parameter of the variational distribution
        #  qE (float): initial expectation of the variational distribution
        #  learnTheta (binary): list with binary matrices with dim (D[m],K)

        Theta_list = [None] * self.M
        for m in range(self.M):
            
            # Initialise constant node
            Kconst = learnTheta[m]==0
            if Kconst.sum() == 0:
                ConstThetaNode = None
            else:
                # ConstThetaNode = Theta_Constant_Node(dim=(self.D[m],s.sum(Kconst),), value=s.repeat(qE[m][:,Kconst][None,:], self.D[m], 0), N_cells=1.)
                ConstThetaNode = Theta_Constant_Node(dim=(self.D[m],s.sum(Kconst),), value=qE[m][:,Kconst], N_cells=1)
                Theta_list[m] = ConstThetaNode

            # Initialise non-constant node
            Klearn = learnTheta[m]==1
            if Klearn.sum() == 0:
                LearnThetaNode = None
            else:
                # FOR NOW WE JUST TAKE THE FIRST ROW BECAUSE IT IS EXPANDED. IT IS UGLY AS HELL
                LearnThetaNode = Theta_Node(dim=(s.sum(Klearn),), pa=pa[m][Klearn], pb=pb[m][Klearn], qa=qa[m][Klearn], qb=qb[m][Klearn], qE=qE[m][0,Klearn])
                Theta_list[m] = LearnThetaNode

            # Initialise mixed node
            if (ConstThetaNode is not None) and (LearnThetaNode is not None):
                Theta_list[m] = Mixed_Theta_Nodes(LearnTheta=LearnThetaNode, ConstTheta=ConstThetaNode, idx=learnTheta[m])

        self.Theta = Multiview_Mixed_Node(self.M, *Theta_list)
        self.nodes["Theta"] = self.Theta

    def initThetaLearn(self, pa, pb, qa, qb, qE):
        # Method to initialise the theta node
        Theta_list = [None] * self.M
        for m in range(self.M):
            Theta_list[m] = Theta_Node(dim=(self.K,), pa=pa[m], pb=pb[m], qa=qa[m], qb=qb[m], qE=qE[m][0,:])
        self.Theta = Multiview_Variational_Node(self.M, *Theta_list)
        self.nodes["Theta"] = self.Theta

    def initThetaConst(self, value):
        # Method to initialise a constant theta node
        Theta_list = [None] * self.M
        for m in range(self.M):
            Theta_list[m] = Theta_Constant_Node(dim=(self.D[m],self.K,), value=value[m], N_cells=1.)
        self.Theta = Multiview_Constant_Node(self.M, *Theta_list)
        self.nodes["Theta"] = self.Theta

    def initMuZ(self, clusters=None, pmean=0, pvar=1, qmean=0, qvar=1, qE=None):
        if clusters is None:
            print "Not implemented"
            exit()
            clusters = s.zeros(self.N, int)
        self.MuZ = MuZ_Node(pmean, pvar, qmean, qvar, clusters, self.K)
        # self.Clusters = Constant_Node(pmean, pvar, qmean, qvar, clusters, self.K)
        self.nodes['MuZ'] = self.MuZ

    def initExpectations(self, *nodes):
        # Method to initialise the expectations of some nodes
        for node in nodes:
            self.nodes[node].updateExpectations()

    def getNodes(self):
        return { k:v for (k,v) in self.nodes.items()}
