import scipy as s

from .variational_nodes import Variational_Node
from .basic_nodes import Constant_Node

"""
This module defines nodes that are a mix of variational and constant.
Examples:
- For some factors the spike and slab sparsity parameter (theta) is fixed, whereas for others is learnt.
- For some factors the weights can be made sparse whereas for others we might not want to do any sparsity (i.e. mean or covariates)
"""

class Mixed_ThetaW_Nodes_mk(Variational_Node, Constant_Node):
    """
    Class for a mixture of LearningThetaW and ConstantThetaW nodes.
    For a total of K factors, some (Klearn) will learn Theta whereas for the others (Kconst) it will be constant
        K = Klearn + Kconst
    """
    def __init__(self, LearnTheta, ConstTheta, idx):
        # Inputs:
        # - LearnTheta: ThetaW_Node_mk with dimensions (Klearn,)
        # - ConstTheta: ThetaW_Constant_Node_mk with dimensions (D,Kconst) or (Kconst,1) - NOT IMPLEMENTED YET -
        # - idx: list or numpy array indicating which factors are LearnTheta(idx=1. or idx=True) and which are ConstTheta(idx=0. or idx=False)
        self.constTheta = ConstTheta
        self.learnTheta = LearnTheta

        self.K = ConstTheta.dim[1] + LearnTheta.dim[0]
        self.D = ConstTheta.dim[0]

        self.idx = idx

    def addMarkovBlanket(self, **kargs):
        # SHOULD WE ALSO ADD MARKOV BLANKET FOR CONSTHTETA???
        self.learnTheta.addMarkovBlanket(**kargs)

    def getExpectations(self):
        # TODO fix for intercept etc:
        # 1 - make sure that getExpectations takes an expand argument which it doesnt use (to be compatible with other theta nodes )
        # 2 - make sure that instead of expanding here from theta_learn we just use expand=True as argument, to make it compatible to both per group and not
        # Get expectations from ConstTheta nodes (D,Kconst)
        Econst = self.constTheta.getExpectations().copy()

        # Get expectations from LearnTheta nodes and expand to (D,Kconst)
        Elearn = self.learnTheta.getExpectations().copy()
        Elearn["E"] = s.repeat(Elearn["E"][None,:], self.D, 0)
        Elearn["lnE"] = s.repeat(Elearn["lnE"][None,:], self.D, 0)
        Elearn["lnEInv"] = s.repeat(Elearn["lnEInv"][None,:], self.D, 0)

        # Concatenate expectations to (D,K)
        E = s.concatenate((Econst["E"], Elearn["E"]), axis=1)
        lnE = s.concatenate((Econst["lnE"], Elearn["lnE"]), axis=1)
        lnEInv = s.concatenate((Econst["lnEInv"], Elearn["lnEInv"]), axis=1)

        # Permute to the right order given by self.idx
        idx = s.concatenate((s.nonzero(1-self.idx)[0],s.where(self.idx)[0]), axis=0)
        E, lnE, lnEinv = E[:,idx], lnE[:,idx], lnEInv[:,idx]
        return dict({'E': E, 'lnE': lnE, 'lnEInv':lnEInv})

    def getExpectation(self):
        return self.getExpectations()['E']

    def updateExpectations(self):
        self.learnTheta.updateExpectations()

    def updateParameters(self):
        # the argument contains the indices of the non_annotated factors
        self.learnTheta.updateParameters(s.nonzero(self.idx)[0])

    def calculateELBO(self):
        return self.learnTheta.calculateELBO()

    def removeFactors(self, *idx):
        for i in idx:
            if self.idx[idx] == 1:
                self.learnTheta.removeFactors(s.where(i == s.nonzero(self.idx)[0])[0])
            else:
                self.constTheta.removeFactors(s.where(i == s.nonzero(1-self.idx)[0])[0])
            self.idx = self.idx[s.arange(self.K)!=i]
            self.K -= 1

    def sample(self, dist='P'):
        # sample from constant and learnt nodes
        const_samp = self.constTheta.sample().copy()
        learn_samp = self.learnTheta.sample().copy()
        learn_samp = s.repeat(learn_samp[None,:], self.D, 0)

        # concatenate and reorder
        samp = s.concatenate((const_samp, learn_samp), axis=1)
        idx = s.concatenate((s.nonzero(1-self.idx)[0],s.where(self.idx)[0]), axis=0)
        samp = samp[:,idx]

        # save and return
        self.samp = samp
        return self.samp

# TODO could remove
class Mixed_ThetaZ_Nodes_k(Variational_Node, Constant_Node):
    """
    Class for a mixture of LearningThetaZ and ConstantThetaZ nodes.
    For a total of K factors, some (Klearn) will learn Theta whereas for the others (Kconst) it will be constant
        K = Klearn + Kconst
    """
    def __init__(self, LearnTheta, ConstTheta, idx):
        # Inputs:
        # - LearnTheta: ThetaZ_Node_k with dimensions (Klearn,)
        # - ConstTheta: ThetaZ_Constant_Node_k with dimensions (D,Kconst) or (Kconst,1) - NOT IMPLEMENTED YET -
        # - idx: list or numpy array indicating which factors are LearnTheta(idx=1. or idx=True) and which are ConstTheta(idx=0. or idx=False)
        self.constTheta = ConstTheta
        self.learnTheta = LearnTheta

        self.K = ConstTheta.dim[1] + LearnTheta.dim[0]
        self.N = ConstTheta.dim[0]

        self.idx = idx

    def addMarkovBlanket(self, **kargs):
        # SHOULD WE ALSO ADD MARKOV BLANKET FOR CONSTHTETA???
        self.learnTheta.addMarkovBlanket(**kargs)

    def getExpectations(self):

        # Get expectations from ConstTheta nodes (N,Kconst)
        Econst = self.constTheta.getExpectations().copy()

        # Get expectations from LearnTheta nodes and expand to (N,Kconst)
        Elearn = self.learnTheta.getExpectations().copy()
        Elearn["E"] = s.repeat(Elearn["E"][None,:], self.N, 0)
        Elearn["lnE"] = s.repeat(Elearn["lnE"][None,:], self.N, 0)
        Elearn["lnEInv"] = s.repeat(Elearn["lnEInv"][None,:], self.N, 0)

        # Concatenate expectations to (N,K)
        E = s.concatenate((Econst["E"], Elearn["E"]), axis=1)
        lnE = s.concatenate((Econst["lnE"], Elearn["lnE"]), axis=1)
        lnEInv = s.concatenate((Econst["lnEInv"], Elearn["lnEInv"]), axis=1)

        # Permute to the right order given by self.idx
        idx = s.concatenate((s.nonzero(1-self.idx)[0],s.where(self.idx)[0]), axis=0)
        E, lnE, lnEinv = E[:,idx], lnE[:,idx], lnEInv[:,idx]
        return dict({'E': E, 'lnE': lnE, 'lnEInv':lnEInv})

    def getExpectation(self):
        return self.getExpectations()['E']

    def updateExpectations(self):
        self.learnTheta.updateExpectations()

    def updateParameters(self):
        # the argument contains the indices of the non_annotated factors
        self.learnTheta.updateParameters(s.nonzero(self.idx)[0])

    def calculateELBO(self):
        return self.learnTheta.calculateELBO()

    def removeFactors(self, *idx):
        for i in idx:
            if self.idx[idx] == 1:
                self.learnTheta.removeFactors(s.where(i == s.nonzero(self.idx)[0])[0])
            else:
                self.constTheta.removeFactors(s.where(i == s.nonzero(1-self.idx)[0])[0])
            self.idx = self.idx[s.arange(self.K)!=i]
            self.K -= 1

    def sample(self, dist='P'):
        # sample from constant and learnt nodes
        const_samp = self.constTheta.sample().copy()
        learn_samp = self.learnTheta.sample().copy()
        learn_samp = s.repeat(learn_samp[None,:], self.D, 0)

        # concatenate and reorder
        samp = s.concatenate((const_samp, learn_samp), axis=1)
        idx = s.concatenate((s.nonzero(1-self.idx)[0],s.where(self.idx)[0]), axis=0)
        samp = samp[:,idx]

        # save and return
        self.samp = samp
        return self.samp
