"""
Module to build a bioFAM model
"""

from sys import path
from time import time,sleep
import numpy as np
import scipy as s
#from joblib import Parallel, delayed

from biofam.core.BayesNet import *
from biofam.build_model.init_model import initModel
from biofam.build_model.utils import *

class buildModel(object):
    def __init__(self, data, data_opts, model_opts, dimensionalities):
        self.data = data
        self.data_opts = data_opts
        self.model_opts = model_opts
        self.dim = dimensionalities
        pass

    def createMarkovBlankets(self):
        """ Define the markov blankets """
        pass

    def createSchedule(self):
        """ Define the schedule of updates """
        pass

    def build_nodes(self):
        """ Build all nodes """
        pass

    def get_nodes(self):
        """ Get all nodes """
        assert hasattr(self, 'net'), "BayesNet has not been created"
        return self.net.getNodes()


class buildBiofam(buildModel):
    def  __init__(self, data, data_opts, model_opts, dimensionalities):
        buildModel.__init__(self, data, data_opts, model_opts, dimensionalities)

        # create an instance of initModel
        self.init_model = initModel(self.dim, self.data, self.model_opts["likelihoods"])

        # Build all nodes
        self.build_nodes()

        # Define markov blankets
        self.createMarkovBlankets()

        # Define training schedule
        self.createSchedule()

        # Create BayesNet class
        self.createBayesNet()

    def build_nodes(self):
        """ Method to build all nodes """

        # Build general nodes
        self.build_Y()
        self.build_Z()
        self.build_W()
        self.build_Tau()

        # define ARD sparsity per sample group (on Z)
        if self.model_opts['ard_z']:
            self.build_AlphaZ()

        # define ARD sparsity per feature_group, or view (on W)
        if self.model_opts['ard_w']:
            self.build_AlphaW()

        # define feature-wise spike and slab sparsity in Z
        if self.model_opts['sl_z']:
            self.build_ThetaZ()

        # define feature-wise spike and slab sparsity in W
        if self.model_opts['sl_w']:
            self.build_ThetaW()

    def build_Y(self):
        """ Build node Y for the observations """
        self.init_model.initY()

    def build_Z(self):
        """ Build node Z for the factors or latent variables """
        if self.model_opts['sl_z']:
            self.init_model.initSZ(qmean_T1="random")
        else:
            # TODO change Z node so that we dont use a multivariate prior when no covariance structure
            self.init_model.initZ(qmean="random")

    def build_W(self):
        """ Build node W for the weights """
        if self.model_opts['sl_w']:
            self.init_model.initSW(qmean_S1='random')
        else:
            self.init_model.initW(qmean='random')

    def build_Tau(self):
        # TODO sort out how to choose where to use Tau
        self.init_model.initTau(self.data_opts['samples_groups'])

    def build_AlphaZ(self):
        """ Build node AlphaZ for the ARD prior on the factors """

        # ARD prior per sample group
        self.init_model.initAlphaZ(self.data_opts['samples_groups'])

    def build_AlphaW(self):
        """ Build node AlphaW for the ARD prior on the weights"""

        # ARD prior per factor and feature_group (view)
        self.init_model.initAlphaW()

    def build_ThetaZ(self):
        # TODO use mixed theta node instead when fixed in update and init -> should be ok then for intercept
        """ Build node ThetaZ for the Spike and Slab prior on the factors """

        # Initialise hyperparameters for the ThetaZ prior
        initTheta_a = 1.
        initTheta_b = 1.

        self.init_model.initThetaZ(self.data_opts['samples_groups'], qa=initTheta_a, qb=initTheta_b)

    def build_ThetaW(self):
        """ Build node ThetaW for the Spike and Slab prior on the weights """

        # Initialise hyperparameters for the ThetaW prior
        initTheta_a = 1.
        initTheta_b = 1.

        self.init_model.initThetaW(qa=initTheta_a, qb=initTheta_b)

    def createMarkovBlankets(self):
        """ Define the markov blankets """

        # Fetch all nodes
        nodes = self.init_model.getNodes()

        # Define basic connections
        nodes['Y'].addMarkovBlanket(Z=nodes["Z"], W=nodes["W"], Tau=nodes["Tau"])
        nodes['Z'].addMarkovBlanket(Y=nodes["Y"], W=nodes["W"], Tau=nodes["Tau"])
        nodes['W'].addMarkovBlanket(Y=nodes["Y"], Z=nodes["Z"], Tau=nodes["Tau"])
        nodes['Tau'].addMarkovBlanket(Y=nodes["Y"], W=nodes["W"], Z=nodes["Z"])

        # Add ThetaZ in the markov blanket of Z and viceversa if Spike and Slab prior on Z
        if self.model_opts['sl_z']:
            nodes['Z'].addMarkovBlanket(ThetaZ=nodes["ThetaZ"])
            nodes["ThetaZ"].addMarkovBlanket(Z=nodes["Z"])

        # Add ThetaW in the markov blanket of W and viceversa if Spike and Slab prior on Z
        if self.model_opts['sl_w']:
            nodes['W'].addMarkovBlanket(ThetaW=nodes["ThetaW"])
            nodes["ThetaW"].addMarkovBlanket(W=nodes["W"])

        # Add AlphaZ in the markov blanket of Z and viceversa if ARD prior on Z
        if self.model_opts['ard_z']:
            nodes['AlphaZ'].addMarkovBlanket(Z=nodes['Z'])
            nodes['Z'].addMarkovBlanket(AlphaZ=nodes['AlphaZ'])

        # Add AlphaW in the markov blanket of W and viceversa if ARD prior on W
        if self.model_opts['ard_w']:
            nodes['AlphaW'].addMarkovBlanket(W=nodes['W'])
            nodes['W'].addMarkovBlanket(AlphaW=nodes['AlphaW'])

    def createSchedule(self):
        """ Define the schedule of updates """

        # TO-DO:
        # - RICARD: I THINK THE SCHEDULE OF UPDATES SHOULD NOT BE INSIDE BUILD_MODEL
        # - ALLOW SCHEDULE TO BE PROVIDED AS TRAIN_OPTIONS
        # - IF PROVIDED, SO A SANITY CHECKS THAT THE CORRECT NODES CAN BE FOUND AND THERE ARE NO DUPLICATED

        # Define basic schedule of updates
        schedule = ['Y', 'W', 'Z', 'Tau']
        # schedule = ['Y', 'Z', 'W', 'Tau']

        # Insert ThetaW after W if Spike and Slab prior on W
        if self.model_opts['sl_w']:
            ix = schedule.index("W")
            schedule.insert(ix+1, 'ThetaW')

        # Insert ThetaZ after Z if Spike and Slab prior on Z
        if self.model_opts['sl_z']:
            ix = schedule.index("Z")
            schedule.insert(ix+1, 'ThetaZ')

        # Insert AlphaW after W if ARD prior on W
        if self.model_opts['ard_w']:
            ix = schedule.index("W")
            schedule.insert(ix+1, 'AlphaW')

        # Insert AlphaZ after Z if ARD prior on Z
        if self.model_opts['ard_z']:
            ix = schedule.index("Z")
            schedule.insert(ix+1, 'AlphaZ')

        self.schedule = schedule

    def createBayesNet(self):
        """ Method to create the BayesNet class """
        self.net = BayesNet(dim=self.dim, nodes=self.init_model.getNodes())
