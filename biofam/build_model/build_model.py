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
        self.init_model.initY(noise_on=self.model_opts['noise_on'])

    def build_Z(self):
        """ Build node Z for the factors or latent variables """

        # if learning the intercept, add a covariate of ones to Z
        if self.model_opts["learn_intercept"]:
            self.model_opts['covariates'] = s.ones((self.dim["N"],1))
            self.model_opts['scale_covariates'] = [False]

        if self.model_opts['sl_z']:
            self.init_model.initSZ()
        else:
            # TODO change Z node so that we dont use a multivariate prior when no covariance structure
            self.init_model.initZ()

    def build_W(self):
        """ Build node W for the weights """
        if self.model_opts['sl_w']:
            self.init_model.initSW()
        else:
            # TODO change Z node so that we dont use a multivariate prior when no covariance structure
            # TODO could also make sure that SZ and SW can have a Sigma covariance prior
            self.init_model.initW()

    def build_Tau(self):
        # TODO sort out how to choose where to use Tau
        self.init_model.initTau(on=self.model_opts['noise_on'])

    def build_AlphaZ(self):
        """ Build node AlphaZ for the ARD prior on the factors """

        # ARD prior per sample group
        if self.data_opts['sample_groups'] is not None:
            # TODO check whats going on here
            self.init_model.initAlphaZ_groups(self.data_opts['sample_groups'])

        # ARD prior per factor
        else:
            self.init_model.initAlphaZ()

    def build_AlphaW(self):
        """ Build node AlphaW for the ARD prior on the weights"""

        # ARD prior per factor and feature_group (view)
        # RICARD: SYMMETRISE AS IN ALPHAZ??
        self.init_model.initAlphaW()

    def build_ThetaZ(self):
        """ Build node ThetaZ for the Spike and Slab prior on the factors """

        # Initialise hyperparameters for the ThetaZ prior
        initTheta_a = 1.
        initTheta_b = 1.#0.001  #0.001 #1.

        # Specify for which factors to learn ThetaZ
        learnTheta_ix = np.ones(self.dim['K'])

        # Do not learn ThetaZ for the intercept factor
        if self.model_opts["learn_intercept"]:
            learnTheta_ix[0] = 0

        self.init_model.initThetaZ_Mixed(learnTheta_ix, qa=initTheta_a, qb=initTheta_b)

    def build_ThetaW(self):
        """ Build node ThetaW for the Spike and Slab prior on the weights """

        # Initialise hyperparameters for the ThetaW prior
        initTheta_a = 1.
        initTheta_b = 1.#.001  #0.001 #1.
        learnTheta_ix = [np.ones(self.dim['K'])] * self.dim['M']

        # TODO: this for loop cannot possibly work change
        if self.model_opts["learn_intercept"]:
            for ix in learnTheta_ix:
                ix[0] = 0
        self.init_model.initThetaW_Mixed(learnTheta_ix, qa=initTheta_a, qb=initTheta_b)

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
        # schedule = ['Y', 'Z', 'W', 'Tau']
        schedule = ['Y', 'W', 'Z', 'Tau']

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













## IGNORE BELOW, UNDER CONSTRUCTION ###


class buildSpatialBiofam(buildBiofam):

    def __init__(self, data_opts, model_opts):
        assert 'data_x' in data_opts, 'positions not found in data options'
        super(buildSpatialBiofam, self).__init__(data_opts, model_opts)

    def build_all(self):
        # build the nodes which are missing in build_all
        self.build_Sigma()
        super(buildSimulationBiofam, self).build_all()

    def build_Sigma(self):
        if self.model_opts['cov_on'] == 'samples':
            # TODO add a if statement to check if there is a sigma_clust argument to see if blockSigma is needed
            if self.data_opts['dataClust'] is None:
                self.init_model.initSigmaZ(self.data_opts['data_x'])
            else:
                self.init_model.initSigmaBlockZ(self.data_opts['data_x'], clust=self.data_opts['dataClust'])
        else:
            params = [None] * M
            for m in range(M):
                if self.data_opts['view_has_covariance_prior'][m]:
                    params[m]={'X':self.data_opts['data_x'][m],
                    'sigma_clust':self.data_opts['dataClust'],'n_diag':0}
                else:
                    params[m]={'pa': 1e-14, 'pb':1e-14, 'qa':1., 'qb':1.}
            self.init_model.initMixedSigmaAlphaW(view_has_covariance_prior,params)

    def createMarkovBlankets(self):
        super(buildSpatialBiofam, self).createMarkovBlankets()
        nodes = self.init_model.getNodes()

        # create the markov blanket for the sigma nodes
        if self.model_opts['cov_on'] == 'samples':
            nodes["Sigma"].addMarkovBlanket(Z=nodes["Z"])
            nodes["Z"].addMarkovBlanket(Sigma=nodes["Sigma"])

        if self.model_opts['cov_on'] == 'features':
            nodes["Sigma"].addMarkovBlanket(W=nodes["W"])
            nodes["W"].addMarkovBlanket(Sigma=nodes["Sigma"])

    def createSchedule(self):
        super(buildSpatialBiofam, self).createSchedule()

        # add the sigma node at the right position in the schedule
        if self.model_opts['cov_on'] == 'samples':
            ix = self.find_node(self.schedule, 'Z')[0][0]
            np.insert(self.schedule, ix + 1, 'Sigma')

        if self.model_opts['cov_on'] == 'features':
            ix = self.find_node(self.schedule, 'W')[0][0]
            np.insert(self.schedule, ix + 1, 'Sigma')


class buildSimulationBiofam(buildBiofam):
    def __init__(self, model_opts):
        M = model_opts['M']
        N = model_opts['N']
        D = model_opts['D']

        self.init_model = initNewModel(dim, data, model_opts["likelihoods"])

    def build_all():
        # TODO add somewhere
        # if 'spatialFact' in model_opts:
        #     n_diag = (1-model_opts['spatialFact']) * K
        # else:
        #     n_diag = K
        # if model_opts["sample_X"]:
        #     if model_opts["transpose_sparsity"]:
        #        dataX = [s.random.normal(0, 1, [D[m], 2]) for m in range(M)]
        #        dataClust = [None] * M
        #        view_has_covariance_prior = [True] * M
        #     else:
        #        dataX = s.random.normal(0, 1, [N, 2])
        #        dataClust  = None
        # data = [np.ones([N, D[m]]) * np.nan for m in range(M)]
        # if model_opts["transpose_sparsity"]:
        #
        #    if dataX is None :
        #        view_has_covariance_prior = [None] * M
        #    else:
        #        view_has_covariance_prior = [dataX[m] is not None for m in range(M)]
        #
        #    #for test purpose
        #    if (dataX is None) and (model_opts["sample_X"]):
        #        dataX = [s.random.normal(0, 1, [D[m], 2]) for m in range(M)]
        #        dataClust = [None] * M
        #        view_has_covariance_prior = [True] * M
        #
        # #for test purpose
        # else:
        #    if (dataX is None) and (model_opts["sample_X"]):
        #        dataX = s.random.normal(0, 1, [N, 2])
        #        dataClust = None

        super(buildSimulationBiofam, self).build_all()
        self.build_Sigma()
        # build other stuff

    def build_Sigma():
        self.init_model.initSigma()
        self.init_model['Sigma'].addMarkovBlanket(self.init_model['Z'])

    def build_Tau():
        # for simulations we should change the parameters
        pass
    def build_ThetaZ():
        # TODO reimplement that so we can simulate different levels of sparsity
        if 'sparsity' in model_opts:
            assert 0. <= model_opts['sparsity'] <= 1., 'sparsty level must be between 0 and 1'
            priorTheta_a = 10.
            priorTheta_b = ((1 - model_opts['sparsity']) / model_opts['sparsity']) * priorTheta_a
        elif 'noise' in model_opts:
            priorTheta_a = 10.
            priorTheta_b = 10.

    def build_Tau(self):
        #TODO enable different levels of Noise    if 'noise' in model_opts:  # simulation
        pa = 20.; pb = pa * model_opts['noise']
        pa = 1e-14; pb = 1e-14

    def build_AlphaW(self):
        # also requires different prior hyperparameters
        pass
