"""
Module to build and initialise a bioFAM model
"""

from sys import path
from time import time,sleep
import numpy as np
import scipy as s
#from joblib import Parallel, delayed

from biofam.core.BayesNet import *
from biofam.build_model.init_model import initModel
from biofam.build_model.utils import *

# # QUESTION: do we need this ?
class buildModel(object):
    def __init__(self):
        self.net = None

# TODO could also split in  different model: pne basic parent with no ARD, one below which would be groups biofam where you choose what ARD to use
# TODO create schedule here too ?
class buildBiofam(buildModel):

    def  __init__(self, data_opts, model_opts):
        # defining model's dimensionalities
        self.data = data_opts['data']
        M = len(self.data)
        N = self.data[0].shape[0]
        D = s.asarray([self.data[m].shape[1] for m in range(M)])
        K = model_opts["factors"]
        self.dim = {'M': M, 'N': N, 'D': D, 'K': K}

        # update for learning learn_intercept
        if model_opts["learn_intercept"]:
            self.dim["K"] += 1
            if 'covariates_files' in model_opts:
                model_opts['covariatesFiles'] = s.insert(model_opts['covariatesFiles'], obj=0, values=1, axis=1)
                model_opts['scale_covariates'].insert(0,False)
            else:
                model_opts['covariates_files'] = s.ones((dim["N"],1))
                model_opts['scale_covariates'] = [False]

        # save options in the object
        self.data_opts = data_opts
        self.model_opts = model_opts

        # create an instance of initModel and start buidling
        self.init_model = initModel(self.dim, self.data, self.model_opts["likelihoods"])
        self.net = None
        self.build_all()

    def build_all(self):
        self.build_Y()
        self.build_Z()
        self.build_W()
        self.build_Tau()

        # define ARDs
        if self.model_opts['ard_z']:
            self.build_AlphaZ()
        if self.model_opts['ard_w']:
            self.build_AlphaW()

        # define thetas
        if self.model_opts['sl_z']:
            self.build_ThetaZ()
        if self.model_opts['sl_w']:
            self.build_ThetaW()

        self.createMarkovBlankets()
        self.createSchedule()
        self.net = BayesNet(dim=self.dim, nodes=self.init_model.getNodes())

    def build_Y(self):
        # TODO enable a noise on cells instead of features with some option
        self.init_model.initY(sample_wise_noise=False)

    def build_Z(self):
        # Initialise latent variables
        # TODO change to 'if sl_Z in model_opts ' ? would enable to not have to
        # add all the possible options to model_opt
        if self.model_opts['sl_z']:
            self.init_model.initSZ()
        else:
            # TODO change Z node so that we dont use a multivariate prior when no covariance structure
            self.init_model.initZ()

    def build_W(self):
        #Initialise weights
        if self.model_opts['sl_w']:
            self.init_model.initSW()
        else:
            # TODO change Z node so that we dont use a multivariate prior when no covariance structure
            # TODO could also make sure that SZ and SW can have a Sigma covariance prior
            self.init_model.initW()

    def build_Tau(self):
        # TODO sort out dimensionality of Tau
        self.init_model.initTau(transposed=False)

    def build_AlphaZ(self):
        # cases to distinguish are whether there are groups or not and whether we should account for them
        if self.data_opts['sample_groups'] is not None:
            self.init_model.initAlphaZ_groups(self.data_opts['sample_groups'])
        else:
            self.init_model.initAlphaZ_k()

    def build_AlphaW(self):
        self.init_model.initAlphaW_mk()

    def build_ThetaZ(self):
        initTheta_a = 1.
        initTheta_b = 0.001 #1.
        learnTheta_ix = np.ones(self.dim['K'])

        # TODO: this for loop cannot possibly work change
        if self.model_opts["learn_intercept"]:
            for ix in learnTheta_ix:
                ix[0] = 0

        self.init_model.initThetaMixedZ_k(learnTheta_ix, qa=initTheta_a, qb=initTheta_b)

    def build_ThetaW(self):
        # Initialise sparsity on the weights
        initTheta_a = 1.
        initTheta_b = 1.
        learnTheta_ix = [np.ones(self.dim['K'])] * self.dim['M']

        # TODO: this for loop cannot possibly work change
        if self.model_opts["learn_intercept"]:
            for ix in learnTheta_ix:
                ix[0] = 0

        self.init_model.initThetaMixedW_mk(learnTheta_ix, qa=initTheta_a, qb=initTheta_b)

    def createMarkovBlankets(self):
        nodes = self.init_model.getNodes()

        # basic connections
        nodes['Y'].addMarkovBlanket(Z=nodes["Z"], W=nodes["W"], Tau=nodes["Tau"])
        nodes['Z'].addMarkovBlanket(Y=nodes["Y"], W=nodes["W"], Tau=nodes["Tau"])
        nodes['W'].addMarkovBlanket(Y=nodes["Y"], Z=nodes["Z"], Tau=nodes["Tau"])
        nodes['Tau'].addMarkovBlanket(Y=nodes["Y"], W=nodes["W"], Z=nodes["Z"])

        # adding theta nodes if spike and slab
        if self.model_opts['sl_z']:
            nodes['Z'].addMarkovBlanket(ThetaZ=nodes["ThetaZ"])
            nodes["ThetaZ"].addMarkovBlanket(Z=nodes["Z"])
        if self.model_opts['sl_w']:
            nodes['W'].addMarkovBlanket(ThetaW=nodes["ThetaW"])
            nodes["ThetaW"].addMarkovBlanket(W=nodes["W"])

        # adding alpha nodes if ARD
        if self.model_opts['ard_z']:
            nodes['AlphaZ'].addMarkovBlanket(Z=nodes['Z'])
            nodes['Z'].addMarkovBlanket(AlphaZ=nodes['AlphaZ'])
        if self.model_opts['ard_w']:
            nodes['AlphaW'].addMarkovBlanket(W=nodes['W'])
            nodes['W'].addMarkovBlanket(AlphaW=nodes['AlphaW'])

    def createSchedule(self):
        schedule = np.array(['Y', 'W', 'Z', 'Tau'])

        # insert thetaW after W
        if self.model_opts['sl_w']:
            ix = self.find_node(schedule, 'W')[0][0]
            np.insert(schedule, ix + 1, 'ThetaW')

        if self.model_opts['sl_z']:
            ix = self.find_node(schedule, 'Z')[0][0]
            np.insert(schedule, ix + 1, 'ThetaZ')

        if self.model_opts['ard_w']:
            ix = self.find_node(schedule, 'W')[0][0]
            np.insert(schedule, ix + 1, 'AlphaW')

        if self.model_opts['ard_z']:
            ix = self.find_node(schedule, 'Z')[0][0]
            np.insert(schedule, ix + 1, 'AlphaZ')

        self.schedule = schedule

    def find_node(self, arr, nd):
        truth_table = (arr == nd)
        indices = np.where(truth_table)
        if len(indices) == 0:
            print('Warning, node', nd, 'not found')
        if len(indices) >1 :
            print('Warning, node', nd, 'found on multiple locations')
        return indices


class buildSpatialBiofam(buildBiofam):
    def __init__(self):
        pass

    def build_AlphaW(self):
        # TODO adapt and same for alphaZ
        if dataX is not None:
            params = [None] * M
            for m in range(M):
                if view_has_covariance_prior[m]:
                    params[m]={'X':dataX[m],'sigma_clust':dataClust[m],'n_diag':n_diag}
                else:
                    params[m]={'pa': pa, 'pb':pb, 'qa':1., 'qb':1.}
            init.initMixedSigmaAlphaW_mk(view_has_covariance_prior,params)

        if dataX is not None:
            # TODO add a if statement to check if there is a sigma_clust argument to see if blockSigma is needed
            if dataClust is None:
                init.initSigmaZ_k(dataX, n_diag=n_diag)
            else:
                init.initSigmaBlockZ_k(dataX, clust=dataClust, n_diag=n_diag)
    def createMarkovBlankets(self):
        # TODO adapt this code
        if dataX is not None:
            nodes["SigmaAlphaW"].addMarkovBlanket(W=nodes["W"])
            nodes["W"].addMarkovBlanket(SigmaAlphaW=nodes["SigmaAlphaW"])

        if dataX is not None:
            nodes["SigmaZ"].addMarkovBlanket(Z=nodes["Z"])
            nodes["Z"].addMarkovBlanket(SigmaZ=nodes["SigmaZ"])


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

        self.super().build_all()
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
