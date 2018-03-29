"""
Module to build and initialise a bioFAM model
"""

import scipy as s
from sys import path
from time import time,sleep
import pandas as pd
import numpy as np
#from joblib import Parallel, delayed

from biofam.core.BayesNet import *
from init_model import initNewModel
from biofam.build_model.utils import *


def build_model(model_opts, data=None):
    """Method to build a bioFAM model"""

    print("\n")
    print("#"*24)
    print("## Building the model ##")
    print("#"*24)
    print("\n")
    sleep(1)

    # Define dimensionalities
    if data is None:
        M = model_opts['M']
        N = model_opts['N']
        D = model_opts['D']
        data = [np.ones([N, D[m]]) * np.nan for m in range(M)]
    else:
        M = len(data)
        N = data[0].shape[0]
        D = s.asarray([ data[m].shape[1] for m in range(M) ])

    K = model_opts["K"]
    dim = {'M':M, 'N':N, 'D':D, 'K':K }

    ###########################
    ## Do some sanity checks ##
    ###########################

    # If learnIntercept is True, add one extra factor as a covariate with constant 1s
    if model_opts["learnIntercept"]:
        dim["K"] += 1
        if model_opts['covariates'] is not None:
            model_opts['covariates'] = s.insert(model_opts['covariates'], obj=0, values=1, axis=1)
            model_opts['scale_covariates'].insert(0,False)
        else:
            model_opts['covariates'] = s.ones((dim["N"],1))
            model_opts['scale_covariates'] = [False]

    #####################################
    ## Define and initialise the nodes ##
    #####################################

    # Initialise the model
    init = initNewModel(dim, data, model_opts["likelihood"])

    if model_opts['transpose_sparsity']:

        # Initialise latent variables
        priorZ_mean_S0 = 0.; priorZ_meanS1 = 0.; priorZ_varS0 = 1.; priorZ_varS1 = 1.; priorZ_theta = 1.;
        initZ_meanS0 = 0.; initZ_meanS1 = 0.; initZ_varS0 = 1.; initZ_varS1 = 1; initZ_theta = 1.;
        initZ_qEZ_S0 = 0.; initZ_qEZ_S1 = 0.; initZ_qES = 1.

        init.initSZ(pmean_S0=priorZ_mean_S0, pmean_S1=priorZ_meanS1, pvar_S0=priorZ_varS0, pvar_S1=priorZ_varS1,
                    ptheta=priorZ_theta,
                    qmean_S0=initZ_meanS0, qmean_S1=initZ_meanS1, qvar_S0=initZ_varS0, qvar_S1=initZ_varS1,
                    qtheta=initZ_theta,
                    qEZ_S0=initZ_qEZ_S0, qEZ_S1=initZ_qEZ_S1, qES=initZ_qES)

        pmean = 0.;
        pvar = 1.;
        qmean = "random";
        qvar = 1.

        init.initW(pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, covariates=model_opts['covariates'],
                   scale_covariates=model_opts['scale_covariates'])

        # Initialise sparsity on latent variables
        learnTheta_ix = np.ones(K)        
        initTheta_a = 1.
        initTheta_b = 1.
        initTheta_E = 1.
        priorTheta_a = 1.
        priorTheta_b = 1.
        # TO-DOOOOOOOOOO
        # if model_opts["learnIntercept"]:
        #     learnTheta[:,0] = 0. # Remove sparsity from the weight vector that will capture the feature-wise means
        
        if model_opts["learnIntercept"]:
            for ix in learnTheta_ix:
                ix[0] = 0
        init.initThetaMixedZ_k(learnTheta_ix, pa=priorTheta_a, pb=priorTheta_b, qa=initTheta_a, qb=initTheta_b,
                               qE=initTheta_E)

    else:

        # Initialise latent variables
        pmean = 0.; pvar = 1.; qmean = "random"; qvar = 1.
        init.initZ(pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, covariates=model_opts['covariates'], scale_covariates=model_opts['scale_covariates'])

        # Initialise weights
        priorSW_mean_S0=0.; priorSW_meanS1=0.; priorSW_varS0=1.; priorSW_varS1=1.; priorSW_theta=1.;
        initSW_meanS0=0.; initSW_meanS1=0.; initSW_varS0=1.; initSW_varS1=1; initSW_theta=1.;
        initSW_qEW_S0=0.; initSW_qEW_S1=0.; initSW_qES=1.;

        # TO-DOOOOOOOOOO: fIX LEARN INTERCEPT
        # if model_opts["learnIntercept"]:
        # for m in range(M):
        #     if model_opts["likelihood"][m]=="gaussian":
        #         initSW_meanS1[m][:,0] = data[m].mean(axis=0)
        #         initSW_varS1[m][:,0] = 1e-5
        #         initSW_theta[m][:,0] = 1.

        init.initSW(pmean_S0=priorSW_mean_S0, pmean_S1=priorSW_meanS1, pvar_S0=priorSW_varS0, pvar_S1=priorSW_varS1, ptheta=priorSW_theta,
            qmean_S0=initSW_meanS0, qmean_S1=initSW_meanS1, qvar_S0=initSW_varS0, qvar_S1=initSW_varS1, qtheta=initSW_theta,
            qEW_S0=initSW_qEW_S0, qEW_S1=initSW_qEW_S1, qES=initSW_qES)
        
        # Initialise sparsity on weights
        learnTheta = [ s.ones((D[m], K)) for m in range(M) ]
        priorTheta_a = 1.
        priorTheta_b = 1.
        initTheta_a = 1.
        initTheta_b = 1.
        initTheta_E = 1.
        # TO-DOOOOOOOOOO
        # if model_opts["learnIntercept"]:
        #     for m in range(M):
        #     learnTheta[m][:,0] = 0. # Remove sparsity from the weight vector that will capture the feature-wise means
        learnTheta_ix = [np.ones(K)]*M
        if model_opts["learnIntercept"]:
            for ix in learnTheta_ix:
                ix[0] = 0
        init.initThetaMixedW_mk(learnTheta_ix, pa=priorTheta_a, pb=priorTheta_b, qa=initTheta_a,  qb=initTheta_b, qE=initTheta_E)


    # Initialise ARD on weights
    pa=1e-14; pb=1e-14; qa=1.; qb=1.; qE=1.
    init.initAlphaW_mk(pa=pa, pb=pb, qa=qa, qb=qb)

    # Initialise ARD on latent variables
    if model_opts['transpose_sparsity']:
        qa = 1.; qb = 1.
        init.initAlphaZ_k(pa=pa, pb=pb, qa=qa, qb=qb)

    # Initialise precision of noise
    pa=1e-14; pb=1e-14; qa=1.; qb=1.; qE=1.
    init.initTau(pa=pa, pb=pb, qa=qa, qb=qb, transposed=model_opts["transpose_noise"])


    # Observed data
    init.initY(transpose_noise=model_opts["transpose_noise"])

    ############################################
    ## Define the markov blanket of each node ##
    ############################################

    nodes = init.getNodes()
    if model_opts['transpose_sparsity']:
        nodes["SZ"].addMarkovBlanket(W=nodes["W"], Tau=nodes["Tau"], AlphaZ=nodes["AlphaZ"], Y=nodes["Y"], ThetaZ=nodes["ThetaZ"])
        nodes["AlphaW"].addMarkovBlanket(W=nodes["W"])
        nodes["AlphaZ"].addMarkovBlanket(SZ=nodes["SZ"])
        nodes["W"].addMarkovBlanket(SZ=nodes["SZ"], Tau=nodes["Tau"], AlphaW=nodes["AlphaW"], Y=nodes["Y"])
        nodes["Y"].addMarkovBlanket(SZ=nodes["SZ"], W=nodes["W"], Tau=nodes["Tau"])
        nodes["Tau"].addMarkovBlanket(SZ=nodes["SZ"], W=nodes["W"], Y=nodes["Y"])
        nodes["ThetaZ"].addMarkovBlanket(SZ=nodes["SZ"])
    else:
        nodes["Z"].addMarkovBlanket(SW=nodes["SW"], Tau=nodes["Tau"], Y=nodes["Y"])
        nodes["ThetaW"].addMarkovBlanket(SW=nodes["SW"])
        nodes["AlphaW"].addMarkovBlanket(SW=nodes["SW"])
        nodes["SW"].addMarkovBlanket(Z=nodes["Z"], Tau=nodes["Tau"], Alpha=nodes["AlphaW"], Y=nodes["Y"], Theta=nodes["ThetaW"])
        nodes["Y"].addMarkovBlanket(Z=nodes["Z"], SW=nodes["SW"], Tau=nodes["Tau"])
        nodes["Tau"].addMarkovBlanket(Z=nodes["Z"], SW=nodes["SW"], Y=nodes["Y"])

    #################################
    ## Initialise Bayesian Network ##
    #################################

    net = BayesNet(dim=dim, nodes=init.getNodes())

    return net
