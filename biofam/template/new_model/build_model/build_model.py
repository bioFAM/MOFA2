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
    # TODO : enable using covariance matrix for prior on W instead of spike and slab

    print ("\n")
    print ("#"*24)
    print ("## Building the model ##")
    print ("#"*24)
    print ("\n")
    sleep(1)

    # Define dimensionalities

    # TODO : not the right entry point for X and cluster files if transpose (multiple files)

    if data is None:

        M = model_opts['M']
        N = model_opts['N']
        D = model_opts['D']

        if model_opts["transpose"]:
            X = [None] * M
            sigma_clust = [None] * M
            view_has_covariance_prior = [True] * M

            for m in range(M):
                if model_opts['positions_samples_file'][m] is not None:
                    X[m] = np.loadtxt(model_opts['positions_samples_file'][m])
                else:
                    X[m] = s.random.norm(0, 1, [D[m], 2])

        else:
            if model_opts['positions_samples_file'] is not None:
                X = np.loadtxt(model_opts['positions_samples_file'])
            else:
                X = s.random.norm(0, 1, [N, 2])
            sigma_clust = None

        data = [np.ones([N, D[m]]) * np.nan for m in range(M)]

    else:

        M = len(data)
        N = data[0].shape[0]
        D = s.asarray([data[m].shape[1] for m in range(M)])

        if model_opts["transpose"]:
            X = [None] * M
            sigma_clust = [None] * M
            view_has_covariance_prior = [None] * M

            # TODO change that, super dirty
            for m in range(M):
                if model_opts['positions_samples_file'][m] is not None:
                    view_has_covariance_prior[m] = True
                    try:
                        X[m] = np.loadtxt(model_opts['positions_samples_file'], delimiter=',')
                    except:
                        X[m] = np.loadtxt(model_opts['positions_samples_file'], delimiter=' ')
                else:
                    view_has_covariance_prior[m] = False

                # load sigma cluster if among arguments
                if model_opts['sigma_cluster_file'][m] is not None:
                    sigma_clust[m] = np.loadtxt(model_opts['sigma_cluster_file'])
                else:
                    pass

            if model_opts['permute_samples'] == 1:
                for m in range(M):
                    perm = np.random.permutation(D[m])
                    X[m]= X[m][perm, :]

        else:
            # TODO change that, super dirty
            for m in range(M):
                if model_opts['positions_samples_file'] is not None:
                    try:
                        X = np.loadtxt(model_opts['positions_samples_file'], delimiter=',')
                    except:
                        X = np.loadtxt(model_opts['positions_samples_file'], delimiter=' ')
                else:
                    pass

                # load sigma cluster if among arguments
                if model_opts['sigma_cluster_file'] is not None:
                    sigma_clust = np.loadtxt(model_opts['sigma_cluster_file'])
                else:
                    pass

            if model_opts['permute_samples'] == 1:
                perm = np.random.permutation(N)
                X = X[perm, :]

    K = model_opts["K"]

    if 'spatialFact' in model_opts:
        n_diag = model_opts['spatialFact'] * K
    else:
        n_diag = 0

    dim = {'M': M, 'N': N, 'D': D, 'K': K}
    print(dim)

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

    # Initialise latent variables
    if model_opts['transpose']:
        priorZ_mean_T0 = 0.
        priorZ_meanT1 = 0.
        priorZ_varT0 = 1.
        priorZ_varT1 = 1.
        priorZ_theta = 1. # Question : fix s to 1 in transpose model ?
        initZ_meanT0 = 0.
        initZ_meanT1 = 0.
        initZ_varT0 = 1.
        initZ_varT1 = 1
        initZ_theta = 1.
        initZ_qEZ_T0 = 0.
        initZ_qEZ_T1 = 0.
        initZ_qET = 1.
        init.initTZ(pmean_T0=priorZ_mean_T0, pmean_T1=priorZ_meanT1, pvar_T0=priorZ_varT0, pvar_T1=priorZ_varT1,
                    ptheta=priorZ_theta,
                    qmean_T0=initZ_meanT0, qmean_T1=initZ_meanT1, qvar_T0=initZ_varT0, qvar_T1=initZ_varT1,
                    qtheta=initZ_theta,
                    qEZ_T0=initZ_qEZ_T0, qEZ_T1=initZ_qEZ_T1, qET=initZ_qET)
    else:
        pmean = 0.;
        pcov = 1.;
        qmean = "random";
        qvar = 1.
        init.initZ(pmean=pmean, pcov=pcov, qmean=qmean, qvar=qvar, covariates=model_opts['covariates'],
                   scale_covariates=model_opts['scale_covariates'])

    #Initialise weights
    if model_opts['transpose']:
        pmean = 0.;
        pvar = 1.;
        qmean = "random";
        qvar = 1.
        init.initW(pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, covariates=model_opts['covariates'],
                   scale_covariates=model_opts['scale_covariates'])
    else:
        priorW_mean_S0 = 0.;
        priorW_meanS1 = 0.;
        priorW_varS0 = 1.;
        priorW_varS1 = 1.;
        priorW_theta = 1.;
        initW_meanS0 = 0.;
        initW_meanS1 = 0.;
        initW_varS0 = 1.;
        initW_varS1 = 1;
        initW_theta = 1.;
        initW_qEW_S0 = 0.;
        initW_qEW_S1 = 0.;
        initW_qES = 1.;

        # TO-DOOOOOOOOOO: fIX LEARN INTERCEPT
        # if model_opts["learnIntercept"]:
        # for m in range(M):
        #     if model_opts["likelihood"][m]=="gaussian":
        #         initSW_meanS1[m][:,0] = data[m].mean(axis=0)
        #         initSW_varS1[m][:,0] = 1e-5
        #         initSW_theta[m][:,0] = 1.

        init.initSW(pmean_S0=priorW_mean_S0, pmean_S1=priorW_meanS1, pvar_S0=priorW_varS0, pvar_S1=priorW_varS1, ptheta=priorW_theta,
            qmean_S0=initW_meanS0, qmean_S1=initW_meanS1, qvar_S0=initW_varS0, qvar_S1=initW_varS1, qtheta=initW_theta,
            qEW_S0=initW_qEW_S0, qEW_S1=initW_qEW_S1, qES=initW_qES)


    # Initialise ARD or covariance prior structure on W for each view

    if model_opts['transpose']:
        params = [None] * M
        for m in range(M):
            if view_has_covariance_prior[m]:
                params[m]={'X':X[m],'sigma_clust':sigma_clust[m],'n_diag':n_diag}
            else:
                params[m]={'pa':1e-14, 'pb':1e-14, 'qa':1., 'qb':1., "qE":1.}

        init.initMixedSigmaAlphaW_mk(view_has_covariance_prior,params)
    else:
        # TODO do sth here for siulations
        pa = 1e-14; pb = 1e-14; qa = 1.; qb = 1.; qE = 1.
        init.initAlphaW_mk(pa=pa, pb=pb, qa=qa, qb=qb)


    # Initialise ARD or covariance prior structure on Z

    if model_opts['transpose']:
        pa = 1e-14; pb = 1e-14; qa = 1.; qb = 1.; qE = 1.
        init.initAlphaZ_k(pa=pa, pb=pb, qa=qa, qb=qb)
    else:
        if model_opts['covariance_samples'] is not None:
            # TODO add a if statement to check if there is a sigma_clust argument to see if blockSigma is needed
            if sigma_clust is None:
                init.initSigmaZ_k(X, n_diag=n_diag)
            else:
                init.initSigmaBlockZ_k(X, clust=sigma_clust, n_diag=n_diag)
        else:
            pa = 1e-14; pb = 1e-14; qa = 1.; qb = 1.; qE = 1.
            init.initAlphaZ_k(pa=pa, pb=pb, qa=qa, qb=qb)


    # Initialise precision of noise
    # TODO do sth here for siulations
    pa=1e-14; pb=1e-14; qa=1.; qb=1.; qE=1.
    init.initTau(pa=pa, pb=pb, qa=qa, qb=qb)

    if model_opts["transpose"]:
        # Initialise sparsity on the factors
        learnTheta = s.ones((1, K))
        priorTheta_a = 1.
        priorTheta_b = 1.
        initTheta_a = 1.
        initTheta_b = 1.
        initTheta_E = 1.
        # TO-DOOOOOOOOOO
        # if model_opts["learnIntercept"]:
        #     learnTheta[:,0] = 0. # Remove sparsity from the weight vector that will capture the feature-wise means
        learnTheta_ix = np.ones(K)
        if model_opts["learnIntercept"]:
            for ix in learnTheta_ix:
                ix[0] = 0
        init.initThetaMixedZ_k(learnTheta_ix, pa=priorTheta_a, pb=priorTheta_b, qa=initTheta_a, qb=initTheta_b,
                               qE=initTheta_E)

    else:
        # Initialise sparsity on the weights
        learnTheta = [s.ones((D[m], K)) for m in xrange(M)]
        priorTheta_a = 1.
        priorTheta_b = 1.
        initTheta_a = 1.
        initTheta_b = 1.
        initTheta_E = 1.
        # TO-DOOOOOOOOOO
        # if model_opts["learnIntercept"]:
        #     for m in range(M):
        #     learnTheta[m][:,0] = 0. # Remove sparsity from the weight vector that will capture the feature-wise means
        learnTheta_ix = [np.ones(K)] * M
        if model_opts["learnIntercept"]:
            for ix in learnTheta_ix:
                ix[0] = 0
        init.initThetaMixedW_mk(learnTheta_ix, pa=priorTheta_a, pb=priorTheta_b, qa=initTheta_a, qb=initTheta_b,
                                qE=initTheta_E)

    # Observed data
    init.initY()

    ############################################
    ## Define the markov blanket of each node ##
    ############################################

    nodes = init.getNodes()
    if model_opts['transpose']:
        nodes["ThetaZ"].addMarkovBlanket(TZ=nodes["TZ"])
        nodes["AlphaZ"].addMarkovBlanket(TZ=nodes["TZ"])
        nodes["TZ"].addMarkovBlanket(AlphaZ=nodes["AlphaZ"], ThetaZ=nodes["ThetaZ"], Y=nodes["Y"], W=nodes["W"],
                                     Tau=nodes["Tau"])

        if model_opts["covariance_samples"]:
            nodes["SigmaAlphaW"].addMarkovBlanket(W=nodes["W"])
            nodes["W"].addMarkovBlanket(SigmaAlphaW=nodes["SigmaAlphaW"])
        else:
            nodes["AlphaW"].addMarkovBlanket(W=nodes["W"])
            nodes["W"].addMarkovBlanket(AlphaW=nodes["AlphaW"])
        nodes["W"].addMarkovBlanket(Y=nodes["Y"], TZ=nodes["TZ"], Tau=nodes["Tau"])

        nodes["Y"].addMarkovBlanket(TZ=nodes["TZ"], W=nodes["W"], Tau=nodes["Tau"])
        nodes["Tau"].addMarkovBlanket(TZ=nodes["TZ"], W=nodes["W"], Y=nodes["Y"])

    else:
        if model_opts["covariance_samples"]:
            nodes["SigmaZ"].addMarkovBlanket(Z=nodes["Z"])
            nodes["Z"].addMarkovBlanket(SigmaZ=nodes["SigmaZ"])
        else:
            nodes["AlphaZ"].addMarkovBlanket(Z=nodes["Z"])
            nodes["Z"].addMarkovBlanket(AlphaZ=nodes["AlphaZ"])
        nodes["Z"].addMarkovBlanket(Y=nodes["Y"], SW=nodes["SW"], Tau=nodes["Tau"])

        nodes["ThetaW"].addMarkovBlanket(SW=nodes["SW"])
        nodes["AlphaW"].addMarkovBlanket(SW=nodes["SW"])
        nodes["SW"].addMarkovBlanket(AlphaW=nodes["AlphaW"], ThetaW=nodes["ThetaW"], Y=nodes["Y"], Z=nodes["Z"],
                                     Tau=nodes["Tau"])
        nodes["Y"].addMarkovBlanket(Z=nodes["Z"], SW=nodes["SW"], Tau=nodes["Tau"])
        nodes["Tau"].addMarkovBlanket(Z=nodes["Z"], SW=nodes["SW"], Y=nodes["Y"])






    #################################
    ## Initialise Bayesian Network ##
    #################################

    net = BayesNet(dim=dim, nodes=init.getNodes())

    return net
