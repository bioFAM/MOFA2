"""
Module to build and initialise a bioFAM model
"""

from sys import path
from time import time,sleep
import numpy as np
import scipy as s
#from joblib import Parallel, delayed

from biofam.core.BayesNet import *
from init_model import initNewModel
from biofam.build_model.utils import *

#TODO : remove the 2 TODO before giving to britta

def build_model(model_opts, data=None, dataX=None, dataClust=None, dataCovariates=None):
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

        if model_opts["sample_X"]:
            if model_opts["transpose_sparsity"]:
               dataX = [s.random.normal(0, 1, [D[m], 2]) for m in range(M)]
               dataClust = [None] * M
               view_has_covariance_prior = [True] * M
            else:
               dataX = s.random.normal(0, 1, [N, 2])
               dataClust  = None

        data = [np.ones([N, D[m]]) * np.nan for m in range(M)]

    else:

        M = len(data)
        N = data[0].shape[0]
        D = s.asarray([data[m].shape[1] for m in range(M)])

        if model_opts["transpose_sparsity"]:

           if dataX is None :
               view_has_covariance_prior = [None] * M
           else:
               view_has_covariance_prior = [dataX[m] is not None for m in range(M)]

           #for test purpose
           if (dataX is None) and (model_opts["sample_X"]):
               dataX = [s.random.normal(0, 1, [D[m], 2]) for m in range(M)]
               dataClust = [None] * M
               view_has_covariance_prior = [True] * M

        #for test purpose
        else:
           if (dataX is None) and (model_opts["sample_X"]):
               dataX = s.random.normal(0, 1, [N, 2])
               dataClust = None

    K = model_opts["K"]

    if 'spatialFact' in model_opts:
        n_diag = (1-model_opts['spatialFact']) * K
    else:
        n_diag = K

    dim = {'M': M, 'N': N, 'D': D, 'K': K}
    print(dim)


    ###########################
    ## Do some sanity checks ##
    ###########################

    #TODO : add covariates by passing dataCovariates in argument of build_model (loading before in entry_point)
    # If learnIntercept is True, add one extra factor as a covariate with constant 1s
    if model_opts["learnIntercept"]:
        dim["K"] += 1
        if model_opts['covariatesFiles'] is not None:
            model_opts['covariatesFiles'] = s.insert(model_opts['covariatesFiles'], obj=0, values=1, axis=1)
            model_opts['scale_covariates'].insert(0,False)
        else:
            model_opts['covariatesFiles'] = s.ones((dim["N"],1))
            model_opts['scale_covariates'] = [False]

    #####################################
    ## Define and initialise the nodes ##
    #####################################

    # Initialise the model
    init = initNewModel(dim, data, model_opts["likelihood"])

    # Initialise latent variables
    if model_opts['transpose_sparsity']:
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
        init.initSZ(pmean_T0=priorZ_mean_T0, pmean_T1=priorZ_meanT1, pvar_T0=priorZ_varT0, pvar_T1=priorZ_varT1,
                    ptheta=priorZ_theta,
                    qmean_T0=initZ_meanT0, qmean_T1=initZ_meanT1, qvar_T0=initZ_varT0, qvar_T1=initZ_varT1,
                    qtheta=initZ_theta)
    else:
        pmean = 0.;
        pcov = 1.;
        qmean = "random";
        qvar = 1.;
        precompute_pcovinv = dataX is None
        # Warning : precompute_pcovinv must be True if and only if "AlphaZ" and "SigmaZ" not in init.nodes
        # (currently, we chosed that if we have no dataX, we do not have "AlphaZ" since it is redundant with "SigmaAlphaW"
        # and if we have dataX, we have necessarily "SigmaAlphaZ")
        init.initZ(pmean=pmean, pcov=pcov, qmean=qmean, qvar=qvar,  covariates=model_opts['covariatesFiles'], #covariates=dataCovariates,
                   scale_covariates=model_opts['scale_covariates'], precompute_pcovinv=precompute_pcovinv)


    #Initialise weights
    if model_opts['transpose_sparsity']:
        pmean = 0.;
        pcov = 1.;
        qmean = "random";
        qvar = 1.
        precompute_pcovinv = [(dataX is None) and not(model_opts["ARD_per_view"])] * M
        # Warning : precompute_pcovinv must be True if and only if "AlphaW" and "SigmaAlphaW" not in init.nodes
        init.initW(pmean=pmean, pcov=pcov, qmean=qmean, qvar=qvar, covariates=model_opts['covariatesFiles'], #covariates=dataCovariates,
                   scale_covariates=model_opts['scale_covariates'], precompute_pcovinv=precompute_pcovinv)
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

        init.initSW(pmean_S0=priorW_mean_S0, pmean_S1=priorW_meanS1, pvar_S0=priorW_varS0, pvar_S1=priorW_varS1, ptheta=priorW_theta,
            qmean_S0=initW_meanS0, qmean_S1=initW_meanS1, qvar_S0=initW_varS0, qvar_S1=initW_varS1, qtheta=initW_theta)


    # Initialise ARD or covariance prior structure on W for each view

    # values of the parameters of the prior Gamma distribution
    if 'noise' in model_opts: #simulation
        #pa = 0.01; pb = 0.01
        pa = 1 ; pb = 1.
    else:
        pa = 1e-14 ; pb = 1e-14

    if model_opts['transpose_sparsity']:
        if dataX is not None:
            params = [None] * M
            for m in range(M):
                if view_has_covariance_prior[m]:
                    params[m]={'X':dataX[m],'sigma_clust':dataClust[m],'n_diag':n_diag}
                else:
                    params[m]={'pa': pa, 'pb':pb, 'qa':1., 'qb':1.}
            init.initMixedSigmaAlphaW_mk(view_has_covariance_prior,params)
        else:
            if model_opts["ARD_per_view"]:
                qa = 1.; qb = 1.
                init.initAlphaW_mk(pa=pa, pb=pb, qa=qa, qb=qb)
    else:
        # TODO do sth here for siulations
        qa = 1.; qb = 1.
        init.initAlphaW_mk(pa=pa, pb=pb, qa=qa, qb=qb)


    # Initialise ARD or covariance prior structure on Z

    # values of the parameters of the prior Gamma distribution
    if 'noise' in model_opts:  # simulation
        #pa = 0.01; pb = 0.01
        pa = 1; pb = 1.
    else:
        pa = 1e-14; pb = 1e-14

    if model_opts['transpose_sparsity']:
        qa = 1.; qb = 1.
        init.initAlphaZ_k(pa=pa, pb=pb, qa=qa, qb=qb)
    else:
        if dataX is not None:
            # TODO add a if statement to check if there is a sigma_clust argument to see if blockSigma is needed
            if dataClust is None:
                init.initSigmaZ_k(dataX, n_diag=n_diag)
            else:
                init.initSigmaBlockZ_k(dataX, clust=dataClust, n_diag=n_diag)
        else:
            pass
            #qa = 1.; qb = 1.
            #init.initAlphaZ_k(pa=pa, pb=pb, qa=qa, qb=qb)

    #Initialise sparsity on weights (or factors if transpose)

    # values of the parameters of the prior Beta distribution
    if 'sparsity' in model_opts:
        assert 0. <= model_opts['sparsity'] <= 1., 'sparsty level must be between 0 and 1'
        priorTheta_a = 10.
        priorTheta_b = ((1 - model_opts['sparsity']) / model_opts['sparsity']) * priorTheta_a
    elif 'noise' in model_opts:
        priorTheta_a = 10.
        priorTheta_b = 10.
    else:
        priorTheta_a = 1.
        priorTheta_b = 1.


    if model_opts["transpose_sparsity"]:
        # Initialise sparsity on the factors
        initTheta_a = 1.
        initTheta_b = 0.001 #1.
        #initTheta_E = 1. #mandatory and used if ThetaMixed contains a ThetaConstNode

        # learnTheta = s.ones((D, K))
        # TO-DOOOOOOOOOO
        # if model_opts["learnIntercept"]:
        #     learnTheta[:,0] = 0. # Remove sparsity from the weight vector that will capture the feature-wise means
        learnTheta_ix = np.ones(K)
        if model_opts["learnIntercept"]:
            for ix in learnTheta_ix:
                ix[0] = 0

        init.initThetaMixedZ_k(learnTheta_ix, pa=priorTheta_a, pb=priorTheta_b, qa=initTheta_a, qb=initTheta_b)
        #init.initThetaMixedZ_k(learnTheta_ix, pa=priorTheta_a, pb=priorTheta_b, qa=initTheta_a, qb=initTheta_b,qE=initTheta_E)

    else:
        # Initialise sparsity on the weights
        initTheta_a = 1.
        initTheta_b = 1.
        #initTheta_E = 1. #mandatory and used if ThetaMixed contains a ThetaConstNode

        # learnTheta = [s.ones((D[m], K)) for m in range(M)]
        # TO-DOOOOOOOOOO
        # if model_opts["learnIntercept"]:
        #     for m in range(M):
        #     learnTheta[m][:,0] = 0. # Remove sparsity from the weight vector that will capture the feature-wise means
        learnTheta_ix = [np.ones(K)] * M
        if model_opts["learnIntercept"]:
            for ix in learnTheta_ix:
                ix[0] = 0
        init.initThetaMixedW_mk(learnTheta_ix, pa=priorTheta_a, pb=priorTheta_b, qa=initTheta_a, qb=initTheta_b)
        #init.initThetaMixedW_mk(learnTheta_ix, pa=priorTheta_a, pb=priorTheta_b, qa=initTheta_a, qb=initTheta_b, qE=initTheta_E)


    # Noise


    # values of the parameters of the prior Gamma distribution
    if 'noise' in model_opts:  # simulation
        pa = 20.; pb = pa * model_opts['noise']
    else:
        pa = 1e-14; pb = 1e-14

    qa=1.; qb=1.
    init.initTau(pa=pa, pb=pb, qa=qa, qb=qb, transposed=model_opts["transpose_noise"])


    # Observed data
    init.initY(transpose_noise=model_opts["transpose_noise"])


    ############################################
    ## Define the markov blanket of each node ##
    ############################################

    nodes = init.getNodes()
    if model_opts['transpose_sparsity']:
        nodes["ThetaZ"].addMarkovBlanket(SZ=nodes["SZ"])
        nodes["AlphaZ"].addMarkovBlanket(SZ=nodes["SZ"])
        nodes["SZ"].addMarkovBlanket(AlphaZ=nodes["AlphaZ"], ThetaZ=nodes["ThetaZ"], Y=nodes["Y"], W=nodes["W"],Tau=nodes["Tau"])

        if dataX is not None:
            nodes["SigmaAlphaW"].addMarkovBlanket(W=nodes["W"])
            nodes["W"].addMarkovBlanket(SigmaAlphaW=nodes["SigmaAlphaW"])
        else:
            if model_opts["ARD_per_view"]:
                nodes["AlphaW"].addMarkovBlanket(W=nodes["W"])
                nodes["W"].addMarkovBlanket(AlphaW=nodes["AlphaW"])
        nodes["W"].addMarkovBlanket(Y=nodes["Y"], SZ=nodes["SZ"], Tau=nodes["Tau"])

        nodes["Y"].addMarkovBlanket(SZ=nodes["SZ"], W=nodes["W"], Tau=nodes["Tau"])
        nodes["Tau"].addMarkovBlanket(SZ=nodes["SZ"], W=nodes["W"], Y=nodes["Y"])

    else:
        if dataX is not None:
            nodes["SigmaZ"].addMarkovBlanket(Z=nodes["Z"])
            nodes["Z"].addMarkovBlanket(SigmaZ=nodes["SigmaZ"])
        else:
            pass
            #nodes["AlphaZ"].addMarkovBlanket(Z=nodes["Z"])
            #nodes["Z"].addMarkovBlanket(AlphaZ=nodes["AlphaZ"])
        nodes["Z"].addMarkovBlanket(Y=nodes["Y"], SW=nodes["SW"], Tau=nodes["Tau"])

        nodes["ThetaW"].addMarkovBlanket(SW=nodes["SW"])
        nodes["AlphaW"].addMarkovBlanket(SW=nodes["SW"])
        nodes["SW"].addMarkovBlanket(AlphaW=nodes["AlphaW"], ThetaW=nodes["ThetaW"], Y=nodes["Y"], Z=nodes["Z"],Tau=nodes["Tau"])
        nodes["Y"].addMarkovBlanket(Z=nodes["Z"], SW=nodes["SW"], Tau=nodes["Tau"])
        nodes["Tau"].addMarkovBlanket(Z=nodes["Z"], SW=nodes["SW"], Y=nodes["Y"])

    #################################
    ## Initialise Bayesian Network ##
    #################################

    net = BayesNet(dim=dim, nodes=init.getNodes())

    return net
