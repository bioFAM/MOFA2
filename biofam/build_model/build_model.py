"""
Module with functions to build and initialise the model
"""

import scipy as s
from sys import path
from time import time,sleep
import pandas as pd
import numpy as np
#from joblib import Parallel, delayed

from biofam.core.BayesNet import *
from init_nodes import *
from utils import *

<<<<<<< HEAD
def build_model(data, model_opts):
    """Method to run a single trial of a MOFA model
    data: 
=======
def runSingleTrial(data, data_opts, model_opts, train_opts, seed=None, trial=1, verbose=False, NetType='Sampler'):
    """Method to run a single trial of a MOFA model
    data:
    data_opts
>>>>>>> d0c0d1e1275991f7bcc5f30147dc58b6c6dbb4e3
    model_opts:
    train_opts:

    PARAMETERS
    ----------
    """

    print ("\n")
    print ("#"*24)
    print ("## Building the model ##")
    print ("#"*24)
    print ("\n")
    sleep(1)

    # Define dimensionalities
    M = len(data)
    N = data[0].shape[0]
    D = s.asarray([ data[m].shape[1] for m in range(M) ])
    K = model_opts["K"]

    dim = {'M':M, 'N':N, 'D':D, 'K':K }

    ###########################
    ## Do some sanity checks ##
    ###########################

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
    init = initModel(dim, data, model_opts["likelihood"], seed=seed)

    # Latent variables
    pmean = 0.; pvar = 1.; qmean = "random"; qvar = 1.
    init.initZ(pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, covariates=model_opts['covariates'], scale_covariates=model_opts['scale_covariates'])

    # Weights
    priorSW_mean_S0=0.; priorSW_meanS1=0.; priorSW_varS0=1.; priorSW_varS1=1.; priorSW_theta=1.; 
    initSW_meanS0=0.; initSW_meanS1=0.; initSW_varS0=1.; initSW_varS1=1;, initSW_theta=1.;
    initSW_qEW_S0=0.; initSW_qEW_S1=0.; initSW_qES=1.;
    
    # TO-DOOOOOOOOOO: fIX LEARN INTERCEPT
    # if model_opts["learnIntercept"]:
    # for m in range(M):
    #     if model_opts["likelihood"][m]=="gaussian":
    #         initSW_meanS1[m][:,0] = data[m].mean(axis=0)
    #         initSW_varS1[m][:,0] = 1e-5
    #         initSW_theta[m][:,0] = 1.

    init.initSW(self, pmean_S0=priorSW_mean_S0, pmean_S1=priorSW_meanS1, pvar_S0=priorSW_varS0, pvar_S1=priorSW_varS1, ptheta=priorSW_theta, 
        qmean_S0=initSW_meanS0, qmean_S1=initSW_meanS1, qvar_S0=initSW_varS0, qvar_S1=initSW_varS1, qtheta=initSW_theta, 
        qEW_S0=initSW_qEW_S0, qEW_S1=initSW_qEW_S1, qES=initSW_qES)
    

    # ARD on weights
    pa=1e-14; pb=1e-14; qa=1.; qb=1.; qE=1.
    init.initAlphaW_mk(pa=pa, pb=pb, qa=qa, qb=qb)

    # Precision of noise
    pa=1e-14; pb=1e-14; qa=1.; qb=1.; qE=1.
    init.initTau(pa=pa, pb=pb, qa=qa, qb=qb)

    # Sparsity on the weights
    learnTheta = [ s.ones((self.D[m],self.K)) for m in xrange(self.M) ]
    priorTheta_a = 1.
    priorTheta_b = 1.
    initTheta_a = 1.
    initTheta_b = 1.
    initTheta_E = 1.
    # TO-DOOOOOOOOOO
    # if model_opts["learnIntercept"]:
    #     for m in range(M): 
    #     learnTheta[m][:,0] = 0. # Remove sparsity from the weight vector that will capture the feature-wise means 
    init.initThetaMixed(pa=priorTheta_a, pb=priorTheta_b, qa=initTheta_a,  qb=initTheta_b, qE=initTheta_E, learnTheta=learnTheta)

    # Observed data
    init.initY()

    ############################################
    ## Define the markov blanket of each node ##
    ############################################

    nodes = init.getNodes()
    nodes["Z"].addMarkovBlanket(SW=nodes["SW"], Tau=nodes["Tau"], Y=nodes["Y"])
    nodes["Theta"].addMarkovBlanket(SW=nodes["SW"])
    nodes["AlphaW"].addMarkovBlanket(SW=nodes["SW"])
    nodes["SW"].addMarkovBlanket(Z=nodes["Z"], Tau=nodes["Tau"], Alpha=nodes["AlphaW"], Y=nodes["Y"], Theta=nodes["Theta"])
    nodes["Y"].addMarkovBlanket(Z=nodes["Z"], SW=nodes["SW"], Tau=nodes["Tau"])
    nodes["Tau"].addMarkovBlanket(Z=nodes["Z"], SW=nodes["SW"], Y=nodes["Y"])

    #################################
    ## Initialise Bayesian Network ##
    #################################

    net = BayesNet(dim=dim, nodes=init.getNodes())

    return net


def train_model(bayesnet, train_opts):

    # QC on the Bayesian Network
    assert type(bayesnet)==BayesNet, "'bayesnet' has to be a BayesNet class"

    # Define training options
    bayesnet.setTrainOptions(train_opts)

    ####################
    ## Start training ##
    ####################

    print ("\n")
    print ("#"*45)
    print ("## Training the model with seed %d ##" % (train_opts['seed']))
    print ("#"*45)
    print ("\n")
    sleep(1)
<<<<<<< HEAD

    bayesnet.iterate()
=======

    if NetType == 'VB':
        net.iterate()
    if NetType == 'Sampler':
        net.simulate()

    return net

def runMultipleTrials(data, data_opts, model_opts, train_opts, keep_best_run, seed=None, verbose=True):

    """Method to run multiple trials of a MOFA model

    PARAMETERS
    -----
    data:
    data_opts
    model_opts:
    train_opts:
    seed:
    trial:
    verbose:
    """
    # trained_models = Parallel(n_jobs=train_opts['cores'], backend="threading")(
    # trained_models = Parallel(n_jobs=train_opts['cores'])(
    #     delayed(runSingleTrial)(data,data_opts,model_opts,train_opts,seed,i) for i in range(1,train_opts['trials']+1))
    trained_models = [ runSingleTrial(data,data_opts,model_opts,train_opts,seed,i) for i in range(1,train_opts['trials']+1) ]
>>>>>>> d0c0d1e1275991f7bcc5f30147dc58b6c6dbb4e3

    print("\n")
    print("#"*43)
    print("## Training finished ##")
    print("#"*43)
    print("\n")

<<<<<<< HEAD
    return (bayesnet)
    
    

=======
    #####################
    ## Process results ##
    #####################

    # Select the trial with the best lower bound or keep all models
    if train_opts['trials'] > 1:
        if keep_best_run:
            lb = map(lambda x: x.getTrainingStats()["elbo"][-1], trained_models)
            save_models = [ trials[s.argmax(lb)] ]
            outfiles = [ data_opts['outfile'] ]
        else:
            save_models = trained_models
            tmp = os.path.splitext(data_opts['outfile'])
            outfiles = [ tmp[0]+"_"+str(t)+tmp[1]for t in range(train_opts['trials']) ]
    else:
        save_models = trained_models
        outfiles = [ data_opts['outfile'] ]

    ##################
    ## Save results ##
    ##################

    sample_names = data[0].index.tolist()
    feature_names = [  data[m].columns.values.tolist() for m in range(len(data)) ]
    for t in range(len(save_models)):
        print("Saving model %d in %s...\n" % (t,outfiles[t]))
        saveModel(save_models[t], outfile=outfiles[t], view_names=data_opts['view_names'],
            sample_names=sample_names, feature_names=feature_names, train_opts=train_opts, model_opts=model_opts)
>>>>>>> d0c0d1e1275991f7bcc5f30147dc58b6c6dbb4e3
