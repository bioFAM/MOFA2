"""
Module to train a bioFAM model
"""

import scipy as s
import pandas as pd
import numpy as np

from biofam.core.BayesNet import BayesNet


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

    bayesnet.iterate()


    print("\n")
    print("#"*43)
    print("## Training finished ##")
    print("#"*43)
    print("\n")

    return (bayesnet)
