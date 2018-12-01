"""
Module to train a bioFAM model
"""

import scipy as s
import pandas as pd
import numpy as np

from biofam.core.BayesNet import BayesNet


def train_model(model):

    # Sanity check on the Bayesian Network
    assert isinstance(model, BayesNet), "'bayesnet' has to be a BayesNet class"


    ####################
    ## Start training ##
    ####################

    print ("\n")
    print ("#"*40)
    print ("## Training the model with seed %d ##" % (model.options['seed']))
    print ("#"*40)
    print ("\n")

    model.iterate()

    print("\n")
    print("#"*23)
    print("## Training finished ##")
    print("#"*23)
    print("\n")
