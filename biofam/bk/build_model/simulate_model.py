"""
Module to simulate a bioFAM model
"""

import scipy as s
import pandas as pd
import numpy as np

from biofam.core.BayesNet import BayesNet


def simulate_model(bayesnet):

    # QC on the Bayesian Network
    assert type(bayesnet)==BayesNet, "'bayesnet' has to be a BayesNet class"

    ####################
    ## Start training ##
    ####################

    print ("\n")
    print ("#"*45)
    print ("## Simulating the model")
    print ("#"*45)
    print ("\n")

    bayesnet.simulate()


    print("\n")
    print("#"*43)
    print("## Simulation finished ##")
    print("#"*43)
    print("\n")
