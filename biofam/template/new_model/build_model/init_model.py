
"""
Module to initalise a bioFAM model
"""

import scipy as s
import scipy.stats as stats
from sys import path
import sklearn.decomposition

from biofam.core.nodes import *
from biofam.build_model.init_model import initModel


class initModelNewModel(initModel):
    def __init__(self, dim, data, lik):
        """
        PARAMETERS
        ----------
         dim: dictionary
            keyworded dimensionalities: N for the number of samples, M for the number of views, K for the number of latent variables, D for the number of features per view (a list)
         data: list of length M with ndarrays of dimensionality (N,Dm):
            observed data
         lik: list of length M with strings
            likelihood for each view
        """
        super(initModelNewModel, self).__init__(dim, data, lik)

    # implement intialisation of new nodes here 
