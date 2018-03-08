
"""
Module to initialise the nodes in the new model
"""

import scipy as s
import scipy.stats as stats
from sys import path
import sklearn.decomposition


from biofam.core.nodes import *

# specific for spatialFA
from biofam.build_model.init_model import initModel



class initNewModel(initModel):
    def __init__(self, dim, data, lik):
        """
        PARAMETERS
        ----------
         dim: dictionary
            keyworded dimensionalities: N for the number of samples, M for the number of views, K for the number of latent variables, D for the number of features per view (a list)
         data: list of ndarrays of length M:
            observed data
         lik: list of strings
            likelihood for each view
        """
        super(initNewModel, self).__init__(dim, data, lik)
