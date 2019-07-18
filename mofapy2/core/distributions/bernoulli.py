import scipy as s
from .basic_distributions import Distribution

from mofapy2.core.utils import *

class Bernoulli(Distribution):
    """
    Class to define Bernoulli distributions

    Equations:

    """
    def __init__(self, dim, theta, E=None):
        Distribution.__init__(self, dim)

        # Initialise parameters
        theta = s.ones(dim) * theta
        self.params = { 'theta':theta }

        # Initialise expectations
        if E is None:
            self.updateExpectations()
        else:
            self.expectations = { 'E':s.ones(dim)*E }

        # Check that dimensionalities match
        self.CheckDimensionalities()

    def updateExpectations(self):
        E = self.params['theta']
        self.expectations = { 'E':E }

    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        return s.prod( self.params['theta']**x * (1-self.params['theta'])**(1-x) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        return s.sum( x*s.log(self.params['theta']) + (1-x)*s.log(1-self.params['theta']) )

    def sample(self):
        return s.random.binomial(1, self.params['theta'])
