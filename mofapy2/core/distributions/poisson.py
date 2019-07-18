import scipy as s
import scipy.stats as stats
from .basic_distributions import Distribution

from mofapy2.core.utils import *

class Poisson(Distribution):
    """
    Class to define Poisson distributions.

    Equations:
    p(x|theta) = theta**x * exp(-theta) * 1/theta!
    log p(x|a,b) = x*theta - theta - log(x!)
    E[x] = theta
    var[x] = theta
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
        assert x.dtype == int, "x has to be an integer array"
        theta = self.params['theta'].flatten()
        x = x.flatten()
        # return s.prod (stats.poisson.pmf(x,theta) )
        return s.prod( s.divide(theta**x * s.exp(-theta),s.misc.factorial(x)) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        assert x.dtype == int, "x has to be an integer array"
        theta = self.params['theta'].flatten()
        x = x.flatten()
        # return s.log( s.prod (stats.poisson.pmf(x,theta) ))
        return s.sum( x*s.log(theta) - theta - s.log(s.misc.factorial(x)) )
