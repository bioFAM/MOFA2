import scipy as s
import scipy.special as special
import scipy.stats as stats
from .basic_distributions import Distribution

from mofapy2.core.utils import *

class Binomial(Distribution):
    """
    Class to define Binomial distributions

    Equations:
    p(x|N,theta) = binom(N,x) * theta**(x) * (1-theta)**(N-x)
    log p(x|N,theta) = log(binom(N,x)) + x*theta + (N-x)*(1-theta)
    E[x] = N*theta
    var[x] = N*theta*(1-theta)
    """
    def __init__(self, dim, N, theta, E=None):
        Distribution.__init__(self, dim)

        # Initialise parameters
        theta = s.ones(dim)*theta
        N = s.ones(dim)*N
        self.params = { 'theta':theta, 'N':N }

        # Initialise expectations
        if E is None:
            self.updateExpectations()
        else:
            E = s.ones(dim)*E
            self.expectations = { 'E':E }

        # Check that dimensionalities match
        self.CheckDimensionalities()

    def updateExpectations(self):
        E = self.params["N"] * self.params["N"]
        self.expectations = { 'E':E }

    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        assert x.dtype == int, "x has to be an integer array"
        # return s.prod( stats.binom.pmf(x, self.params["N"], self.theta) )
        return s.prod( special.binom(self.params["N"],x) * self.params["theta"]**x * (1-self.params["theta"])**(self.params["N"]-x) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        assert x.dtype == int, "x has to be an integer array"
        # print s.sum (stats.binom.logpmf(x, self.params["N"], self.theta) )
        return s.sum( s.log(special.binom(self.params["N"],x)) + x*s.log(self.params["theta"]) + (self.params["N"]-x)*s.log(1-self.params["theta"]) )

    def sample(self):
        return s.random.binomial(self.params['N'], self.params['theta'])
