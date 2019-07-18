import scipy as s
import scipy.stats as stats
from .basic_distributions import Distribution

from mofapy2.core.utils import *

class UnivariateGaussian(Distribution):
    """
    Class to define univariate Gaussian distributions

    Equations:
    Class for a univariate Gaussian distributed node
    p(x|mu,sigma^2) = 1/sqrt(2*pi*sigma^2) * exp(-0.5*(x-mu)^2/(sigma^2) )
    log p(x|mu,sigma^2) =
    E[x] = mu
    var[x] = sigma^2
    H[x] = 0.5*log(sigma^2) + 0.5*(1+log(2pi))

    """
    def __init__(self, dim, mean, var, E=None, E2=None):
        Distribution.__init__(self, dim)

        # Initialise parameters
        mean = s.ones(dim) * mean
        var = s.ones(dim) * var
        self.params = { 'mean':mean, 'var':var }

        # Initialise expectations
        self.expectations = {}
        if E is None:
            self.updateExpectations()
        else:
            self.expectations['E'] = s.ones(dim)*E

        if E2 is not None:
            self.expectations['E2'] = s.ones(dim)*E2

        # Check that dimensionalities match
        self.CheckDimensionalities()

    def updateExpectations(self):
        # Update first and second moments using current parameters
        E = self.params['mean']
        E2 = E**2 + self.params['var']
        self.expectations = { 'E':E, 'E2':E2 }

    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        # print stats.norm.pdf(x, loc=self.mean, scale=s.sqrt(self.var))
        return s.sum( (1/s.sqrt(2*s.pi*self.params['var'])) * s.exp(-0.5*(x-self.params['mean'])**2/self.params['var']) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        # return s.log(stats.norm.pdf(x, loc=self.mean, scale=s.sqrt(self.var)))
        return s.sum( -0.5*s.log(2*s.pi) - 0.5*s.log(self.params['var']) -0.5*(x-self.params['mean'])**2/self.params['var'] )

    def entropy(self):
        return s.sum( 0.5*s.log(self.params['var']) + 0.5*(1+s.log(2*s.pi)) )

    def sample(self):
        return s.random.normal(self.params['mean'], np.sqrt(self.params['var']))
