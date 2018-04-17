import scipy as s
import numpy.linalg as linalg
import scipy.stats as stats
from .basic_distributions import Distribution

from biofam.core.utils import *


class MultivariateGaussian(Distribution):
    """
    Class to define multivariate Gaussian distribution.
    This class can store N multivate Gaussian Distributions of dimensionality D each.

    Equations:
    p(X|Mu,Sigma) = 1/(2pi)^{D/2} * 1/(|Sigma|^0.5) * exp( -0.5*(X-Mu)^{T} Sigma^{-1} (X-Mu) )
    log p(X|Mu,Sigma) = -D/2*log(2pi) - 0.5log(|Sigma|) - 0.5*(X-Mu)^{T} Sigma^{-1} (X-Mu)
    E[X] = Mu
    E[X^2] = E[X]E[X] + Cov[X] = MuMu + Sigma
    cov[X] = Sigma
    H[X] = 0.5*log(|Sigma|) + D/2*(1+log(2pi))

    Dimensionalities:
    - X: (N,D)
    - Mu: (N,D)
    - Sigma: (N,D,D)
    - E[X]: (N,D)
    - E[X^2]: (N,D,D)
    """
    def __init__(self, dim, mean, cov, E=None):
        Distribution.__init__(self, dim)

        # Check dimensions are correct
        assert len(dim) == 2, "You need to specify two dimensions for 'dim': (number of distributions, dimensionality) "
        assert not (dim[0]==1 and dim[1]==1), "A 1x1 multivariate normal reduces to a Univariate normal "

        ## Initialise the parameters ##

        # Initialise the mean
        # If 'mean' is a scalar, broadcast it to all dimensions
        if isinstance(mean,(int,float)): mean = s.ones( (dim[0],dim[1]) ) * mean
        # If 'mean' has dim (D,) and we have N distributions, broadcast it to all N distributions
        if len(mean.shape)==1 and mean.shape[0]==dim[1]: mean = s.repeat(mean,dim[0],0)
        assert sum(mean.shape) > 2, "The mean has to be a matrix with shape (N,D) "

        # Initialise the covariance
        # If 'cov' is a matrix and not a tensor, broadcast it along the zeroth axis
        if len(cov.shape) == 2: cov = s.repeat(cov[None,:,:],dim[0],0)
        assert (cov.shape[1]==cov.shape[2]) and (sum(cov.shape[1:])>1), "The covariance has to be a tensor with shape (N,D,D)"

        # Check that the dimensionalities of 'mean' and 'cov' match
        # TODO sort out
        # assert cov.shape[1] == mean.shape[1] == dim[1], "Error in the dimensionalities"
        # assert cov.shape[0] == mean.shape[0] == dim[0], "Error in the dimensionalities"

        self.params = {'mean':mean, 'cov':cov }

        # Initialise expectations
        if E is None:
            self.updateExpectations()
        else:
            self.expectations = { 'E':E }


    def updateExpectations(self):
        # Update first and second moments using current parameters
        E = self.params['mean']

        # self.E2 = s.empty( (self.dim[0],self.dim[1],self.dim[1]) )
        # for i in range(self.dim[0]):
        #     self.E2[i,:,:] = s.outer(self.E[i,:],self.E[i,:]) + self.cov[i,:,:]

        E2 = self.params['cov'].copy()
        # TODO sort out index
        # import pdb; pdb.set_trace()
        for i in range(self.dim[1]):
            E2[i,:,:] += s.outer(E[:,i],E[:,i])

        self.expectations = {'E':E, 'E2':E2}

    def density(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        return s.sum( stats.multivariate_normal.pdf(x, mean=self.params['mean'][n,:], cov=self.params['cov'][n,:,:]) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        l = 0.
        D = self.dim[1]
        for n in range(self.dim[0]):
            qterm = (x[n,:]-self.params['mean'][n,:]).T.dot(linalg.det(self.params['cov'][n,:,:])).dot(x[n,:]-self.params['mean'][n,:])
            l += -0.5*D*s.log(2*s.pi) - 0.5*s.log(linalg.det(self.params['cov'][n,:,:])) -0.5*qterm
        return l
        # return s.sum( s.log(stats.multivariate_normal.pdf(x, mean=self.mean[n,:], cov=self.cov[n,:,:])) )

    def removeDimensions(self, axis, idx):
        # Method to remove undesired dimensions
        # - axis (int): axis from where to remove the elements
        # - idx (numpy array): indices of the elements to remove
        assert axis <= len(self.dim)
        assert s.all(idx < self.dim[axis])
        self.params["mean"] = s.delete(self.params["mean"], axis=1, obj=idx)
        # self.params["cov"] = s.delete(self.params["cov"], axis=1, obj=idx)
        # self.params["cov"] = s.delete(self.params["cov"], axis=2, obj=idx)
        self.params["cov"] = s.delete(self.params["cov"], axis=0, obj=idx)
        self.expectations["E"] = s.delete(self.expectations["E"], axis=1, obj=idx)
        # self.expectations["E2"] = s.delete(self.expectations["E2"], axis=1, obj=idx)
        # self.expectations["E2"] = s.delete(self.expectations["E2"], axis=2, obj=idx)
        self.expectations["E2"] = s.delete(self.expectations["E2"], axis=0, obj=idx)
        self.dim = (self.dim[0],self.dim[1]-len(idx))

    def sample(self):
        return s.random.multivariate_normal(self.params['mean'], self.params['cov'])

    # def entropy(self):
        # CHECK THIs Is CORRECT
        # tmp = sum( [ logdet(self.cov[i,:,:]) for i in range(self.dim[0]) ] )
        # return ( 0.5*(tmp + (self.dim[0]*self.dim[1])*(1+s.log(2*pi)) ).sum() )