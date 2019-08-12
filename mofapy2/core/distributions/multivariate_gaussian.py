import scipy as s
import numpy.linalg as linalg
import scipy.stats as stats
from .basic_distributions import Distribution

from mofapy2.core.utils import *

# TODO : remove a loop

class MultivariateGaussian(Distribution):
    """
    Class to define multivariate Gaussian distribution.
    This class can store:
    - if axis_cov=1 : N multivariate Gaussian Distributions of dimensionality D each
       (each line of the X matrix is a multivariate Gaussian)
    - if axis_cov=0 : D multivariate Gaussian Distributions of dimensionality N each
       (each column of the X matrix is a multivariate Gaussian)

    Equations (for axis_cov=1) :
    p(X|Mu,Sigma) = 1/(2pi)^{D/2} * 1/(|Sigma|^0.5) * exp( -0.5*(X-Mu)^{T} Sigma^{-1} (X-Mu) )
    log p(X|Mu,Sigma) = -D/2*log(2pi) - 0.5log(|Sigma|) - 0.5*(X-Mu)^{T} Sigma^{-1} (X-Mu)
    E[X] = Mu
    E[X^2] = E[X]E[X] + Cov[X] = MuMu + Sigma
    cov[X] = Sigma
    H[X] = 0.5*log(|Sigma|) + D/2*(1+log(2pi))

    Dimensionalities :
    - X: (N,D)
    - Mu: (N,D)
    - Sigma: (N,D,D) if axis_cov=1, (D,N,N) if axis_cov=0, as a list of N matrices (D,D)
    - E[X]: (N,D)
    - E[X^2]: (N,D,D) if axis_cov=1, (D,N,N) if axis_cov=0, , as a list of N matrices (D,D)
    """
    def __init__(self, dim, mean, cov, axis_cov=1, E=None):
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
        if isinstance(cov,s.ndarray):
            if len(cov.shape) == 2:
                assert ((axis_cov == 0)or(axis_cov == 1)), "Error : axis_cov is the index of the dimension of the covariance matrix, either 0 or 1"
                if axis_cov == 1 :
                    cov = [cov] * dim[0]
                else:
                    cov = [cov] * dim[1]

        assert (cov[0].shape[0] == cov[0].shape[1]) and (sum(cov[0].shape) > 1), "The covariance has to be a tensor with shape (N,D,D) or (D,N,N)"

        # Check that the dimensionalities of 'mean' and 'cov' match
        if axis_cov == 1:
            assert cov[0].shape[0] == mean.shape[1] == dim[1], "The covariance has to be a tensor with shape (N,D,D)"
            assert len(cov) == mean.shape[0] == dim[0], "The covariance has to be a tensor with shape (N,D,D)"
        else:
            assert cov[0].shape[0] == mean.shape[0] == dim[0], "The covariance has to be a tensor with shape (D,N,N)"
            assert len(cov) == mean.shape[1] == dim[1], "The covariance has to be a tensor with shape (D,N,N)"

        self.axis_cov = axis_cov
        self.params = {'mean':mean, 'cov':cov }

        # Initialise expectations
        if E is None:
            self.updateExpectations()
        else:
            self.expectations = { 'E':E }

    def updateExpectations(self):
        # Update first and second moments, and expectation of X*X^T, using current parameters
        E = self.params['mean']

        # self.E2 = s.empty( (self.dim[0],self.dim[1],self.dim[1]) )
        # for i in range(self.dim[0]):
        #     self.E2[i,:,:] = s.outer(self.E[i,:],self.E[i,:]) + self.cov[i,:,:]

        # computing the expectation of X*X.T
        # Work but not useful now !
        #EXXT = self.params['cov'].copy()
        # TODO sort out index

        #if self.axis_cov == 1:
        #    for i in range(self.dim[0]):
        #        EXXT[i,:,:] += s.outer(E[i,:],E[i,:])
        #else:
        #    for i in range(self.dim[1]):
        #        EXXT[i,:,:] += s.outer(E[:,i],E[:,i])

        # from the expectation of X*X.T to the expectation of X^2
        # Work but not useful now !
        #E2 = np.zeros((self.dim[0], self.dim[1]))
        #if self.axis_cov == 1:
        #    for i in range(self.dim[0]):
        #        E2[i, :] = np.diag(EXXT[i, :, :]).flatten()  # extracting the diagonal
        #else:
        #    for i in range(self.dim[1]):
        #        E2[:, i] = np.diag(EXXT[i, :, :]).flatten()  # extracting the diagonal

        self.expectations = {'E': E}
        #self.expectations = {'E':E, 'E2':E2, 'EXXT':EXXT}

    #def density(self, x):
    #    assert x.shape == self.dim, "Problem with the dimensionalities"
    #    return s.sum( stats.multivariate_normal.pdf(x, mean=self.params['mean'][n,:], cov=self.params['cov'][n,:,:]) )

    def loglik(self, x):
        assert x.shape == self.dim, "Problem with the dimensionalities"
        l = 0.

        if self.axis_cov == 1:
            D = self.dim[1]
            for n in range(self.dim[0]):
                qterm = (x[n,:]-self.params['mean'][n,:]).T.dot(linalg.det(self.params['cov'][n])).dot(x[n,:]-self.params['mean'][n,:])
                l += -0.5*D*s.log(2*s.pi) - 0.5*s.log(linalg.det(self.params['cov'][n])) -0.5*qterm
            # return s.sum( s.log(stats.multivariate_normal.pdf(x, mean=self.mean[n,:], cov=self.cov[n,:,:])) )

        else:
            N = self.dim[0]
            for d in range(self.dim[1]):
                qterm = (x[:, d] - self.params['mean'][:, d]).dot(linalg.det(self.params['cov'][d])).dot((x[:, d] - self.params['mean'][:, d]).T)
                l += -0.5 * N * s.log(2 * s.pi) - 0.5 * s.log(linalg.det(self.params['cov'][d])) - 0.5 * qterm

        return l

    def removeDimensions(self, axis, idx):
        # Method to remove undesired dimensions
        # - axis (int): axis from where to remove the elements
        # - idx (numpy array): indices of the elements to remove
        assert axis <= len(self.dim)
        assert s.all(idx < self.dim[axis])

        self.params["mean"] = s.delete(self.params["mean"], axis=axis, obj=idx)
        self.expectations["E"] = s.delete(self.expectations["E"], axis=axis, obj=idx)
        #self.expectations["E2"] = s.delete(self.expectations["E2"], axis=axis, obj=idx)

        if self.axis_cov == 1: #cov has shape (a,b,b) when mean has shape (a,b)
            if axis == 0:
                for i in idx: del self.params["cov"][i]
            else:
                self.params["cov"] = s.delete(self.params["cov"], axis=0, obj=idx)
                self.params["cov"] = s.delete(self.params["cov"], axis=1, obj=idx)

        else: #cov has shape (b,a,a) when mean has shape (a,b)
            if axis == 0:
                self.params["cov"] = s.delete(self.params["cov"], axis=0, obj=idx)
                self.params["cov"] = s.delete(self.params["cov"], axis=1, obj=idx)
            else:
                for i in idx: del self.params["cov"][i]

        self.updateDim(axis=axis, new_dim=self.dim[axis] - len(idx))

    def sample(self):
        if axis_cov==1:
            samples = []
            for n in range(self.dim[0]):
                samples.append(s.random.multivariate_normal(self.params['mean'][n,:], self.params['cov'][n]))
            samples = np.array(samples)
        else:
            samples = []
            for d in range(self.dim[1]):
                samples.append(s.random.multivariate_normal(self.params['mean'][:,d], self.params['cov'][d]))
            samples = np.array(samples).T
        return samples

    # def entropy(self):
        # CHECK THIs Is CORRECT
        # tmp = sum( [ logdet(self.cov[i,:,:]) for i in range(self.dim[0]) ] )
        # return ( 0.5*(tmp + (self.dim[0]*self.dim[1])*(1+s.log(2*pi)) ).sum() )
