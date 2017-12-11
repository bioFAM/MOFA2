import scipy as s
from .basic_distributions import Distribution
from .bernoulli import Bernoulli
from .univariate_gaussian import UnivariateGaussian

from biofam.core.utils import *


class BernoulliGaussian(Distribution):
    """
    Class to define a Bernoulli-Gaussian distributions (for more information see Titsias and Gredilla, 2014)

    The best way to characterise a joint Bernoulli-Gaussian distribution P(w,s) is by considering
    its factorisation p(w|s)p(s) where s is a bernoulli distribution and w|s=0 and w|s=1 are normal distributions

    Equations:
    p(w,s) = Normal(w|mean,var) * Bernoulli(s|theta)
    FINISH EQUATIONS

    ROOM FOR IMPROVEMENT: i think the current code is inefficient because you have to keep track of the params
    and expectations in both the factorised distributions and the joint one. I think
    perhaps is better not to define new Bernoulli and UnivariateGaussian but work directly with the joint model

    """
    def __init__(self, dim, mean_S0, mean_S1, var_S0, var_S1, theta, EW_S0=None, EW_S1=None, ES=None):
        Distribution.__init__(self,dim)
        self.S = Bernoulli(dim=dim, theta=theta, E=ES)
        self.W_S0 = UnivariateGaussian(dim=dim, mean=mean_S0, var=var_S0, E=EW_S0)
        self.W_S1 = UnivariateGaussian(dim=dim, mean=mean_S1, var=var_S1, E=EW_S1)

        # Collect parameters
        self.params = { 'mean_S0':self.W_S0.params['mean'],
                        'mean_S1':self.W_S1.params['mean'],
                        'var_S0':self.W_S0.params['var'],
                        'var_S1':self.W_S1.params['var'],
                        'theta':self.S.params['theta'] }

        # Collect expectations
        self.updateExpectations()

    def setParameters(self,**params):
        # Setter function for parameters
        self.S.setParameters(theta=params['theta'])
        self.W_S0.setParameters(mean=params['mean_S0'], var=params['var_S0'])
        self.W_S1.setParameters(mean=params['mean_S1'], var=params['var_S1'])
        self.params = params

    def updateParameters(self):
        # Method to update the parameters of the joint distribution based on its constituent distributions
        self.params = { 'theta':self.S.params["theta"],
                        'mean_S0':self.W_S0.params["mean"], 'var_S0':self.W_S0.params["var"],
                        'mean_S1':self.W_S1.params["mean"], 'var_S1':self.W_S1.params["var"] }

    def updateExpectations(self):
        # Method to calculate the expectations based on the current estimates for the parameters

        # Update expectations of the constituent distributions
        self.S.updateExpectations()
        self.W_S0.updateExpectations()
        self.W_S1.updateExpectations()

        # Calculate expectations of the joint distribution
        ES = self.S.getExpectation()
        EW = self.W_S1.getExpectation()
        E = ES * EW
        ESWW = ES * (s.square(EW) + self.params["var_S1"])
        # ESWW = self.params["theta"] * (self.params["mean_S1"]**2 + self.params["var_S1"])
        EWW = ES*(s.square(EW)+self.params["var_S1"]) + (1-ES)*self.params["var_S0"]
        # EWW = self.params["theta"]*(self.params["mean_S1"]**2+self.params["var_S1"]) + (1-self.params["theta"])*self.params["var_S0"]

        # Collect expectations
        self.expectations = {'E':E, 'ES':ES, 'EW':EW, 'ESWW':ESWW, 'EWW':EWW }

    def removeDimensions(self, axis, idx):
        # Method to remove undesired dimensions
        # - axis (int): axis from where to remove the elements
        # - idx (numpy array): indices of the elements to remove
        assert axis <= len(self.dim)
        assert s.all(idx < self.dim[axis])
        self.S.removeDimensions(axis,idx)
        self.W_S0.removeDimensions(axis,idx)
        self.W_S1.removeDimensions(axis,idx)
        self.updateParameters()
        self.updateExpectations()

    def updateDim(self, axis, new_dim):
        # Function to update the dimensionality of a particular axis
        self.S.updateDim(axis,new_dim)
        self.W_S0.updateDim(axis,new_dim)
        self.W_S1.updateDim(axis,new_dim)
        dim = list(self.dim)
        dim[axis] = new_dim
        self.dim = tuple(dim)
