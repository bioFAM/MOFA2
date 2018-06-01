
"""
This module is used to define the class containing the entire Bayesian Network,
and the corresponding attributes/methods to train the model, set algorithmic options, calculate lower bound, etc.
"""

from __future__ import division
from time import time
import os
import scipy as s
import pandas as pd
import sys

from biofam.core.nodes.variational_nodes import Variational_Node
from .utils import nans

class BayesNet(object):
    def __init__(self, dim, nodes):
        """ Initialisation of a Bayesian network

        PARAMETERS
        ----------
        dim: dict
            keyworded dimensionalities, ex. {'N'=10, 'M'=3, ...}
        nodes: dict
            dictionary with all nodes where the keys are the name of the node and the values are instances the 'Node' class
        """

        self.dim = dim
        self.nodes = nodes
        self.options = None  # TODO rename to train_options everywhere

        # Training and simulations flag
        self.trained = False
        self.simulated = False

    def setTrainOptions(self, train_opts):
        """ Method to store training options """

        # Sanity checks
        assert "maxiter" in train_opts, "'maxiter' not found in the training options dictionary"
        assert "start_drop" in train_opts, "'start_drop' not found in the training options dictionary"
        assert "freq_drop" in train_opts, "'freq_drop' not found in the training options dictionary"
        assert "verbose" in train_opts, "'verbose' not found in the training options dictionary"
        assert "tolerance" in train_opts, "'tolerance' not found in the training options dictionary"
        assert "forceiter" in train_opts, "'forceiter' not found in the training options dictionary"
        assert "schedule" in train_opts, "'schedule' not found in the training options dictionary"
        assert "start_sparsity" in train_opts, "'start_sparsity' not found in the training options dictionary"

        self.options = train_opts

        # TODO for testing purpose
        self.options['step_tau'] = 1.
        self.options['forgetting_rate'] = 1.
        self.options['batch_size'] = .2
        # self.options['batch_size'] = None

    def getParameters(self, *nodes):
        """ Method to collect all parameters of a given set of nodes

        PARAMETERS
        ----------
        nodes: iterable
            name of the nodes (all nodes by default)
        """

        if len(nodes) == 0: nodes = self.nodes.keys()
        params = {}
        for node in nodes:
            tmp = self.nodes[node].getParameters()
            if tmp != None: params[node] = tmp
        return params

    def getExpectations(self, only_first_moments=False, *nodes):
        """Method to collect all expectations of a given set of nodes

        PARAMETERS
        ----------
        only_first_moments: bool
            get only first moments? (Default is False)
        nodes: list
            name of the nodes (Default is all nodes)
        """

        if len(nodes) == 0: nodes = self.nodes.keys()
        expectations = {}
        for node in nodes:
            if only_first_moments:
                tmp = self.nodes[node].getExpectation()
            else:
                tmp = self.nodes[node].getExpectations()
            expectations[node] = tmp
        return expectations

    def getNodes(self):
        """ Method to return all nodes """
        return self.nodes

    def simulate(self, dist='P'):
        """ Method to simulate from the generative model """
        if 'SW' in self.nodes:
            self.nodes["SW"].sample(dist)
            self.nodes["Z"].sample(dist)
            self.nodes["Z"].samp -= self.nodes["Z"].samp.mean() #centering
        else:
            self.nodes["SZ"].sample(dist)
            self.nodes["W"].sample(dist)
            for m in range(self.dim["M"]):
                self.nodes["W"].nodes[m].samp -= self.nodes["W"].nodes[m].samp.mean() #centering
        self.nodes['Tau'].sample(dist)
        self.nodes['Y'].sample(dist)

        self.simulated = True

    def sampleData(self):
        """ Method to sample data from the prior distributions of the generative model """
        if ~self.simulated:
            self.simulate()
        return self.nodes['Y'].sample(dist='P')

    def saveData(self):
        # TODO some function here to save simulated data
        pass

    def removeInactiveFactors(self, by_r2=None):
        """Method to remove inactive factors

        PARAMETERS
        ----------
        by_r2: float
            threshold to shut down factors based on the coefficient of determination
        """
        drop_dic = {}

        if by_r2 is not None:

            Z = self.nodes['Z'].getExpectation()
            W = self.nodes["W"].getExpectation()

            Y = self.nodes["Y"].getExpectation()

            all_r2 = s.zeros([self.dim['M'], self.dim['K']])
            for m in range(self.dim['M']):

                # Fetch the mask for missing vlaues
                mask = self.nodes["Y"].getNodes()[m].getMask()

                # Calculate predictions and mask them
                Ypred_m = s.dot(Z, W[m].T)
                Ypred_m[mask] = 0.

                # If there is an intercept term, regress it out, as it greatly decreases the fraction of variance explained by the other factors
                # (THIS IS NOT IDEAL...)
                if s.all(Z[:,0]==1.):
                    # calculate the total R2
                    SS = ((Y[m] - Y[m].mean())**2.).sum()
                    Res = ((Ypred_m - Y[m])**2.).sum()
                    R2_tot = 1. - Res/SS

                    # intercept_m = s.outer(Z[:,0], W[m][:,0].T)
                    # intercept_m[mask] = 0. # DO WE NEED TO DO THIS???
                    # Ypred_m -= intercept_m

                    all_r2[:,0] = 1.
                    for k in range(1,self.dim['K']):
                        # adding the intercept to the predictions with factor k
                        Ypred_mk = s.outer(Z[:,k], W[m][:,k]) + s.outer(Z[:,0], W[m][:,0])
                        Ypred_mk[mask] = 0.  # should not be necessary anymore as we do masked data - this ?
                        Res = ((Y[m] - Ypred_mk)**2.).sum()
                        R2_k = 1. - Res/SS
                        all_r2[m,k] = R2_k/R2_tot
                # No intercept term
                else:
                    # calculate the total R2
                    SS = ((Y[m] - Y[m].mean())**2.).sum()
                    Res = ((Ypred_m - Y[m])**2.).sum()
                    R2_tot = 1. - Res/SS

                    for k in range(self.dim['K']):
                        Ypred_mk = s.outer(Z[:,k], W[m][:,k])
                        Ypred_mk[mask] = 0.
                        Res = ((Y[m] - Ypred_mk)**2.).sum()
                        R2_k = 1. - Res/SS
                        all_r2[m,k] = R2_k/R2_tot

            if by_r2 is not None:
                drop_dic["by_r2"] = s.where( (all_r2>by_r2).sum(axis=0) == 0)[0]
                if len(drop_dic["by_r2"]) > 0:
                    drop_dic["by_r2"] = [ s.random.choice(drop_dic["by_r2"]) ]

        # Drop the factors
        drop = s.unique(s.concatenate(list(drop_dic.values())))
        if len(drop) > 0:
            for node in self.nodes.keys():
                self.nodes[node].removeFactors(drop)
        self.dim['K'] -= len(drop)

        if self.dim['K']==0:
            print("All factors shut down, no structure found in the data.")
            exit()

        pass

    def step_size(self, iter):
        # return the step size for the considered iterration
        tau = self.options['step_tau']
        kappa = self.options['forgetting_rate']
        return (iter + tau)**(-kappa)

    def sample_mini_batch(self):
        # TODO if multiple group, sample indices in each group evenly ? prob yes
        S_pc = self.options['batch_size']
        S = S_pc * self.dim['N']
        ix = s.random.choice(range(self.dim['N']), size= S, replace=False)
        return ix

    def iterate(self):
        """Method to start iterating and updating the variables using the VB algorithm"""

        # Define some variables to monitor training
        nodes = list(self.getVariationalNodes().keys())
        elbo = pd.DataFrame(data = nans((self.options['maxiter'], len(nodes)+1 )), columns = nodes+["total"] )
        activeK = nans((self.options['maxiter']))

        # Start training
        ro = None
        ix = None
        stochastic = True
        for i in range(self.options['maxiter']):
            if stochastic: # TODO
                ro = self.step_size(i)
                ix = self.sample_mini_batch()

            t = time();
            # Remove inactive latent variables
            if (i >= self.options["start_drop"]) and (i % self.options['freq_drop']) == 0:
                if any(self.options['drop'].values()):
                    self.removeInactiveFactors(**self.options['drop'])
                activeK[i] = self.dim["K"]

            # Update node by node, with E and M step merged
            for node in self.options['schedule']:
                # print "Node: " + str(node)
                # t = time()
                if (node=="ThetaW" or node=="ThetaZ") and i<self.options['start_sparsity']:
                    continue
                self.nodes[node].update(ix, ro)
                # print "time: " + str(time()-t)

            # Calculate Evidence Lower Bound
            if (i+1) % self.options['elbofreq'] == 0:
                elbo.iloc[i] = self.calculateELBO()

                # Print first iteration
                if i==0:
                    print("Iteration 1: time=%.2f ELBO=%.2f, Factors=%d, Covariates=%d" % (time() - t, elbo.iloc[i]["total"], (~self.nodes['Z'].covariates).sum(),self.nodes['Z'].covariates.sum()))
                    if self.options['verbose']:
                        print("".join([ "%s=%.2f  " % (k,v) for k,v in elbo.iloc[i].drop("total").iteritems() ]) + "\n")

                else:
                    # Check convergence using the ELBO
                    delta_elbo = elbo.iloc[i]["total"]-elbo.iloc[i-self.options['elbofreq']]["total"]

                    # Print ELBO monitoring
                    print("Iteration %d: time=%.2f ELBO=%.2f, deltaELBO=%.4f, Factors=%d, Covariates=%d" % (i+1, time()-t, elbo.iloc[i]["total"], delta_elbo, (~self.nodes['Z'].covariates).sum(), self.nodes['Z'].covariates.sum() ))
                    if delta_elbo<0 and not stochastic:
                        print("Warning, lower bound is decreasing..."); print('\a')
                        #import os; os.system('play --no-show-progress --null --channels 1 synth %s sine %f' % (0.01, 440))

                    if self.options['verbose']:
                        print("".join([ "%s=%.2f  " % (k,v) for k,v in elbo.iloc[i].drop("total").iteritems() ]) + "\n")

                    # Assess convergence
                    if (0 <= delta_elbo < self.options['tolerance']) and (not self.options['forceiter']):
                        activeK = activeK[:(i+1)]
                        elbo = elbo[:(i+1)]
                        print ("Converged!\n")
                        break

            # Do not calculate lower bound
            else:
                print("Iteration %d: time=%.2f, K=%d\n" % (i+1,time()-t,self.dim["K"]))

            # Flush (we need this to print when running on the cluster)
            sys.stdout.flush()

        # Finish by collecting the training statistics
        self.train_stats = { 'activeK':activeK, 'elbo':elbo["total"].values, 'elbo_terms':elbo.drop("total",1) }
        self.trained = True

    def getVariationalNodes(self):
        """ Method to return all variational nodes """
        # TODO problem with dictionnary comprehension here
        to_ret = {}
        for node in self.nodes.keys():
            if isinstance(self.nodes[node],Variational_Node):
                to_ret[node] =self.nodes[node]

        return to_ret
        # return { node:self.nodes[node] for node in self.nodes.keys() if isinstance(self.nodes[node],Variational_Node)}
        # return { k:v for k,v in self.nodes.items() if isinstance(v,Variational_Node) }

    def getTrainingStats(self):
        """ Method to return training statistics """
        return self.train_stats

    def getTrainingOpts(self):
        """ Method to return training options """
        return self.options

    def getTrainingData(self):
        """ Method to return training data """
        return self.nodes["Y"].getValues()

    def calculateELBO(self, *nodes):
        """Method to calculate the Evidence Lower Bound of the model"""
        if len(nodes) == 0: nodes = self.getVariationalNodes().keys()
        elbo = pd.Series(s.zeros(len(nodes)+1), index=list(nodes)+["total"])
        for node in nodes:
            elbo[node] = float(self.nodes[node].calculateELBO())
            elbo["total"] += elbo[node]
        return elbo
