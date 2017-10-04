
"""
This module is used to define the class containing the entire Bayesian Network,
and the corresponding attributes/methods to train the model, set algorithmic options, calculate lower bound, etc.

To-do:
- default values for schedule and options
"""

from __future__ import division
from time import time
import os
import scipy as s
import pandas as pd
import sys

from biofam.nodes.variational_nodes import Variational_Node
from .utils import corr, nans



class BayesNet(object):
    def __init__(self, dim, nodes, schedule, options, trial=1):
        """ Initialisation of a Bayesian network

        PARAMETERS
        ----------
        dim: dict
            keyworded dimensionalities, ex. {'N'=10, 'M'=3, ...}
        nodes: dict
            dictionary with all nodes where the keys are the name of the node and the values are instances the 'Node' class
        schedule: iterable
            list or tuple with the names of the nodes to be updated in the given order. Nodes not present in schedule will not be updated
        options: dict
            training options, such as maximum number of iterations, training options, etc.
        trial: int
            this is an auxiliary variable for parallelised running of multiple trials
        """

        self.dim = dim
        self.nodes = nodes
        self.schedule = schedule
        self.options = options
        self.trial = trial

        # Training flag
        self.trained = False

    def removeInactiveFactors(self, by_norm=None, by_pvar=None, by_cor=None, by_r2=None):
        """Method to remove inactive factors

        PARAMETERS
        ----------
        by_norm: float
            threshold to shut down factors based on the norm of the latent variable
            CURRENTLY NOT IMPLEMENTED
        by_pvar: float
            threshold to shut down factors based on the proportion of variance explained
            CURRENTLY NOT IMPLEMENTED
        by_cor: float
            threshold to shut down factors based on the correlation between latent variables
            CURRENTLY NOT IMPLEMENTED
        by_r2: float
            threshold to shut down factors based on the coefficient of determination
        """
        drop_dic = {}

        # Shut down based on norm of latent variable vectors
        #   Advantages: independent of likelihood type, works with pseudodata
        #   Disadvantages: it does not take into account the weights, covariates are never removed.
        # if by_norm is not None:
        #     Z = self.nodes["Z"].getExpectation()
        #     # Z = Z + 1e-6*stats.norm.rvs(loc=0, scale=1, size=(Z.shape[0],Z.shape[1])) # Add some noise to remove structure (see XX)
        #     drop_dic["by_norm"] = s.where((Z**2).mean(axis=0) < by_norm)[0]
        #     if len(drop_dic["by_norm"]) > 0:
        #         drop_dic["by_norm"] = [ s.random.choice(drop_dic["by_norm"]) ]

        # Shut down based on coefficient of determination with respect to the residual variance
        #   Advantages: it takes into account both weights and latent variables, is based on how well the model fits the data
        #   Disadvantages: slow, doesnt work with non-gaussian data
        if by_r2 is not None:
            Z = self.nodes['Z'].getExpectation()
            Y = self.nodes["Y"].getExpectation()
            W = self.nodes["SW"].getExpectation()
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
                    Ypred_m_intercept = s.outer(Z[:,0], W[m][:,0].T) 
                    Ypred_m_intercept[mask] = 0. # DO WE NEED TO DO THIS???
                    Ypred_m -= Ypred_m_intercept
                    all_r2[:,0] = 1.
                    SS = (Ypred_m**2.).sum()
                    for k in range(1,self.dim['K']):
                        Ypred_mk = s.outer(Z[:,k], W[m][:,k])
                        Ypred_mk[mask] = 0.
                        Res = ((Ypred_m - Ypred_mk)**2.).sum()
                        all_r2[m,k] = 1. - Res/SS
                # No intercept term
                else:
                    SS = (Ypred_m**2.).sum()
                    for k in range(self.dim['K']):
                        Ypred_mk = s.outer(Z[:,k], W[m][:,k])
                        Ypred_mk[mask] = 0.
                        Res = ((Ypred_m - Ypred_mk)**2.).sum()
                        all_r2[m,k] = 1. - Res/SS

            if by_r2 is not None:
                drop_dic["by_r2"] = s.where( (all_r2>by_r2).sum(axis=0) == 0)[0]
                if len(drop_dic["by_r2"]) > 0:
                    drop_dic["by_r2"] = [ s.random.choice(drop_dic["by_r2"]) ]

        # Shut down based on the proportion of residual variance explained by each factor
        # IT DOESNT WORK, THERE IS SOME ERROR TO BE FIXED
        #   Good: it is the proper way of doing it,
        #   Bad: slow, does it work with non-gaussian data?
        # if by_pvar is not None:
        #     Z = self.nodes["Z"].getExpectation()
        #     Y = self.nodes["Y"].getExpectation()
        #     tau = self.nodes["Tau"].getExpectation()
        #     alpha = self.nodes["Alpha"].getExpectation()

        #     factor_pvar = s.zeros((self.dim['M'],self.dim['K']))
        #     for m in range(self.dim['M']):
        #         residual_var = (s.var(Y[m],axis=0) - 1/tau[m]).sum()
        #         for k in range(self.dim["K"]):
        #             factor_var = (self.dim["D"][m]/alpha[m][k])# * s.var(Z[:,k])
        #             factor_pvar[m,k] = factor_var / residual_var
        #     drop_dic['by_pvar'] = s.where( (factor_pvar>by_pvar).sum(axis=0) == 0)[0]

        # Shut down factors that are highly correlated
        # (Q) Which of the two factors should we remove? Maybe the one that explains less variation
        # if by_cor is not None:
        #     Z = self.nodes["Z"].getExpectation()
        #     r = s.absolute(corr(Z.T,Z.T))
        #     s.fill_diagonal(r,0)
        #     r *= s.tri(*r.shape)
        #     drop_dic["by_cor"] = s.where(r>by_cor)[0]
        #     if len(drop_dic["by_cor"]) > 0:
        #         # Drop just one latent variable, chosen randomly
        #         drop_dic["by_cor"] = [ s.random.choice(drop_dic["by_cor"]) ]

        # Drop the factors
        drop = s.unique(s.concatenate(drop_dic.values()))
        if len(drop) > 0:
            for node in self.nodes.keys():
                self.nodes[node].removeFactors(drop)
        self.dim['K'] -= len(drop)

        if self.dim['K']==0:
            print("Shut down all components, no structure found in the data.")
            exit()

        pass

    def iterate(self):
        """Method to start iterating and updating the variables using the VB algorithm"""

        # Define some variables to monitor training
        nodes = list(self.getVariationalNodes().keys())
        elbo = pd.DataFrame(data = nans((self.options['maxiter'], len(nodes)+1 )), columns = nodes+["total"] )
        activeK = nans((self.options['maxiter']))
        
        # Start training
        for i in range(self.options['maxiter']):
            t = time();

            # Remove inactive latent variables
            if (i >= self.options["startdrop"]) and (i % self.options['freqdrop']) == 0:
                if any(self.options['drop'].values()):
                    self.removeInactiveFactors(**self.options['drop'])
                activeK[i] = self.dim["K"]

            # Update node by node, with E and M step merged
            for node in self.schedule:
                # print "Node: " + str(node)
                # t = time()
                if node=="Theta" and i<self.options['startSparsity']:
                    continue
                self.nodes[node].update()
                # print "time: " + str(time()-t)

            # Calculate Evidence Lower Bound
            if (i+1) % self.options['elbofreq'] == 0:
                elbo.iloc[i] = self.calculateELBO()

                # Print first iteration
                if i==0:
                    print("Trial %d, Iteration 1: time=%.2f ELBO=%.2f, Factors=%d, Covariates=%d" % (self.trial, time()-t,elbo.iloc[i]["total"], (~self.nodes["Z"].covariates).sum(), self.nodes["Z"].covariates.sum() ))
                    if self.options['verbose']:
                        print("".join([ "%s=%.2f  " % (k,v) for k,v in elbo.iloc[i].drop("total").iteritems() ]) + "\n")

                else:
                    # Check convergence using the ELBO
                    delta_elbo = elbo.iloc[i]["total"]-elbo.iloc[i-self.options['elbofreq']]["total"]

                    # Print ELBO monitoring
                    print("Trial %d, Iteration %d: time=%.2f ELBO=%.2f, deltaELBO=%.4f, Factors=%d, Covariates=%d" % (self.trial, i+1, time()-t, elbo.iloc[i]["total"], delta_elbo, (~self.nodes["Z"].covariates).sum(), self.nodes["Z"].covariates.sum() ))
                    if delta_elbo<0: print("Warning, lower bound is decreasing..."); print('\a')
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

    def getParameters(self, *nodes):
        """Method to collect all parameters of a given set of nodes (all by default)

        PARAMETERS
        ----------
        nodes: iterable
            name of the nodes
        """

        if len(nodes) == 0: nodes = self.nodes.keys()
        params = {}
        for node in nodes:
            tmp = self.nodes[node].getParameters()
            if tmp != None: params[node] = tmp
        return params

    def getExpectations(self, only_first_moments=False, *nodes):
        """Method to collect all expectations of a given set of nodes (all by default)
        
        PARAMETERS
        ----------
        only_first_moments: bool
            get only first moments?
        nodes: list
            name of the nodes
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

    def getVariationalNodes(self):
        """ Method to return all variational nodes """
        return { k:v for k,v in self.nodes.items() if isinstance(v,Variational_Node) }

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
