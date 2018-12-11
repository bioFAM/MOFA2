
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
import numpy.ma as ma
import math
import resource

from biofam.core.nodes.variational_nodes import Variational_Node
from .utils import corr, nans, infer_platform

import warnings
warnings.filterwarnings("ignore")

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
        assert "gpu_mode" in train_opts, "'gpu_mode' not found in the training options dictionary"
        assert "stochastic" in train_opts, "'stochastic' not found in the training options dictionary"

        self.options = train_opts

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

    def calculate_variance_explained(self):

        # Collect relevant expectations
        Z = self.nodes['Z'].getExpectation()
        W = self.nodes["W"].getExpectation()
        Y = self.nodes["Y"].getExpectation()

        # Get groups
        groups = self.nodes["AlphaZ"].groups if "AlphaZ" in self.nodes else s.array([0]*self.dim['N'])

        r2 = [ s.zeros([self.dim['M'], self.dim['K']])] * self.dim['P']
        for m in range(self.dim['M']):
            mask = self.nodes["Y"].getNodes()[m].getMask()
            Ypred_m = s.dot(Z, W[m].T)
            Ypred_m[mask] = 0.
            for g in range(self.dim['P']):
                gg = groups==g
                SS = s.square(Y[m][gg,:]).sum()
                for k in range(self.dim['K']):
                    Ypred_mk = s.outer(Z[gg,k], W[m][:,k])
                    Ypred_mk[mask[gg,:]] = 0.
                    Res_k = ((Y[m][gg,:] - Ypred_mk)**2.).sum()
                    r2[g][m,k] = 1. - Res_k/SS
        return r2

    def calculate_total_variance_explained(self):

        # Collect relevant expectations
        Z = self.nodes['Z'].getExpectation()
        W = self.nodes["W"].getExpectation()
        Y = self.nodes["Y"].getExpectation()

        r2 = s.zeros(self.dim['M'])
        for m in range(self.dim['M']):
            mask = self.nodes["Y"].getNodes()[m].mask

            Ypred_m = s.dot(Z, W[m].T)
            Ypred_m[mask] = 0.

            Res = ((Y[m].data - Ypred_m)**2.).sum()
            SS = s.square(Y[m]).sum()

            r2[m] = 1. - Res/SS

        return r2

    def removeInactiveFactors(self, min_r2=None):
        """Method to remove inactive factors

        PARAMETERS
        ----------
        min_r2: float
            threshold to shut down factors based on a minimum variance explained per group and view
        """
        drop_dic = {}

        if min_r2 is not None:
            r2 = self.calculate_variance_explained()

            tmp = [ s.where( (r2[g]>min_r2).sum(axis=0) == 0)[0] for g in range(self.dim['P']) ]
            drop_dic["min_r2"] = list(set.intersection(*map(set,tmp)))
            if len(drop_dic["min_r2"]) > 0:
                drop_dic["min_r2"] = [ s.random.choice(drop_dic["min_r2"]) ]

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

    def step_size(self, i):
        # return the step size for the considered iterration
        tau = self.options['tau']
        kappa = self.options['forgetting_rate']
        return (i + tau)**(-kappa)

    def step_size2(self, i):
        # return the step size for the considered iterration
        tau = self.options['tau']
        kappa = self.options['forgetting_rate']
        return tau / ((1 + kappa * i)**(3./4.))

    def sample_mini_batch_replace(self):
        # TODO if multiple group, sample indices in each group evenly ? prob yes
        S_pc = self.options['batch_size']
        S = S_pc * self.dim['N']
        ix = s.random.choice(range(self.dim['N']), size=S, replace=False)
        return ix

    def sample_mini_batch_no_replace(self, i):
        """ Method to define mini batches"""

        # TODO :
        # - if multiple group, sample indices in each group evenly ? prob yes
        # - ??? shuffle the data at the beginnign of every epoch

        # Sample mini-batch indices and define epoch
        n_batches = math.ceil(1./self.options['batch_size'])
        S = self.options['batch_size'] * self.dim['N']
        batch_ix = i % n_batches
        epoch = int(i / n_batches)
        if batch_ix == 0:
            print("## Epoch %d ##" % epoch)
            print("-------------------------------------------------------------------------------------------")
            self.shuffled_ix = s.random.choice(range(self.dim['N']), size= self.dim['N'], replace=False)

        min = int(S * batch_ix)
        max = int(S * (batch_ix + 1))
        if max > self.dim['N']:
            max = self.dim['N']

        ix = self.shuffled_ix[min:max]

        # Define mini-batch for each node
        self.nodes['Y'].define_mini_batch(ix)
        self.nodes['Tau'].define_mini_batch(ix)
        if 'AlphaZ' in self.nodes:
            self.nodes['AlphaZ'].define_mini_batch(ix)
        if 'ThetaZ' in self.nodes:
            self.nodes['ThetaZ'].define_mini_batch(ix)  

        return ix, epoch
    
    def iterate(self):
        """Method to start iterating and updating the variables using the VB algorithm"""

        # Define some variables to monitor training
        nodes = list(self.getVariationalNodes().keys())
        elbo = pd.DataFrame(data = nans((self.options['maxiter'], len(nodes)+1 )), columns = nodes+["total"] )
        activeK = nans((self.options['maxiter']))

        # Precompute terms
        for n in self.nodes:
            self.nodes[n].precompute(self.options)

        # Print statistics before training
        if self.options['verbose']: 
            foo = self.calculateELBO()
            r2 = self.calculate_total_variance_explained()

            print('ELBO before training:')
            print("".join([ "%s=%.2f  " % (k,v) for k,v in foo.drop("total").iteritems() ]) + "\nTotal: %.2f\n" % foo["total"])
            print('Schedule of updates: ',self.options['schedule']);
            print("Variance explained:\t" + "   ".join([ "View %s: %.3f%%" % (m,100*r2[m]) for m in range(self.dim["M"])]))
            if self.options['stochastic']:
                print("Using stochastic variational inference with the following parameters:")
                print("- Batch size (fraction of samples): %.2f\n- Forgetting rate: %.2f\n" % (100*self.options['batch_size'], self.options['forgetting_rate']) )
            print("\n")

        ro = 1.
        ix = None
        for i in range(self.options['maxiter']):
            t = time();

            # Sample mini-batch and define step size for stochastic inference
            if self.options['stochastic'] and (i >= self.options["start_stochastic"]-1):
                ix, epoch = self.sample_mini_batch_no_replace(i)
                ro = self.step_size2(epoch)

            # Remove inactive latent variables
            if (i >= self.options["start_drop"]) and (i % self.options['freq_drop']) == 0:
                # if any(self.options['drop'].values()):
                if self.options['drop']["min_r2"] is not None:
                    self.removeInactiveFactors(**self.options['drop'])
                activeK[i] = self.dim["K"]

            # Update node by node, with E and M step merged
            t_updates = time()
            for node in self.options['schedule']:
                if (node=="ThetaW" or node=="ThetaZ") and i<self.options['start_sparsity']:
                    continue
                self.nodes[node].update(ix, ro)
            t_updates = time() - t_updates

            # Calculate Evidence Lower Bound
            if (i+1) % self.options['elbofreq'] == 0:
                t_elbo = time()
                elbo.iloc[i] = self.calculateELBO()
                t_elbo = time() - t_elbo

                # Print first iteration
                if i==0:
                    print("Iteration 1: time=%.2f, ELBO=%.2f, Factors=%d" % (time() - t, elbo.iloc[i]["total"], (self.dim['K'])))
                    if self.options['verbose']:
                        print("".join([ "%s=%.2f  " % (k,v) for k,v in elbo.iloc[i].drop("total").iteritems() ]))
                else:
                    # Check convergence using the ELBO
                    delta_elbo = elbo.iloc[i]["total"]-elbo.iloc[i-self.options['elbofreq']]["total"]

                    # Print ELBO monitoring
                    print("Iteration %d: time=%.2f, ELBO=%.2f, deltaELBO=%.4f, Factors=%d" % (i+1, time()-t, elbo.iloc[i]["total"], delta_elbo, (self.dim['K'])))
                    if delta_elbo<0 and not self.options['stochastic']: print("Warning, lower bound is decreasing...\a")

                    # Print ELBO decomposed by node and variance explained
                    if self.options['verbose']:
                        print("".join([ "%s=%.2f  " % (k,v) for k,v in elbo.iloc[i].drop("total").iteritems() ]))
                        print('Time spent in ELBO computation: %.1f%%' % (100*t_elbo/(t_updates+t_elbo)) )

                    # Assess convergence
                    if (abs(delta_elbo) < self.options['tolerance']) and (not self.options['forceiter']):
                        activeK = activeK[:(i+1)]; elbo = elbo[:(i+1)]
                        print ("Converged!\n"); break

            # Do not calculate lower bound
            else:
                print("Iteration %d: time=%.2f, Factors=%d" % (i+1,time()-t,self.dim["K"]))

            # Print other statistics
            if self.options['verbose']:
                # Memory usage
                print('Peak memory usage: %.2f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / infer_platform() ))
                # Variance explained
                r2 = self.calculate_total_variance_explained()
                print("Variance explained:\t" + "   ".join([ "View %s: %.3f%%" % (m,100*r2[m]) for m in range(self.dim["M"])]))
                # Sparsity levels of the weights
                W = self.nodes["W"].getExpectation()
                foo = [s.mean(s.absolute(W[m])<1e-3) for m in range(self.dim["M"])]
                print("Fraction of zero weights:\t" + "   ".join([ "View %s: %.0f%%" % (m,100*foo[m]) for m in range(self.dim["M"])]))
                # Sparsity levels of the factors
                Z = self.nodes["Z"].getExpectation()
                bar = s.mean(s.absolute(Z)<1e-3)
                print("Fraction of zero samples: %.0f%%" % (100*bar))
                # stochastic inference 
                if self.options['stochastic'] and (i >= self.options["start_stochastic"]-1): 
                    print("Step size: %.4f" % ro)

            # Flush (we need this to print when running on the cluster)
            print("\n")
            sys.stdout.flush()

        # Finish by collecting the training statistics
        self.train_stats = { 'activeK':activeK, 'elbo':elbo["total"].values, 'elbo_terms':elbo.drop("total",1) }
        self.trained = True

    def compute_r2_simple(self):
        # compute r2 for the cnosidered mini bact
        # ----------------------------------------------------------------------
        W = s.concatenate(self.nodes['W'].getExpectation())
        Z = self.nodes['Z'].get_mini_batch()['E']
        Y = s.concatenate(self.nodes['Y'].get_mini_batch(), axis=1)

        Y_mask = ma.getmask(Y)
        Y_dat = Y.data
        Y_dat[Y_mask] = 0.

        pred = Z.dot(W.T)
        pred[Y_mask] = 0.
        SS = s.sum((Y_dat - pred)**2.)
        var = s.sum((Y_dat - Y_dat.mean())**2.)

        r2_batch = 1. - SS/var

        # compute r2 for all data
        # ----------------------------------------------------------------------
        W = s.concatenate(self.nodes['W'].getExpectation())
        Z = self.nodes['Z'].getExpectation()
        Y = s.concatenate(self.nodes['Y'].getExpectation(), axis=1)

        Y_mask = ma.getmask(Y)
        Y_dat = Y.data
        Y_dat[Y_mask] = 0.

        pred = Z.dot(W.T)
        pred[Y_mask] = 0.
        SS = s.sum((Y_dat - pred)**2.)
        var = s.sum((Y_dat - Y_dat.mean())**2.)

        r2_tot = 1. - SS/var

        # print
        # ----------------------------------------------------------------------
        print("batch specific r2 is ", r2_batch)
        print("total r2 is ", r2_tot)
        print()

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

