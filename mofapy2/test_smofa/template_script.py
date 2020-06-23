# script to run SMOFA from python

from smofapy.run.entry_point import entry_point

import pandas as pd
import scipy as s
import numpy as np
import io
import requests # to download the online data
import h5py
import scipy.stats as stats
import matplotlib.pyplot as plt
import math
import os
import scipy.spatial as SS
from smofapy.test import simulate_smofa as simsmofa

###########################
## Simulate data ##
###########################
nfactors = 2
lscales = [0.1, 0.3]
N = 1000
noise_level = 1
Dm = 500
missing= 0
seed = 41284
sim = simsmofa.simulate_data(N=N, seed=seed, views=["0", "1", "2", "3"], D=[Dm] * 4,
                            noise_level=noise_level,
                            K=int(nfactors), lscales=lscales)
sim['data'] = simsmofa.mask_data(sim, perc = missing)
data = sim['data']
sample_cov = sim['sample_cov'].reshape(N,1)


###########################
## Initialise SMOFA model ##
###########################

# initialise the entry point
ent = entry_point()

# Set data options
# - scale_views: if views have significantly different ranges/variances, it is good practice to scale each view to unit variance
ent.set_data_options(scale_views=False)

## Set data ##
#(option 1, list of matrices)
ent.set_data_matrix(data, sample_cov)

#(option 2, data.frame)
#ent.set_data_df(data)

## Set model options ##
# - factors: number of factors. Default is K=10
# - likelihods: likelihoods per view (options are "gaussian","poisson","bernoulli").
# 		Default is None, and they are guessed from the data
# - spikeslab_weights: use spike-slab sparsity prior in the weights? (recommended TRUE)
# - ard_weights: use ARD prior in the weights? (TRUE if using multiple views)
# - GP_factors: use a GP prior on the factors (default TRUE)
# - start_opt: at which iteration to start optimizing the lengthscales of the GP priors
# - n_grid: how many grid points for lengthscale optimization
# - mv_Znode: use a multivariate Z node? (default True)
# - smooth_all: only use smooth factors (default False)

ent.set_model_options(factors = 5, start_opt = 20, n_grid =50, mv_Znode = True)

## (Optional) Set stochastic inference options (not in conjuction with GP prior) ##
# - batch_size: float value indicating the batch size (as a fraction of the total data set: 0.10, 0.25 or 0.50)
# - learning_rate: learning rate (we recommend values from 0.25 to 0.75)
# - forgetting_rate: forgetting rate (we recommend values from 0.25 to 0.5)
# - start_stochastic: first iteration to apply stochastic inference (recommended > 5)
#
#ent.set_stochastic_options(batch_size=0.5, learning_rate=0.75, forgetting_rate=0.5, start_stochastic=10)

## (Optional) Set sparse GP inference options (in conjuction with GP prior) ##
# - n_inducing: number of inducing points
# - idx_inducing: optional argument to specify with points to use as inducing points (as index of original variables)
#
ent.set_sparseGP_options(n_inducing=100, idx_inducing = None) #TODO assert that this is called before training options

## Set training options ##
# - iter: number of iterations
# - convergence_mode: "fast", "medium", "slow".
#		For exploration, the fast mode is good enough.
# - startELBO: initial iteration to compute the ELBO (the objective function used to assess convergence)
# - freqELBO: frequency of computations of the ELBO (the objective function used to assess convergence)
# - dropR2: minimum variance explained criteria to drop factors while training.
# 		Default is None, inactive factors are not dropped during training
# - gpu_mode: use GPU mode?
#		if TRUE, this needs cupy installed and a functional GPU, see https://cupy.chainer.org/
# - verbose: verbose mode?
# - seed: random seed
# - weight_views: weight ELBO to cope with imbalance

ent.set_train_options(seed = 1234567, convergence_mode="fast", iter = 1000,
                      dropR2 = 0.01, weight_views = True)


####################################
## Build and train the MOFA model ##
####################################

# Build the model
ent.build()

# Run the model
ent.run()

##################################################################
## (Optional) do dimensionality reduction from the MOFA factors ##
##################################################################

# ent.umap()
# ent.tsne()

####################
## Save the model ##
####################

outfile="test_SMOFA.hdf5"
ent.save(outfile)


# Model analysis
# def explore(ent, sim):
#     # plot inferred factors
#     Zexp = ent.model.getNodes()['Z'].getExpectation()
#     #    Zexp = (Zexp - Zexp.mean(axis=0)) / Zexp.std(axis=0)
#     for i in range(Zexp.shape[1]):
#         plt.scatter(sample_cov, Zexp[:, i])
#         plt.title("Inferred factors")
#         plt.ylim(-4, 4)
#     plt.show()
#
#     # plot simulated factors
#     Z = sim['Z']
#     #    Z = (Z - Z.mean(axis=0)) / Z.std(axis=0)
#     for i in range(Z.shape[1]):
#         plt.scatter(sample_cov, Z[:, i])
#         plt.title("simulated factors")
#         plt.ylim(-4, 4)
#     plt.show()
#
#     # plot ELBO terms
#     plt.plot(ent.model.getTrainingStats()['elbo_terms'])
#     plt.legend(ent.model.getTrainingStats()['elbo_terms'].columns)
#     plt.title("ELBO")
#     plt.show()
#
#     elbos = ent.model.getTrainingStats()['elbo_terms']
#     for n in elbos.columns:
#         plt.plot(elbos[n])
#         plt.title("ELBO_term" + n)
#         plt.show()
#
#     # plot length scales
#     if 'Sigma' in ent.model.nodes.keys():
#         plt.plot(ent.model.getTrainingStats()['length_scales'])
#         plt.legend(ent.model.getTrainingStats()['length_scales'].columns)
#         plt.title("Lengthscales")
#         plt.show()
#
#         # OUTPUT ONLY FOR SMOFA
#         plt.plot(ent.model.getTrainingStats()['W_scales'])
#         plt.legend(ent.model.getTrainingStats()['W_scales'].columns)
#         plt.title("W-scales")
#         plt.show()
#
#         plt.plot(ent.model.getTrainingStats()['Z_scales'])
#         plt.legend(ent.model.getTrainingStats()['Z_scales'].columns)
#         plt.title("Z-scales")
#         plt.show()
#
#         plt.plot(ent.model.getTrainingStats()['alpha'])
#         plt.legend(ent.model.getTrainingStats()['alpha'].columns)
#         plt.title("alpha")
#         plt.show()
#
#         plt.plot(ent.model.getTrainingStats()['Sigma_inv_diag'])
#         plt.legend(ent.model.getTrainingStats()['Sigma_inv_diag'].columns)
#         plt.title("Elemnts on diagonal of Sigma_Inverse")
#         plt.show()
#
#         plt.plot(ent.model.getTrainingStats()['Sigma_det'])
#         plt.legend(ent.model.getTrainingStats()['Sigma_det'].columns)
#         plt.title("Determinant of Sigma")
#         plt.show()
#
#         plt.plot(ent.model.getTrainingStats()['offDiagaonal_Sigma_inv'])
#         plt.legend(ent.model.getTrainingStats()['offDiagaonal_Sigma_inv'].columns)
#         plt.title("Mean of absolute values of off-diagaonal elements of Sigma_inv")
#         plt.show()


# # helper function to save data
# def save_simulated_data(outfile):
#     filex = h5py.File(outfile,'w')
#     for m in range(M):
#         filex.create_dataset("W"+views[m], data=W[m], compression="gzip")
#     filex.create_dataset("Z", data=Z, compression="gzip")
#     filex.create_dataset("Sigma", data=Sigma, compression="gzip")
#     filex.create_dataset("lengthscale", data=lscales, compression="gzip")
#     filex.close()