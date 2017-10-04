"""
Script to test the model
"""

from __future__ import division
from time import time
import cPickle as pkl
import scipy as s
import os
import scipy.special as special
import scipy.stats as stats
import numpy.linalg  as linalg

# Import manually defined functions
from simulate import Simulate
from biofam.core.BayesNet import BayesNet
from biofam.core.nodes.multiview_nodes import *
from biofam.core.nodes.nongaussian_nodes import *
from biofam.core.nodes import *
from biofam.core.utils import *
from biofam.build_model.build_model import *

###################
## Generate data ##
###################

# import numpy; numpy.random.seed(4)

# Define dimensionalities
M = 3
# N = 3000
N = 100
D = s.asarray([1000,1000,1000])
# D = s.asarray([500,500,500])
# D = s.asarray([50,50,50])
K = 7


## Simulate data  ##
data = {}
tmp = Simulate(M=M, N=N, D=D, K=K)

# data['Z'] = s.zeros((N,K))
# data['Z'][:,0] = s.sin(s.arange(N)/(N/20))
# data['Z'][:,1] = s.cos(s.arange(N)/(N/20))
# data['Z'][:,2] = 2*(s.arange(N)/N-0.5)
# data['Z'][:,3] = stats.norm.rvs(loc=0, scale=1, size=N)
# data['Z'][:,4] = stats.norm.rvs(loc=0, scale=1, size=N)
# data['Z'][:,5] = stats.norm.rvs(loc=0, scale=1, size=N)
data['Z'] = stats.norm.rvs(loc=0, scale=1, size=(N,K))


data['alpha'] = [ s.zeros(K,) for m in xrange(M) ]
data['alpha'][0] = [1,1,1,1e6,1,1e6,1e6]
data['alpha'][1] = [1,1,1e6,1,1e6,1,1e6]
data['alpha'][2] = [1,1e6,1,1,1e6,1e6,1]
data['theta'] = [ s.ones((D[m],K))*1.0 for m in xrange(M) ]
# data['theta'] = [ s.repeat(s.linspace(0.1, 1.0, K)[None,:],D[m],0) for m in xrange(M) ]
data['S'], data['W'], data['W_hat'], _ = tmp.initW_spikeslab(theta=data['theta'], alpha=data['alpha']) #, annotation=False
data['mu'] = [ s.ones(D[m])*10. for m in xrange(M)]
data['tau']= [ stats.uniform.rvs(loc=1,scale=3,size=D[m]) for m in xrange(M) ]
# data['tau']= [ stats.uniform.rvs(loc=0.1,scale=3,size=D[m]) for m in xrange(M) ]

missingness = 0.0
missing_view = 0.10
# Y_warp = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
#   likelihood="warp", missingness=missingness, missing_view=missing_view)
Y_gaussian = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
	likelihood="gaussian", missingness=missingness, missing_view=missing_view)
# Y_poisson = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
# 	likelihood="poisson", missingness=missingness, missing_view=missing_view)
# Y_bernoulli = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
	# likelihood="bernoulli", missingness=missingness, missing_view=missing_view)
# Y_binomial = tmp.generateData(W=data['W'], Z=data['Z'], Tau=data['tau'], Mu=data['mu'],
# 	likelihood="binomial", min_trials=10, max_trials=50, missingness=missingness)

data["Y"] = [ Y_gaussian[0], Y_gaussian[1], Y_gaussian[2] ]
# data["Y"] = ( Y_gaussian[0], Y_bernoulli[1], Y_bernoulli[2] )
# data["Y"] = [ Y_bernoulli[0], Y_bernoulli[1], Y_bernoulli[2] ]
# data["Y"] = ( Y_gaussian[0], Y_poisson[1], Y_bernoulli[2] )
# data["Y"] = ( Y_warp[0], )
# data["Y"] = ( Y_warp[0], Y_warp[1], Y_warp[2] )
# data["Y"] = ( Y_gaussian[0], )

##################
## Data options ##
##################


data_opts = {}
data_opts["outfile"] = "/tmp/test.h5"
data_opts['view_names'] = ["View 1","View 2","View 3"]
data_opts['center'] = [True]*M
data_opts['scale_features'] = [True]*M
# data_opts['covariates'] = stats.bernoulli.rvs(0.5, size=N)[:,None]
data_opts['covariates'] = None
# data_opts['scale_covariates'] = [True]
data_opts['scale_covariates'] = None
data_opts['RemoveIncompleteSamples'] = True


#####################
## Data processing ##
#####################

for m in xrange(M):

  # Center the features
  if data_opts['center'][m]:
    data["Y"][m] = (data["Y"][m] - data["Y"][m].mean(axis=0))

  # Scale the features
  if data_opts['scale_features'][m]:
    data["Y"][m] = data["Y"][m] / s.std(data["Y"][m], axis=0)

#################################
## Initialise Bayesian Network ##
#################################

# Define initial number of latent variables
K = 15

# Define model dimensionalities
dim = {}
dim["M"] = M
dim["N"] = N
dim["D"] = D
dim["K"] = K


##############################
## Define the model options ##
##############################

model_opts = {}

# Define type of covariates
# model_opts['nonsparse_covariates'] = None
if data_opts['covariates'] is not None:
  dim["K"] += data_opts['covariates'].shape[1]
  K = dim["K"]

# Define whether to learn means
model_opts["learnMean"] = True
if model_opts["learnMean"]:
  if data_opts['covariates'] is None:
    data_opts['covariates'] = s.ones((dim["N"],1))
    data_opts['scale_covariates'] = [False]
  else:
    data_opts['covariates'] = s.concatenate((data_opts['covariates'],s.ones((dim["N"],1))), axis=1)
    data_opts['scale_covariates'].append(False)
  dim["K"] += data_opts['covariates'].shape[1]
  K = dim["K"]



# Define likelihoods
model_opts['likelihood'] = ['gaussian']* M
# model_opts['likelihood'] = ['gaussian','poisson','bernoulli']
# model_opts['likelihood'] = ['gaussian','bernoulli','bernoulli']
# model_opts['likelihood'] = ['bernoulli']*M

# Define initial number of factors
model_opts['k'] = K

# Define sparsities
model_opts['ardZ'] = False
model_opts['ardW'] = "mk"

# Define for which factors to learn Theta
model_opts['learnTheta'] = [0.*s.ones(K) for m in xrange(M)]

# Define schedule of updates
model_opts['schedule'] = ["Y","SW","Z","AlphaW","Theta","Tau"]

####################################
## Define priors (P distribution) ##
####################################

# Latent Variables
model_opts["priorZ"] = { 'mean':s.zeros((N,K)) }
if model_opts['ardZ']:
  model_opts["priorZ"]['var'] = s.ones((K,))*s.nan
  model_opts["priorAlphaZ"] = { 'a':s.ones(K)*1e-3, 'b':s.ones(K)*1e-3 }
else:
  model_opts["priorZ"]['var'] = s.ones((K,))*1.


# Weights
model_opts["priorSW"] = { 'Theta':[s.nan]*M, 'mean_S0':[s.nan]*M, 'var_S0':[s.nan]*M, 'mean_S1':[s.nan]*M, 'var_S1':[s.nan]*M } # Not required
if model_opts['ardW'] == "m":
  model_opts["priorAlphaW"] = { 'a':[1e-3]*M, 'b':[1e-3]*M }
elif model_opts['ardW'] == "k":
  model_opts["priorAlphaW"] = { 'a':s.ones(K)*1e-3, 'b':s.ones(K)*1e-3 }
elif model_opts['ardW'] == "mk":
  model_opts["priorAlphaW"] = { 'a':[s.ones(K)*1e-14]*M, 'b':[s.ones(K)*1e-14]*M }


# Theta
model_opts["priorTheta"] = { 'a':[s.ones(K,) for m in xrange(M)], 'b':[s.ones(K,) for m in xrange(M)] }
for m in xrange(M):
  for k in xrange(K):
    if model_opts['learnTheta'][m][k]==0:
      model_opts["priorTheta"]["a"][m][k] = s.nan
      model_opts["priorTheta"]["b"][m][k] = s.nan


# Noise
model_opts["priorTau"] = { 'a':[s.ones(D[m])*1e-14 for m in xrange(M)], 'b':[s.ones(D[m])*1e-14 for m in xrange(M)] }


#############################################
## Define initialisations (Q distribution) ##
#############################################

# Latent variables
# model_opts["initZ"] = { 'mean':"orthogonal", 'var':s.ones((N,K)), 'E':None, 'E2':None }
model_opts["initZ"] = { 'mean':"random", 'var':s.ones((K,)), 'E':None, 'E2':None }
if model_opts['ardZ']:
  model_opts["initAlphaZ"] = { 'a':s.nan, 'b':s.nan, 'E':s.ones(K)*100 }

# ARD of weights
if model_opts['ardW'] == "m":
  # model_opts["initAlphaW"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[1.]*M }
  model_opts["initAlphaW"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[ K*D[m]/(data[m].std(axis=0)**2 - 1./model_opts["initTau"]["E"][m]).sum() for m in xrange(M) ] }
elif model_opts['ardW'] == "k":
  model_opts["initAlphaW"] = { 'a':s.nan*s.ones(K), 'b':s.nan*s.ones(K), 'E':s.ones(K)*100 }
elif model_opts['ardW'] == "mk":
  model_opts["initAlphaW"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[100.*s.ones(K) for m in xrange(M)] }

# Theta
model_opts["initTheta"] = { 'a':[s.ones(K,) for m in xrange(M)], 'b':[s.ones(K,) for m in xrange(M)], 'E':[s.ones((D[m],K)) for m in xrange(M)] }
for m in xrange(M):
  for k in xrange(K):
    if model_opts['learnTheta'][m][k]==0.:
      model_opts["initTheta"]["a"][m][k] = s.nan
      model_opts["initTheta"]["b"][m][k] = s.nan


# Weights
model_opts["initSW"] = {
  'Theta':[s.ones((D[m],K)) for m in xrange(M)],
  'mean_S0':[s.zeros((D[m],K)) for m in xrange(M)],
  'var_S0':[s.nan*s.ones((D[m],K)) for m in xrange(M)],
  'mean_S1':[s.zeros((D[m],K)) for m in xrange(M)], # (TO-DO) allow also random
  'var_S1':[s.ones((D[m],K)) for m in xrange(M)],
  'ES':[None]*M, 'EW_S0':[None]*M, 'EW_S1':[None]*M }



# Noise
model_opts["initTau"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[s.ones(D[m])*100 for m in xrange(M)] }


######################################################
## Modify priors and initialisations for covariates ##
######################################################

if data_opts['covariates'] is not None:
  idx = xrange(data_opts['covariates'].shape[1])

  ## Prior distributions (P) ##
  # Latent variables
  if model_opts['ardZ']:
    model_opts["priorZ"]["mean"][:,idx] = s.nan
    model_opts["priorAlphaZ"]["a"][idx] = s.nan
    model_opts["priorAlphaZ"]["b"][idx] = s.nan
  else:
    model_opts["priorZ"]["var"][idx] = s.nan

  ## Variational distributions (Q) ##
  # Latent variables
  # model_opts["initZ"]["mean"][:,idx] = model_opts["covariates"]
  model_opts["initZ"]["var"][idx] = 0.
  if model_opts['ardZ']:
        model_opts["initAlphaZ"]["E"][idx] = s.nan

  ###########################
  ## Non-sparse covariates ##
  ###########################

  # if any(model_opts['nonsparse_covariates']):
  #   idx = s.where(model_opts['nonsparse_covariates'])
  #   model_opts['learnTheta'][:,idx] = 0.

  #   # WHERE IS THE LOOP OVER M?????

  #   ## Prior distributions (P) ##
  #   # Theta
  #   model_opts["priorTheta"]['a'][m][idx] = s.nan
  #   model_opts["priorTheta"]['b'][m][idx] = s.nan

  #   ## Variational distributions (Q) ##
  #   for m in range(M):
  #     # Weights
  #     model_opts["initSW"]["Theta"][m][:,idx] = 1.
  #     # Theta
  #     model_opts["initTheta"]["a"][m][idx] = s.nan
  #     model_opts["initTheta"]["b"][m][idx] = s.nan
  #     model_opts["initTheta"]["E"][m][idx] = 1.


if model_opts["learnMean"]:
  for m in range(M):

    # Weights
    if model_opts['likelihood'][m]=="gaussian":
      model_opts["initSW"]["mean_S1"][m][:,0] = data["Y"][m].mean(axis=0)
      model_opts["initSW"]["var_S1"][m][:,0] = 1e-5

    # Theta
    model_opts['learnTheta'][m][0] = 0.
    model_opts["initSW"]["Theta"][m][:,0] = 1.
    model_opts["priorTheta"]['a'][m][0] = s.nan
    model_opts["priorTheta"]['b'][m][0] = s.nan
    model_opts["initTheta"]["a"][m][0] = s.nan
    model_opts["initTheta"]["b"][m][0] = s.nan
    model_opts["initTheta"]["E"][m][:,0] = 1.


#############################
## Define training options ##
#############################

train_opts = {}
train_opts['elbofreq'] = 1
train_opts['maxiter'] = 1000
# train_opts['tolerance'] = 1E-2
train_opts['tolerance'] = 0.01
train_opts['forceiter'] = True
train_opts['drop'] = { "by_norm":None, "by_pvar":None, "by_cor":None, "by_r2":0.01 }
train_opts['startdrop'] = 3
train_opts['freqdrop'] = 1
train_opts['savefreq'] = s.nan
train_opts['savefolder'] = s.nan
train_opts['startSparsity'] = 500
train_opts['verbosity'] = 2
train_opts['verbose'] = True
train_opts['trials'] = 1
train_opts['cores'] = 1


####################
## Start training ##
####################

keep_best_run = False

runMultipleTrials(data["Y"], data_opts, model_opts, train_opts, keep_best_run)
# runSingleTrial(data["Y"], data_opts, model_opts, train_opts)
