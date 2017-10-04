import argparse
import pandas as pd
import scipy as s
from time import sleep

from .build_model import *

def entry_point():

  banner = """
  ###########################################################
  ###                 __  __  ___  _____ _                ### 
  ###                |  \/  |/ _ \|  ___/ \               ### 
  ###                | |\/| | | | | |_ / _ \              ### 
  ###                | |  | | |_| |  _/ ___ \             ### 
  ###                |_|  |_|\___/|_|/_/   \_\            ### 
  ###                                                     ###
  ########################################################### """

  print(banner)
  sleep(2)

  # Read arguments
  p = argparse.ArgumentParser( description='Basic run script for BioFAM' )

  # I/O
  p.add_argument( '--inFiles',           type=str, nargs='+', required=True,                  help='Input data files (including extension)' )
  p.add_argument( '--outFile',           type=str, required=True,                             help='Output data file (hdf5 format)' )
  p.add_argument( '--delimiter',         type=str, default=" ",                               help='Delimiter for input files' )
  p.add_argument( '--covariatesFile',    type=str, default=None,                               help='Input data file for covariates' )
  p.add_argument( '--header_cols',       action='store_true',                                 help='Do the input files contain column names?' )
  p.add_argument( '--header_rows',       action='store_true',                                 help='Do the input files contain row names?' )

  # Data options
  p.add_argument( '--center_features',   action="store_true",                                 help='Center the features to zero-mean?' )
  p.add_argument( '--scale_features',    action="store_true",                                 help='Scale the features to unit variance?' )
  p.add_argument( '--scale_views',       action="store_true",                                 help='Scale the views to unit variance?' )
  p.add_argument( '--scale_covariates',  type=int, nargs='+', default=0,                      help='' )
  p.add_argument( '--maskAtRandom',      type=float,nargs="+", default=None,                  help='Fraction of data to mask per view')
  p.add_argument( '--maskNSamples',      type=int,nargs="+", default=None,                    help='Number of patients to mask per view')
  p.add_argument( '--RemoveIncompleteSamples', action="store_true",                           help='Remove samples with incomplete views?' )

  # Model options
  p.add_argument( '--factors',           type=int, default=10,                                help='Initial number of latent variables')
  p.add_argument( '--likelihoods',       type=str, nargs='+', required=True,                  help='Likelihood per view, current options are bernoulli, gaussian, poisson')
  p.add_argument( '--views',             type=str, nargs='+', required=True,                  help='View names')
  p.add_argument( '--schedule',          type=str, nargs="+", default=None,                   help='Update schedule, default is ( Y SW Z AlphaW Theta Tau )' )
  p.add_argument( '--learnTheta',        type=int, nargs="+", default=1,                      help='Learn the sparsity parameter from the spike-and-slab (theta)?' )
  p.add_argument( '--initTheta',         type=float, nargs="+", default=1. ,                  help='Initialisation for the sparsity parameter of the spike-and-slab (theta)')
  p.add_argument( '--learnIntercept',    action='store_true',                                 help='Learn the feature-wise mean?' )

  # Training options
  p.add_argument( '--elbofreq',          type=int, default=1,                                 help='Frequency of computation of ELBO' )
  p.add_argument( '--iter',              type=int, default=5000,                              help='Maximum number of iterations' )
  p.add_argument( '--ntrials',           type=int, default=1,                                 help='Number of trials' )
  p.add_argument( '--startSparsity',     type=int, default=100,                               help='Iteration to activate the spike-and-slab')
  p.add_argument( '--tolerance',         type=float, default=0.01 ,                            help='Tolerance for convergence (based on the change in ELBO)')
  p.add_argument( '--startDrop',         type=int, default=1 ,                                help='First iteration to start dropping factors')
  p.add_argument( '--freqDrop',          type=int, default=1 ,                                help='Frequency for dropping factors')
  p.add_argument( '--dropR2',            type=float, default=None ,                           help='Threshold to drop latent variables based on coefficient of determination' )
  p.add_argument( '--nostop',            action='store_true',                                 help='Do not stop when convergence criterion is met' )
  p.add_argument( '--verbose',           action='store_true',                                 help='Use more detailed log messages?')
  p.add_argument( '--seed',              type=int, default=0 ,                                help='Random seed' )


  args = p.parse_args()

  #############################
  ## Define the data options ##
  #############################

  data_opts = {}

  # I/O
  data_opts['input_files'] = args.inFiles
  data_opts['outfile'] = args.outFile
  data_opts['delimiter'] = args.delimiter

  # View names
  data_opts['view_names'] = args.views

  # Headers
  if args.header_rows:
    data_opts['rownames'] = 0
  else:
    data_opts['rownames'] = None

  if args.header_cols:
    data_opts['colnames'] = 0
  else:
    data_opts['colnames'] = None

  #####################
  ## Data processing ##
  #####################

  M = len(data_opts['input_files'])
  assert M == len(data_opts['view_names']), "Length of view names and input files does not match"

  # Data processing: center features
  if args.center_features:
    data_opts['center_features'] = [ True if l=="gaussian" else False for l in args.likelihoods ]
  else:
    if not args.learnIntercept: print("\nWarning... you are not centering the data and not learning the mean...\n")
    data_opts['center_features'] = [ False for l in args.likelihoods ]

  # Data processing: scale views
  if args.scale_views:
    data_opts['scale_views'] = [ True if l=="gaussian" else False for l in args.likelihoods ]
  else:
    data_opts['scale_views'] = [ False for l in args.likelihoods ]

  # Data processing: scale features
  if args.scale_features:
    assert args.scale_views==False, "Scale either entire views or features, not both"
    data_opts['scale_features'] = [ True if l=="gaussian" else False for l in args.likelihoods ]
  else:
    data_opts['scale_features'] = [ False for l in args.likelihoods ]


  # Data processing: mask values
  if args.maskAtRandom is not None:
    data_opts['maskAtRandom'] = args.maskAtRandom
    assert M == len(data_opts['maskAtRandom']), "Length of MaskAtRandom and input files does not match"
  else:
    data_opts['maskAtRandom'] = [0]*M

  if args.maskNSamples is not None:
    data_opts['maskNSamples'] = args.maskNSamples
    assert M == len(data_opts['maskNSamples']), "Length of MaskAtRandom and input files does not match"
  else:
    data_opts['maskNSamples'] = [0]*M

  # Remove incomplete samples?
  data_opts['RemoveIncompleteSamples'] = args.RemoveIncompleteSamples

  ###############
  ## Load data ##
  ###############

  # Load observations
  data = loadData(data_opts)

  # Remove samples with missing views
  if data_opts['RemoveIncompleteSamples']:
    data = removeIncompleteSamples(data)

  # Calculate dimensionalities
  N = data[0].shape[0]
  D = [data[m].shape[1] for m in range(M)]

  # Load covariates
  if args.covariatesFile is not None:
    data_opts['covariates'] = pd.read_csv(args.covariatesFile, delimiter=" ", header=None).as_matrix()
    print("Loaded covariates from " + args.covariatesFile + "with shape " + str(data_opts['covariates'].shape) + "...")
    data_opts['scale_covariates'] = args.scale_covariates
    if len(data_opts['scale_covariates']) == 1 and data_opts['covariates'].shape[1] > 1:
      data_opts['scale_covariates'] = args.scale_covariates[0] * s.ones(data_opts['covariates'].shape[1])
    elif type(data_opts['scale_covariates'])==list:
      assert len(data_opts['scale_covariates']) == data_opts['covariates'].shape[1], "'scale_covariates' has to be the same length as the number of covariates"
    data_opts['scale_covariates'] = [ bool(x) for x in data_opts['scale_covariates'] ]
    args.factors += data_opts['covariates'].shape[1]
  else:
    data_opts['scale_covariates'] = False
    data_opts['covariates'] = None

  # If we want to learn the mean, we add a constant covariate of 1s
  if args.learnIntercept:
    if data_opts['covariates'] is not None:
      # data_opts['covariates'].insert(0, "mean", s.ones(N,))
      data_opts['covariates'] = s.insert(data_opts['covariates'], obj=0, values=1, axis=1)
      data_opts['scale_covariates'].insert(0,False)
    else:
      data_opts['covariates'] = s.ones((N,1))
      data_opts['scale_covariates'] = [False]
    args.factors += 1

  ##############################
  ## Define the model options ##
  ##############################

  model_opts = {}

  # Define initial number of latent factors
  K = model_opts['k'] = args.factors

  # Define likelihoods
  model_opts['likelihood'] = args.likelihoods
  assert M==len(model_opts['likelihood']), "Please specify one likelihood for each view"
  assert set(model_opts['likelihood']).issubset(set(["gaussian","bernoulli","poisson","warp"]))

  # Define whether to learn the feature-wise means
  model_opts["learnIntercept"] = args.learnIntercept

  # Define for which factors and views should we learn 'theta', the sparsity of the factor
  if type(args.learnTheta) == int:
    model_opts['learnTheta'] = [s.ones(K)*args.learnTheta for m in range(M)]
  elif type(args.learnTheta) == list:
    assert len(args.learnTheta) == M, "--learnTheta has to be a binary vector with length number of views"
    model_opts['learnTheta'] = [ args.learnTheta[m]*s.ones(K) for m in range(M) ]
  else:
     print("--learnTheta has to be either 1 or 0 or a binary vector with length number of views")
     exit()

  # Define schedule of updates
  if args.schedule is None:
    model_opts['schedule'] = ( "Y", "SW", "Z", "AlphaW", "Theta", "Tau" )
  else:
    model_opts['schedule'] = args.schedule


  ####################################
  ## Define priors (P distribution) ##
  ####################################

  # Latent Variables
  model_opts["priorZ"] = { 'mean':s.zeros((N,K)) }
  model_opts["priorZ"]['var'] = s.ones((K,))*1.

  # Weights
  model_opts["priorSW"] = { 'Theta':[s.nan]*M, 'mean_S0':[s.nan]*M, 'var_S0':[s.nan]*M, 'mean_S1':[s.nan]*M, 'var_S1':[s.nan]*M } # Not required
  model_opts["priorAlphaW"] = { 'a':[s.ones(K)*1e-14]*M, 'b':[s.ones(K)*1e-14]*M }

  # Theta
  model_opts["priorTheta"] = { 'a':[s.ones(K,) for m in range(M)], 'b':[s.ones(K,) for m in range(M)] }
  for m in range(M):
    for k in range(K):
      if model_opts['learnTheta'][m][k]==0:
        model_opts["priorTheta"]["a"][m][k] = s.nan
        model_opts["priorTheta"]["b"][m][k] = s.nan

  # Tau
  model_opts["priorTau"] = { 'a':[s.ones(D[m])*1e-14 for m in range(M)], 'b':[s.ones(D[m])*1e-14 for m in range(M)] }


  ##############################################
  ## Define initialisations of Q distribution ##
  ##############################################

  # Latent variables
  model_opts["initZ"] = { 'mean':"random", 'var':s.ones((K,)), 'E':None, 'E2':None }

  # Tau
  model_opts["initTau"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[s.ones(D[m])*100 for m in range(M)] }

  # ARD of weights
  model_opts["initAlphaW"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[s.ones(K)*1. for m in range(M)] }

  # Theta
  model_opts["initTheta"] = { 'a':[s.ones(K,) for m in range(M)], 'b':[s.ones(K,) for m in range(M)], 'E':[s.nan*s.zeros((D[m],K)) for m in range(M)] }
  if type(args.initTheta) == float:
    model_opts['initTheta']['E'] = [s.ones((D[m],K))*args.initTheta for m in range(M)]
  elif type(args.initTheta) == list:
    assert len(args.initTheta) == M, "--initTheta has to be a binary vector with length number of views"
    model_opts['initTheta']['E']= [ args.initTheta[m]*s.ones((D[m],K)) for m in range(M) ]
  else:
     print("--learnTheta has to be either 1 or 0 or a binary vector with length number of views")
     exit()

  for m in range(M):
    for k in range(K):
      if model_opts['learnTheta'][m][k]==0.:
        model_opts["initTheta"]["a"][m][k] = s.nan
        model_opts["initTheta"]["b"][m][k] = s.nan

  # Weights
  model_opts["initSW"] = { 
    'Theta':[ model_opts['initTheta']['E'][m] for m in range(M)],
    'mean_S0':[s.zeros((D[m],K)) for m in range(M)],
    'var_S0':[s.nan*s.ones((D[m],K)) for m in range(M)],
    'mean_S1':[s.zeros((D[m],K)) for m in range(M)],
    # 'mean_S1':[stats.norm.rvs(loc=0, scale=1, size=(D[m],K)) for m in range(M)],
    'var_S1':[s.ones((D[m],K)) for m in range(M)],
    'ES':[None]*M, 'EW_S0':[None]*M, 'EW_S1':[None]*M # It will be calculated from the parameters
  }


  ##########################################################
  ## Modify priors and initialisations for the covariates ##
  ##########################################################

  # Covariates are constant vectors and do not have any prior or variational distribution on Z

  if data_opts['covariates'] is not None:
    idx = range(data_opts['covariates'].shape[1])

    # Ignore prior distributions
    model_opts["priorZ"]["var"][idx] = s.nan

    # Ignore variational distribution
    # model_opts["initZ"]["mean"][:,idx] = model_opts["covariates"]
    model_opts["initZ"]["var"][idx] = 0.


  ###########################################################
  ## Modify priors and initialisations for the mean vector ##
  ###########################################################

  # By definition, the weights of the vector associated with the means should not be sparse, therefore we remove
  # the spike and slab prior by not learning theta and initialisating it to one

  if model_opts["learnIntercept"]:
    for m in range(M):
      
      # Weights
      if args.likelihoods[m]=="gaussian":
        model_opts["initSW"]["mean_S1"][m][:,0] = data[m].mean(axis=0)
        model_opts["initSW"]["var_S1"][m][:,0] = 1e-5

      # Theta
      model_opts['learnTheta'][m][0] = 0.
      model_opts["initSW"]["Theta"][m][:,0] = 1.
      model_opts["priorTheta"]['a'][m][0] = s.nan
      model_opts["priorTheta"]['b'][m][0] = s.nan
      model_opts["initTheta"]["a"][m][0] = s.nan
      model_opts["initTheta"]["b"][m][0] = s.nan
      model_opts["initTheta"]["E"][m][:,0] = 1.


  #################################
  ## Define the training options ##
  #################################

  train_opts = {}

  # Maximum number of iterations
  train_opts['maxiter'] = args.iter

  # Lower bound computation frequency
  train_opts['elbofreq'] = args.elbofreq

  # Verbosity
  train_opts['verbose'] = args.verbose

  # Criteria to drop latent variables while training
  train_opts['drop'] = { "by_norm":None, "by_pvar":None, "by_cor":None, "by_r2":args.dropR2 }
  train_opts['startdrop'] = args.startDrop
  train_opts['freqdrop'] = args.freqDrop

  # Tolerance level for convergence
  train_opts['tolerance'] = args.tolerance

  # Do no stop even when convergence criteria is met
  train_opts['forceiter'] = args.nostop

  # Iteration to activate spike and slab sparsity
  train_opts['startSparsity'] = args.startSparsity

  # Number of trials
  train_opts['trials'] = args.ntrials


  #####################
  ## Train the model ##
  #####################

  # Keep the trial with the highest lower bound?
  keep_best_run = False

  # Go!
  # runSingleTrial(data, data_opts, model_opts, train_opts, seed=None)
  runMultipleTrials(data, data_opts, model_opts, train_opts, keep_best_run, args.seed)
