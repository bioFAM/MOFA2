import argparse
import pandas as pd
import scipy as s
from time import time

from build_model import build_model
from biofam.build_model.train_model import train_model
from biofam.build_model.utils import *

def entry_point():

    # Read arguments
    p = argparse.ArgumentParser( description='Basic run script for BioFAM' )

    # I/O
    p.add_argument( '--inFiles',           type=str, nargs='+', required=True,                  help='Input data files (including extension)' )
    p.add_argument( '--outFile',           type=str, required=True,                             help='Output data file (hdf5 format)' )
    p.add_argument( '--delimiter',         type=str, default=" ",                               help='Delimiter for input files' )
    p.add_argument( '--covariatesFile',    type=str, default=None,                              help='Input data file for covariates' )
    p.add_argument( '--header_cols',       action='store_true',                                 help='Do the input files contain column names?' )
    p.add_argument( '--header_rows',       action='store_true',                                 help='Do the input files contain row names?' )

    # Data options
    p.add_argument( '--center_features',   action="store_true",                                 help='Center the features to zero-mean?' )
    p.add_argument( '--scale_features',    action="store_true",                                 help='Scale the features to unit variance?' )
    p.add_argument( '--scale_views',       action="store_true",                                 help='Scale the views to unit variance?' )
    p.add_argument( '--scale_covariates',  type=int, nargs='+', default=0,                      help='Scale covariates?' )

    # Model options
    p.add_argument( '--factors',           type=int, default=10,                                help='Initial number of latent variables')
    p.add_argument( '--likelihoods',       type=str, nargs='+', required=True,                  help='Likelihood per view, current options are bernoulli, gaussian, poisson')
    p.add_argument( '--views',             type=str, nargs='+', required=True,                  help='View names')
    p.add_argument( '--learnIntercept',    action='store_true',                                 help='Learn the feature-wise mean?' )

    # Training options
    p.add_argument( '--elbofreq',          type=int, default=1,                                 help='Frequency of computation of ELBO' )
    p.add_argument( '--iter',              type=int, default=5000,                              help='Maximum number of iterations' )
    p.add_argument( '--startSparsity',     type=int, default=100,                               help='Iteration to activate the spike-and-slab')
    p.add_argument( '--tolerance',         type=float, default=0.01 ,                           help='Tolerance for convergence (based on the change in ELBO)')
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



    ###############
    ## Load data ##
    ###############

    data = loadData(data_opts)

    # Calculate dimensionalities
    N = data[0].shape[0]
    D = [data[m].shape[1] for m in range(M)]

    # Extract sample and features names
    data_opts['sample_names'] = data[0].index.tolist()
    data_opts['feature_names'] = [ data[m].columns.values.tolist() for m in range(len(data)) ]


    #####################
    ## Load covariates ##
    #####################

    model_opts = {}
    model_opts['learnIntercept'] = args.learnIntercept

    if args.covariatesFile is not None:
        model_opts['covariates'] = pd.read_csv(args.covariatesFile, delimiter=" ", header=None).as_matrix()
        print("Loaded covariates from " + args.covariatesFile + "with shape " + str(model_opts['covariates'].shape) + "...")
        model_opts['scale_covariates'] = 1
        args.factors += model_opts['covariates'].shape[1]
    else:
        model_opts['scale_covariates'] = False
        model_opts['covariates'] = None


    ##############################
    ## Define the model options ##
    ##############################

    # Define initial number of latent factors
    model_opts['K'] = args.factors

    # Define likelihoods
    model_opts['likelihood'] = args.likelihoods
    assert M==len(model_opts['likelihood']), "Please specify one likelihood for each view"
    assert set(model_opts['likelihood']).issubset(set(["gaussian","bernoulli","poisson"]))

    #################################
    ## Define the training options ##
    #################################

    train_opts = {}
    train_opts['maxiter'] = args.iter                     # Maximum number of iterations
    train_opts['elbofreq'] = args.elbofreq                # Lower bound computation frequency
    train_opts['verbose'] = args.verbose                  # Verbosity
    train_opts['drop'] = { "by_r2":args.dropR2 }          # Minimum fraction of variance explained to drop latent variables while training
    train_opts['startdrop'] = args.startDrop              # Initial iteration to start dropping factors
    train_opts['freqdrop'] = args.freqDrop                # Frequency of dropping factors
    train_opts['tolerance'] = args.tolerance              # Tolerance level for convergence
    train_opts['forceiter'] = args.nostop                 # Do no stop even when convergence criteria is met
    train_opts['startSparsity'] = args.startSparsity      # Iteration to activate spike and slab sparsity
    if args.seed is None or args.seed==0:                 # Seed for the random number generator
        train_opts['seed'] = int(round(time()*1000)%1e6)
        s.random.seed(train_opts['seed'])
    train_opts['schedule'] = ( "Y", "SW", "Z", "AlphaW", "Theta", "Tau" ) # Define schedule of updates


    #####################
    ## Build the model ##
    #####################

    model = build_model(data, model_opts)

    #####################
    ## Train the model ##
    #####################

    train_model(model, train_opts)

    ################
    ## Save model ##
    ################

    print("Saving model in %s...\n" % data_opts['outfile'])
    train_opts['schedule'] = '_'.join(train_opts['schedule'])
    saveTrainedModel(model=model, outfile=data_opts['outfile'], train_opts=train_opts, model_opts=model_opts,
        view_names=data_opts['view_names'], sample_names=data_opts['sample_names'], feature_names=data_opts['feature_names'])


if __name__ == '__main__':
  entry_point()
