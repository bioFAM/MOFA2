import argparse
import pandas as pd
import scipy as s
from time import time

from build_model import build_model
from biofam.build_model.train_model import train_model
from biofam.build_model.utils import *

#TODO : enable covariatesFiles

def entry_point():

    # Read arguments
    p = argparse.ArgumentParser( description='Basic run script for BioFAM' )

    # I/O
    p.add_argument( '--inFiles',           type=str, nargs='+', required=True,                  help='Input data files (including extension)' )
    p.add_argument( '--outFile',           type=str, required=True,                             help='Output data file (hdf5 format)' )
    p.add_argument( '--delimiter',         type=str, default=" ",                               help='Delimiter for input files' )
    p.add_argument( '--header_cols',       action='store_true',                                 help='Do the input files contain column names?' )
    p.add_argument( '--header_rows',       action='store_true',                                 help='Do the input files contain row names?' )
    p.add_argument( '--covariatesFiles',   type=str, nargs='+', default=None,                              help='Input data file for covariates')
    p.add_argument( '--X_Files',           type=str, nargs='+', default=None,                              help='Use positions of samples for covariance prior structure per factor')
    p.add_argument( '--sigmaClusterFiles', type=str, nargs='+', default=None,                              help='Use clusters assigned to samples for a block covariance prior structure per factor')
    p.add_argument( '--permute_samples',   type=int, default=0,                                 help='Permute samples positions in the data')

    # Data options
    p.add_argument( '--center_features',    action="store_true",                          help='Center the features to zero-mean?' )
    p.add_argument( '--scale_features',     action="store_true",                          help='Scale the features to unit variance?' )
    p.add_argument( '--scale_views',        action="store_true",                          help='Scale the views to unit variance?' )
    p.add_argument( '--scale_covariates',   type=int, nargs='+', default=0,               help='Scale covariates?' )
    p.add_argument( '--shared_features',    action="store_true", default=False,           help='Features, not samples are shared between views?' )

    # Model options
    p.add_argument('--transpose',           action='store_true', help='Noise and sparsity across the common dimension?')
    p.add_argument('--transpose_noise',     action='store_true', help='Noise in the common dimension?')
    p.add_argument('--transpose_sparsity',  action='store_true', help='Sparsity across the common dimension?')
    p.add_argument( '--factors',           type=int, default=10,                                help='Initial number of latent variables')
    p.add_argument( '--likelihoods',       type=str, nargs='+', required=True,                  help='Likelihood per view, current options are bernoulli, gaussian, poisson')
    p.add_argument( '--views',             type=str, nargs='+', required=True,                  help='View names')
    p.add_argument( '--learnIntercept',    action='store_true',                                 help='Learn the feature-wise mean?' )
    p.add_argument('--ARD_per_view',  action='store_true',default=True, help='ARD prior per view ? (relevant option if transpose_sparsity=1, X_Files=None and sample_X=None)')
    p.add_argument( '--sample_X',         type=int, default=0,                                  help='Sample the positions of the samples to test covariance prior structure per factor' )

    # Training options
    p.add_argument( '--elbofreq',           type=int, default=1,                          help='Frequency of computation of ELBO' )
    p.add_argument( '--iter',               type=int, default=5000,                       help='Maximum number of iterations' )
    p.add_argument( '--startSparsity',      type=int, default=100,                        help='Iteration to activate the spike-and-slab')
    p.add_argument( '--tolerance',          type=float, default=0.01 ,                    help='Tolerance for convergence (based on the change in ELBO)')
    p.add_argument( '--startDrop',          type=int, default=1 ,                         help='First iteration to start dropping factors')
    p.add_argument( '--freqDrop',           type=int, default=1 ,                         help='Frequency for dropping factors')
    p.add_argument( '--dropR2',             type=float, default=None ,                    help='Threshold to drop latent variables based on coefficient of determination' )
    p.add_argument( '--nostop',             action='store_true',                          help='Do not stop when convergence criterion is met' )
    p.add_argument( '--verbose',            action='store_true',                          help='Use more detailed log messages?')
    p.add_argument( '--seed',               type=int, default=0 ,                         help='Random seed' )

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

    # Data options: if features are in rows or in columns
    data_opts['shared_features'] = args.shared_features


    ###############
    ## Load data ##
    ###############

    data = loadData(data_opts)

    # Calculate dimensionalities
    N = data[0].shape[0]
    D = [data[m].shape[1] for m in range(M)]

    # Extract sample and features names
    if (data_opts['shared_features']):
        data_opts['feature_names']  = data[0].index.tolist()
        data_opts['sample_names'] = [ data[m].columns.values.tolist() for m in range(len(data)) ]
    else:
        data_opts['sample_names']  = data[0].index.tolist()
        data_opts['feature_names'] = [ data[m].columns.values.tolist() for m in range(len(data)) ]


    #####################
    ## Load covariates ##
    #####################

    model_opts = {}

    if args.covariatesFiles is not None:
        if args.transpose_sparsity:
            assert M == len(args.covariatesFiles), "Length of view names and covariates input files does not match"
        else:
            assert 1 == len(args.covariatesFiles), "Length of view names and covariates input files does not match"

    model_opts['learnIntercept'] = args.learnIntercept

    # TODO : create a function loadDataCovariates in utils.py as loadData or loadDataX, and call it below
    if args.covariatesFiles is not None:
        model_opts['covariatesFiles'] = pd.read_csv(args.covariatesFiles, delimiter=" ", header=None).as_matrix()
        print("Loaded covariates from " + args.covariatesFiles + "with shape " + str(model_opts['covariatesFiles'].shape) + "...")
        model_opts['scale_covariates'] = 1
        args.factors += model_opts['covariatesFiles'].shape[1]
    else:
        model_opts['scale_covariates'] = False
        model_opts['covariatesFiles'] = None

    ##############################
    ## Load X and cluster files ##
    ##############################

    data_opts['X_Files'] = args.X_Files
    data_opts['sigmaClusterFiles'] = args.sigmaClusterFiles
    data_opts['permute_samples'] = args.permute_samples

    dataX, dataClust = loadDataX(data_opts, transpose = args.transpose_sparsity)

    ##############################
    ## Define the model options ##
    ##############################

    # Choose between MOFA (view = omic) and transposed MOFA (view = group of samples or cell population)
    model_opts['transpose'] = args.transpose
    model_opts['transpose_noise'] = args.transpose_noise
    model_opts['transpose_sparsity'] = args.transpose_sparsity
    # To keep reverse-compatibility
    if model_opts['transpose']:
        model_opts['transpose_noise'] = True
        model_opts['transpose_sparsity'] = True
    if model_opts['transpose']: print("Using features as a shared dimension...")

    # Choose or not to simulate the positions of samples
    model_opts['sample_X'] = args.sample_X


    # Define initial number of latent factors
    model_opts['K'] = args.factors

    # Define likelihoods
    model_opts['likelihood'] = args.likelihoods
    assert M == len(model_opts['likelihood']), "Please specify one likelihood for each view"
    assert set(model_opts['likelihood']).issubset(set(["gaussian", "bernoulli", "poisson"]))

    model_opts["ARD_per_view"] = args.ARD_per_view

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
    else:
        train_opts['seed'] = args.seed
    s.random.seed(train_opts['seed'])

    # Define schedule of updates
    # Think to its importance ?
    if model_opts['transpose_sparsity']:
        if (dataX is not None) or (model_opts['sample_X']):
            train_opts['schedule'] = ( "Y", "SZ", "W", "SigmaAlphaW", "AlphaZ", "ThetaZ", "Tau" )
        else:
            if model_opts["ARD_per_view"]:
                train_opts['schedule'] = ( "Y", "SZ", "W", "AlphaW", "ThetaZ", "Tau")
            else:
                train_opts['schedule'] = ( "Y", "SZ", "W", "AlphaZ", "ThetaZ", "Tau")
    else:
        if (dataX is not None) or (model_opts['sample_X']):
            train_opts['schedule'] = ( "Y", "SW", "Z", "AlphaW", "SigmaZ", "ThetaW", "Tau" )
        else:
            train_opts['schedule'] = ( "Y", "SW", "Z", "AlphaW", "ThetaW", "Tau" )

    print("transposed sparsity : ", model_opts['transpose_sparsity'])
    print("transposed noise : ", model_opts['transpose_noise'])
    print("schedule for train : ", train_opts["schedule"])

    #####################
    ## Build the model ##
    #####################

    model = build_model(model_opts, data=data, dataX=dataX, dataClust=dataClust)

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
        view_names=data_opts['view_names'], sample_names=data_opts['sample_names'], feature_names=data_opts['feature_names'], shared_features=data_opts['shared_features'])


if __name__ == '__main__':
  entry_point()
