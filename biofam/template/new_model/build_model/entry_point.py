import argparse
import pandas as pd
import scipy as s
from time import time

from build_model import build_model
from biofam.build_model.train_model import train_model
from biofam.build_model.utils import *

#TODO : enable covariatesFiles
# TODO change the names coming from spatialFA
def entry_point():

    # Read arguments
    p = argparse.ArgumentParser( description='Basic run script for BioFAM' )

    # I/O: Main files
    p.add_argument('--input-files',           type=str, nargs='+', required=True,  help='Input data files (including extension)')
    p.add_argument('--output-file',           type=str, required=True,             help='Output data file (hdf5 format)')

    # I/O: Files descriptions
    p.add_argument('--views',                 type=str, nargs='+', required=False, help='View names')
    p.add_argument('--groups',                type=str, nargs='+', required=False, help='Group names')
    p.add_argument('--treat-as-views',        action='store_true',                 help='Treat all the input files as views')
    p.add_argument('--treat-as-groups',       action='store_true',                 help='Treat all the input files as groups')

    # Data files format options
    p.add_argument('--delimiter',             type=str, default=" ",               help='Delimiter for input files')
    p.add_argument('--header-cols',           action='store_true',                 help='Treat first row as of input files as column names')
    p.add_argument('--header-rows',           action='store_true',                 help='Treat first column of input files as row names')
    p.add_argument('--samples-in-rows',       action="store_true",                 help='Samples are in rows of the input files (default)')
    p.add_argument('--features-in-rows',      action="store_true", default=False,  help='Features (e.g. genes) are in rows of the input files')

    # Data normalisation options
    p.add_argument('--center-features',       action="store_true",                 help='Center the features to zero-mean')
    p.add_argument('--scale-features',        action="store_true",                 help='Scale the features to unit variance')
    p.add_argument('--scale-views',           action="store_true",                 help='Scale the views to unit variance')
    # p.add_argument('--scale-groups',          action="store_true",                 help='Scale sample groups to unit variance')  # TODO: implement per-group scaling
    p.add_argument('--scale-covariates',      type=int, nargs='+', default=0,      help='Scale covariates' )

    # I/O: Additional files
    p.add_argument('--sample-groups-file',    type=str, default=None,              help='If samples contain groups, file containing the labels of the samples')
    p.add_argument('--covariates-files',      type=str, nargs='+', default=None,   help='Input data file for covariates')
    p.add_argument('--x-files',               type=str, nargs='+', default=None,   help='Use positions of samples for covariance prior structure per factor')
    p.add_argument('--sigma-cluster-files',   type=str, nargs='+', default=None,   help='Use clusters assigned to samples for a block covariance prior structure per factor')

    # Model options
    p.add_argument('--feature-wise-noise',    action='store_true', default=True,   help='Noise parameter per feature (e.g. gene)' )
    p.add_argument('--sample-wise-noise',     action='store_true', default=False,  help='Noise parameter per sample')
    p.add_argument('--feature-wise-sparsity', action='store_true', default=True,   help='Sparsity across features (e.g. genes)')
    p.add_argument('--sample-wise-sparsity',  action='store_true', default=False,  help='Sparsity across samples')
    p.add_argument('--factors',               type=int, default=10,                help='Initial number of latent variables')
    p.add_argument('--likelihoods',           type=str, nargs='+', required=True,  help='Likelihood per view, current options are bernoulli, gaussian, poisson')
    p.add_argument('--learn-intercept',       action='store_true',                 help='Learn the feature-wise mean?')
    p.add_argument('--ard-per-view',          action='store_false',                help='ARD prior per view (relevant option if --sample-wise-sparsity, no x_files, and no sample_x)')  # TODO: add symmetry with ard-per-group?

    # Training options
    p.add_argument('--elbofreq',              type=int, default=1,                 help='Frequency of computation of ELBO')
    p.add_argument('--iter',                  type=int, default=5000,              help='Maximum number of iterations')
    p.add_argument('--start-sparsity',        type=int, default=100,               help='Iteration to activate the spike-and-slab')
    p.add_argument('--tolerance',             type=float, default=0.01 ,           help='Tolerance for convergence (based on the change in ELBO)')
    p.add_argument('--start-drop',            type=int, default=1 ,                help='First iteration to start dropping factors')
    p.add_argument('--freq-drop',             type=int, default=1 ,                help='Frequency for dropping factors')
    p.add_argument('--drop-r2',               type=float, default=None ,           help='Threshold to drop latent variables based on coefficient of determination')
    p.add_argument('--non-stop',              action='store_true',                 help='Do not stop when convergence criterion is met')
    p.add_argument('--verbose',               action='store_true',                 help='Use more detailed log messages')
    p.add_argument('--seed',                  type=int, default=0 ,                help='Random seed' )


    # TODO: put options below into respective categories
    # TODO: adapt options below for having clusters on Z
    
    # Misc
    p.add_argument( '--permute-samples',      action='store_true', default=False,  help='Permute samples positions in the data')
    
    # Simulations
    p.add_argument( '--sample_x',             action='store_true', default=False,  help='Sample the positions of the samples to test covariance prior structure per factor' )

    args = p.parse_args()


    #############################
    ## Define the data options ##
    #############################

    data_opts = {}

    # I/O
    data_opts['input_files'] = args.input_files
    data_opts['output_file'] = args.output_file

    # Views and groups options
    data_opts['view_names']  = args.views
    data_opts['group_names'] = args.groups

    # File format
    data_opts['delimiter'] = args.delimiter
    data_opts['rownames'] = 0 if args.header_rows else None
    data_opts['colnames'] = 0 if args.header_rows else None
    data_opts['features_in_rows'] = args.features_in_rows if args.features_in_rows else not args.samples_in_rows


    #####################
    ## Views & groups ###
    #####################

    if args.sample_groups_file:
        data_opts['sample_groups_file'] = args.sample_groups_file
        data_opts['sample_groups'] = loadDataGroups(data_opts)
    
    if data_opts['group_names'] is None and not args.treat_as_groups:
        # All the files are views
        # Please note that this asymmetry is provided in order to keep reverse-compatibility with MOFA software behaviour
        if data_opts['view_names'] is None:
            data_opts['view_names'] = ["view_".format(i) for i in range(len(data_opts['input_files']))]
        else:
            assert len(data_opts["view_names"]) == len(data_opts['input_files']), "The number of view names and the number of input files do not match."
        if args.sample_groups_file is None:
            # No group annotation provided, then there's only one group of samples
            data_opts['group_names'] = ['group_0']
        else:
            # Group annotation provided in a separate file
            # Group names will be obtained from the annotation file
            print("All files are treated as views since group annotation is provided in {}".format(data_opts['sample_groups_file']))
            assert data_opts['group_names'] is None, "Both group names and a file with sample group names are provided. Please provide group as separate files or define group annotation for samples in a separate file while keeping samples for all the groups in one file per view."
            assert args.treat_as_groups is None, "The option --treat-as-groups is not available since group annotation for samples is provided in {}".format(data_opts['sample_groups_file'])
    elif not args.treat_as_groups:
        # group_names do exist in data_opts
        assert len(data_opts['group_names']) == len(data_opts['input_files']), "The number of group names and the number of input files do not match"
    else:  # data_opts['group_names'] is None
        data_opts['group_names'] = ["group_".format(i) for i in range(len(data_opts["input_files"]))]
        

    
    # data_opts['multi_group'] = not (data_opts['group_names'] is None or data_opts['sample_groups_file'] is None)

    # # Fill missing group info in
    # if data_opts['multi_group']:
    #     if data_opts['group_names'] is not None and data_opts['sample_groups'] is not None:
    #         print("Warning: both group names and sample groups are provided.")
    #     if data_opts['group_names'] is None:
    #         data_opts['group_names'] = list(set(data_opts['sample_groups']))
    #     if data_opts['sample_groups'] is None:
    #         data_opts['sample_groups'] = 

    # data_opts['sample_groups_file'] = args.sample_groups



    #####################
    ## Data processing ##
    #####################

    # Check that the likelihood is provided for every input file
    assert len(args.input_files) == len(args.likelihoods), "Please specify one likelihood for each input file"

    # Data processing: center features
    if args.center_features:
        data_opts['center_features'] = [ True if l=="gaussian" else False for l in args.likelihoods ]
    else:
        if not args.learn_intercept: print("\nWarning... you are not centering the data and not learning the mean...\n")
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
    if not args.samples_in_rows and not args.features_in_rows: args.samples_in_rows = True
    assert args.features_in_rows != args.samples_in_rows, "Please choose if features or samples are in rows"
    data_opts['features_in_rows'] = args.features_in_rows


    ###############
    ## Load data ##
    ###############

    data, sample_groups = loadData(data_opts)

    if 'sample_groups' not in data_opts:
        data_opts['sample_groups'] = sample_groups

    # Calculate dimensionalities
    M = len(set(data_opts['view_names']))
    P = len(set(data_opts['group_names']))
    N = data[0].shape[0]
    D = [data[m].shape[1] for m in range(M)]

    # Extract sample and features names
    # if (data_opts['shared_features']):
    #     data_opts['feature_names']  = data[0].index.tolist()
    #     data_opts['sample_names'] = [ data[m].columns.values.tolist() for m in range(len(data)) ]
    # else:
    
    data_opts['sample_names']  = data[0].index.tolist()
    data_opts['feature_names'] = [ data[m].columns.values.tolist() for m in range(len(data)) ]


    #####################
    ## Load covariates ##
    #####################

    model_opts = {}

    if args.covariates_files is not None:
        if args.transpose_sparsity:
            assert M == len(args.covariates_files), "Length of view names and covariates input files does not match"
        else:
            assert 1 == len(args.covariates_files), "Length of view names and covariates input files does not match"

    model_opts['learn_intercept'] = args.learn_intercept

    # TODO : create a function loadDataCovariates in utils.py as loadData or loadDataX, and call it below
    if args.covariates_files is not None:
        model_opts['covariates_files'] = pd.read_csv(args.covariatesFiles, delimiter=" ", header=None).as_matrix()
        print("Loaded covariates from " + args.covariates_files + "with shape " + str(model_opts['covariates_files'].shape) + "...")
        model_opts['scale_covariates'] = 1
        args.factors += model_opts['covariates_files'].shape[1]
    else:
        model_opts['scale_covariates'] = False
        model_opts['covariates_files'] = None

    ##############################
    ## Load X and cluster files ##
    ##############################

    data_opts['x_files'] = args.x_files
    data_opts['sigma_cluster_files'] = args.sigma_cluster_files
    data_opts['permute_samples'] = args.permute_samples
    # TODO: load file here and add to data_opts
    # same for dataX etc
    
    dataX, dataClust = loadDataX(data_opts)

    ##############################
    ## Define the model options ##
    ##############################

    # Choose between MOFA (view = omic) and transposed MOFA (view = group of samples or cell population)
    # model_opts['transpose'] = args.transpose
    # model_opts['transpose_noise'] = args.transpose_noise
    # model_opts['transpose_sparsity'] = args.transpose_sparsity
    # # To keep reverse-compatibility
    # if model_opts['transpose']:
    #     model_opts['transpose_noise'] = True
    #     model_opts['transpose_sparsity'] = True
    # if model_opts['transpose']: print("Using features as a shared dimension...")
    
    # Sparsity settings
    model_opts['feature_wise_sparsity'] = args.feature_wise_sparsity
    model_opts['sample_wise_sparsity'] = args.sample_wise_sparsity
    
    # Noise settings
    model_opts['feature_wise_noise'] = args.feature_wise_noise
    model_opts['sample_wise_noise'] = args.sample_wise_noise

    # Choose or not to simulate the positions of samples
    model_opts['sample_x'] = args.sample_x


    # Define initial number of latent factors
    model_opts['K'] = args.factors

    # Define likelihoods
    model_opts['likelihood'] = args.likelihoods
    assert set(model_opts['likelihood']).issubset(set(["gaussian", "bernoulli", "poisson"]))

    model_opts["ard_per_view"] = args.ard_per_view

    #################################
    ## Define the training options ##
    #################################

    train_opts = {}
    train_opts['maxiter'] = args.iter                            # Maximum number of iterations
    train_opts['elbofreq'] = args.elbofreq                       # Lower bound computation frequency
    train_opts['verbose'] = args.verbose if args.verbose else 2  # Verbosity
    train_opts['drop'] = { "by_r2":args.drop_r2 }                # Minimum fraction of variance explained to drop latent variables while training
    train_opts['start_drop'] = args.start_drop                   # Initial iteration to start dropping factors
    train_opts['freq_drop'] = args.freq_drop                     # Frequency of dropping factors
    train_opts['tolerance'] = args.tolerance                     # Tolerance level for convergence
    train_opts['forceiter'] = args.non_stop                      # Do no stop even when convergence criteria is met
    train_opts['start_sparsity'] = args.start_sparsity           # Iteration to activate spike and slab sparsity
    if args.seed is None or args.seed==0:                        # Seed for the random number generator
        train_opts['seed'] = int(round(time()*1000)%1e6)
    else:
        train_opts['seed'] = args.seed
    s.random.seed(train_opts['seed'])

    # Define schedule of updates
    # Think to its importance ?
    # if model_opts['transpose_sparsity']:
        # if (dataX is not None) or (model_opts['sample_x']):
        #     train_opts['schedule'] = ( "Y", "SZ", "W", "SigmaAlphaW", "AlphaZ", "ThetaZ", "Tau" )
        # else:
        #     if model_opts["ard_per_view"]:
        #         train_opts['schedule'] = ( "Y", "SZ", "W", "AlphaW", "AlphaZ", "ThetaZ", "Tau")
        #     else:
        #         train_opts['schedule'] = ( "Y", "SZ", "W", "AlphaZ", "ThetaZ", "Tau")
    # else:
    # TODO: reflect noise, sparsity, etc. in the schedule
    if (dataX is not None) or (model_opts['sample_x']):
        train_opts['schedule'] = ( "Y", "SW", "Z", "AlphaW", "SigmaZ", "ThetaW", "Tau" )
    elif model_opts['sample_wise_sparsity']:
        train_opts['schedule'] = ( "Y", "SW", "Z", "AlphaZ", "AlphaW", "ThetaW", "Tau" )
    else:
        train_opts['schedule'] = ( "Y", "SW", "Z", "AlphaW", "ThetaW", "Tau" )
    
    # print("transposed sparsity : ", model_opts['transpose_sparsity'])
    # print("transposed noise : ", model_opts['transpose_noise'])
    print("schedule for train : ", train_opts["schedule"])

    #####################
    ## Build the model ##
    #####################
    
    model = build_model(model_opts, data=data, dataX=dataX, dataClust=dataClust, data_groups=sample_groups)

    #####################
    ## Train the model ##
    #####################

    train_model(model, train_opts)

    ################
    ## Save model ##
    ################

    print("Saving model in %s...\n" % data_opts['output_file'])
    train_opts['schedule'] = '_'.join(train_opts['schedule'])
    saveTrainedModel(model=model, outfile=data_opts['output_file'], train_opts=train_opts, model_opts=model_opts,
        view_names=data_opts['view_names'], group_names=data_opts['group_names'], sample_groups=data_opts['sample_groups'], sample_names=data_opts['sample_names'], feature_names=data_opts['feature_names'])


if __name__ == '__main__':
  entry_point()
