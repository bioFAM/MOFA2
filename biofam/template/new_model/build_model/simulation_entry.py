import argparse
import pandas as pd
import scipy as s
from time import time

from build_model import build_model
from biofam.build_model.simulate_model import simulate_model
from biofam.build_model.utils import *

# TODO make it possible to input multidimensional D

def entry_point():

    # Read arguments
    p = argparse.ArgumentParser( description='Basic run script for BioFAM' )

    # I/O
    p.add_argument( '--outFile',           type=str, required=True,                             help='Output data file (hdf5 format)' )
    p.add_argument( '--outDir',            type=str, required=True,                             help='Output dir for text files' )
    p.add_argument( '--covariatesFile',    type=str, default=None,                              help='Input data file for covariates' )

    # Data options
    p.add_argument( '--scale_covariates',  type=int, nargs='+', default=0,                      help='Scale covariates?' )

    # simulation options
    p.add_argument( '--factors',           type=int, default=10,                                help='Initial number of latent variables')
    p.add_argument( '--likelihoods',       type=str, nargs='+', required=True,                  help='Likelihood per view, current options are bernoulli, gaussian, poisson')
    # p.add_argument( '--SpatialFact',       type=float, default=0.2,                             help='fraction of spatial facts for simulations')  # TODO add that to spatialFA
    p.add_argument( '--noise',             type=float, default=10.,                              help='noise level for simulations')
    p.add_argument( '--N',                 type=int, default=500,                               help='number of samples to simulate for')
    p.add_argument( '--D',                 type=int, default=1000,                               help='number of features per view to simulate for')
    p.add_argument( '--M',                 type=int, default=3,                               help='number of views to simulate from ')
    p.add_argument( '--seed',              type=int, default=0 ,                                help='Random seed' )


    args = p.parse_args()

    # Calculate dimensionalities
    N = args.N
    M = args.M
    D = [args.D for m in range(M)]

    #####################
    ## Load covariates ##
    #####################
    model_opts = {}

    data = None

    # could keep covariates for simulations
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
    model_opts['N'] = N
    model_opts['M'] = M
    model_opts['D'] = D

    # Define likelihoods
    model_opts['likelihood'] = args.likelihoods
    if len(model_opts['likelihood']) == 1:  # TODO check that
        model_opts['likelihood'] = [model_opts['likelihood']] * M
    assert M==len(model_opts['likelihood']), "Please specify one likelihood for each view"
    assert set(model_opts['likelihood']).issubset(set(["gaussian","bernoulli","poisson"]))
    model_opts['learnIntercept'] = False # TODO could use that to simulate intercept
    model_opts['outDir'] = args.outDir
    #####################
    ## Build the model ##
    #####################

    model = build_model(data, model_opts)

    #####################
    ## Simulate the model ##
    #####################

    model.simulate()

    ################
    ## Save model ##
    ################

    print("Saving model in %s...\n" % args.outFile)
    saveSimulatedModel(model=model, outfile=args.outFile, train_opts=None, model_opts=model_opts)


if __name__ == '__main__':
  entry_point()
