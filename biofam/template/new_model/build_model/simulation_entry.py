import argparse

from build_model import build_model
from biofam.build_model.simulate_model import simulate_model
from biofam.build_model.utils import *

#TODO : make it possible for views with different dimensions
#TODO : enable covariatesFiles

def entry_point():

    # Read arguments
    p = argparse.ArgumentParser( description='Basic run script for BioFAM' )

    # I/O
    p.add_argument( '--outFile',           type=str, required=True,                             help='Output data file (hdf5 format)' )
    p.add_argument( '--outDir',            type=str, required=True,                             help='Output dir for text files' )
    p.add_argument( '--covariatesFiles',   type=str, nargs='+', default=None,                   help='Input data file for covariates')

    # Data options
    p.add_argument( '--scale_covariates',  type=int, nargs='+', default=0,                      help='Scale covariates?' )

    # simulation options
    p.add_argument( '--factors',           type=int, default=10,                                help='Initial number of latent variables')
    p.add_argument( '--spatialFact',       type=float, default=0.,                              help='Initial percentage of non-spatial latent variables')
    p.add_argument( '--transpose',         type=int, default=0,                                 help='Use the transpose MOFA (a view is a population of cells, not an omic) ? ')
    p.add_argument( '--sample_X',          type=int, default=0,                                 help='Sample the positions of the samples to test covariance prior structure per factor')

    p.add_argument( '--likelihoods',       type=str, nargs='+', required=True,                  help='Likelihood per view, current options are bernoulli, gaussian, poisson')
    p.add_argument( '--views',             type=str, nargs='+', required=True,                  help='View names')
    p.add_argument( '--noise',             type=float, default=10.,                              help='noise level for simulations')
    p.add_argument('--sparsity',           type=float, default=.1,                              help='sparsity level for simulations')
    p.add_argument( '--N',                 type=int, default=500,                               help='number of samples to simulate for')
    p.add_argument( '--D',                 type=int, nargs='+', required=True,                  help='number of features per view to simulate for')
    p.add_argument( '--M',                 type=int, default=3,                               help='number of views to simulate from ')

    args = p.parse_args()

    # Calculate dimensionalities
    N = args.N
    M = args.M
    D = args.D

    #####################
    ## Load covariates ##
    #####################
    model_opts = {}

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
    ## Define the model options ##
    ##############################

    # Define initial number of latent factors
    model_opts['K'] = args.factors
    model_opts['N'] = N
    model_opts['M'] = M
    model_opts['D'] = D

    # Tell if the multi-view MOFA model is used transposed (1 : Yes, 0 : No)
    model_opts['transpose'] = args.transpose

    # Define likelihoods
    model_opts['likelihood'] = args.likelihoods
    if (len(model_opts['likelihood']) == 1) and (M!=1) :  # TODO check that
        model_opts['likelihood'] = [model_opts['likelihood']] * M
    assert M==len(model_opts['likelihood']), "Please specify one likelihood for each view"
    model_opts['learnIntercept'] = False # TODO could use that to simulate intercept
    model_opts['outDir'] = args.outDir
    model_opts['noise'] = args.noise
    model_opts['sparsity'] = args.sparsity
    model_opts['spatialFact'] = args.spatialFact

    #Use sampled positions of the samples
    model_opts['sample_X'] = args.sample_X

    #####################
    ## Build the model ##
    #####################

    model = build_model(model_opts, data = None)

    #####################
    ## Simulate the model ##
    #####################

    simulate_model(model)

    ################
    ## Save model ##
    ################

    print("Saving model in %s...\n" % args.outFile)
    print(model_opts["outDir"]) #, args.outFile)

    saveSimulatedModel(model=model, outfile=args.outFile, view_names=args.views, train_opts=None, model_opts=model_opts)



if __name__ == '__main__':
  entry_point()
