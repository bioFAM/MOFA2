import argparse
import pandas as pd
import scipy as s
import sys
from time import sleep

from build_model import *
from biofam.build_model.train_model import train_model

"""
TO-DO:
- ADD DETAILED EXPLANATION OF ARGUMENTS
- ADd sanity checks
- Add print messages
- Within each method, check that the pipeline order is met

Pipeline
(1) Parse data options
(2) Parse train options or parse model options
(3) Parse data processing options
(4) Load the data or define priors or define variational
(5) Train

CHANGE PARSE TO SET
SET PRIORS AND VARIATIONAL DIST IN BUILD_MODEL, ONLY BIOFAM BRANCH
"""

class entry_point():
  def __init__(self):
    self.print_banner()
    self.dimensionalities = {}

  def print_banner(self):
    """ Method to print the MOFA banner """

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
    # sleep(2)

  def set_data_options(self,
    inFiles, outFile, views, groups,
    delimiter=" ", header_cols=False, header_rows=False
    ):
    """ Parse I/O data options """

    # TO-DO: sanity checks
    # - Check that input file exists
    # - Check that output directory exists: warning if not
    # TO-DO: Verbosity, print messages

    self.data_opts = {}

    # I/O
    if type(inFiles) is not list:
      inFiles = [inFiles]
    self.data_opts['input_files'] = inFiles
    self.data_opts['output_file'] = outFile
    self.data_opts['delimiter'] = delimiter

    # View names
    if type(views) is not list:
      views = [views]
    if type(groups) is not list:
      views = [groups]
    self.data_opts['view_names'] = views
    self.data_opts['group_names'] = groups

    # Headers
    if header_rows is True:
      self.data_opts['rownames'] = 0
    else:
      self.data_opts['rownames'] = None

    if header_cols is True:
      self.data_opts['colnames'] = 0
    else:
      self.data_opts['colnames'] = None


    self.dimensionalities['M'] = len(set(self.data_opts['view_names']))
    self.dimensionalities['P'] = len(set(self.data_opts['group_names']))

  def set_train_options(self,
    iter=5000, elbofreq=1, ntrials=1, startSparsity=100, tolerance=0.01,
    startDrop=1, freqDrop=1, dropR2=0, nostop=False, verbose=False, seed=None, schedule=None
    ):
    """ Parse training options """

    # TO-DO: verbosity, print more messages

    self.train_opts = {}

    # Maximum number of iterations
    self.train_opts['maxiter'] = int(iter)

    # Lower bound computation frequency
    self.train_opts['elbofreq'] = int(elbofreq)

    # Verbosity
    self.train_opts['verbose'] = verbose

    # Criteria to drop latent variables while training
    self.train_opts['drop'] = { "by_r2":float(dropR2) }
    self.train_opts['start_drop'] = int(startDrop)
    self.train_opts['freq_drop'] = int(freqDrop)
    if (dropR2>0): print("\nDropping factors with minimum threshold of {0}% variance explained".format(dropR2))


    # Tolerance level for convergence
    self.train_opts['tolerance'] = float(tolerance)

    # Do no stop even when convergence criteria is met
    self.train_opts['forceiter'] = nostop

    # Iteration to activate spike and slab sparsity
    self.train_opts['start_sparsity'] = int(startSparsity)

    # Number of trials
    # TODO check
    self.train_opts['trials'] = int(ntrials)

    # Seed
    if seed is None:
      seed = 0
    self.train_opts['seed'] = int(seed)

    # Define schedule of updates
    if schedule is not None:
      print("\nWarning... we recommend using the default training schedule\n")
      self.train_opts['schedule'] = schedule
    else:
      self.train_opts['schedule'] = ( "Y", "W", "Z", "AlphaW", "ThetaW", "Tau" )

  def set_model(self, sl_z=False, sl_w=True, ard_z=False, ard_w=True):
    self.model_def={}

    self.model_def['sl_z'] = sl_z
    self.model_def['sl_w'] = sl_w

    self.model_def['ard_z'] = ard_z
    self.model_def['ard_w'] = ard_w


  def set_model_options(self, factors, likelihoods, schedule=None, learnTheta=True, learn_intercept=False):
    """ Parse model options """

    # TO-DO: SANITY CHECKS AND:
    # - learnTheta should be replaced by sparsity=True

    self.model_opts = {}

    # Define initial number of latent factors
    K = self.dimensionalities["K"] = self.model_opts['factors'] = int(factors)
    M = self.dimensionalities["M"]

    # Define likelihoods
    self.model_opts['likelihoods'] = likelihoods
    if type(self.model_opts['likelihoods']) is not list:
      self.model_opts['likelihoods'] = [self.model_opts['likelihoods']]

    assert len(self.model_opts['likelihoods'])==M, "Please specify one likelihood for each view"
    assert set(self.model_opts['likelihoods']).issubset(set(["gaussian","bernoulli","poisson"]))

    # Define whether to learn the feature-wise means
    self.model_opts["learn_intercept"] = learn_intercept
    if learn_intercept:
      self.model_opts['factors'] += 1
      self.dimensionalities["K"] += 1
      K += 1

    # Define for which factors and views should we learn 'theta', the sparsity of the factor
    if type(learnTheta) is bool:
      self.model_opts['sparsity'] = True
      self.model_opts['learnTheta'] = [s.ones(K) for m in range(M)]
    # elif type(learnTheta) is list:
    #   self.model_opts['sparsity'] = True
    #   assert len(learnTheta)==M, "--learnTheta has to be a binary vector with length number of views"
    #   self.model_opts['learnTheta'] = [ learnTheta[m]*s.ones(K) for m in range(M) ]
    else:
       print("--learnTheta has to be either 1 or 0 or a binary vector with length number of views")
       exit(1)

    # TODO sort that out
    self.data_opts['features_in_rows'] = False

  # TODO merge with data_options ?
  def set_dataprocessing_options(self,
    center_features=False, scale_features=False, scale_views=False,
    maskAtRandom=None, maskNSamples=None, RemoveIncompleteSamples=False
    ):

    """ Parse data processing options """

    # TO-DO: more verbose messages


    # Sanity checks
    M = self.dimensionalities["M"]
    # assert len(self.data_opts['view_names'])==M, "Length of view names and input files does not match"

    # Data processing: center features
    if center_features is True:
      self.data_opts['center_features'] = [ True if l=="gaussian" else False for l in self.model_opts["likelihoods"] ]
    else:
      if not self.model_opts["learn_intercept"]: print("\nWarning... you are not centering the data and not learning the mean...\n")
      self.data_opts['center_features'] = [ False for l in self.model_opts["likelihoods"] ]

    # Data processing: scale views
    if scale_views is True:
      self.data_opts['scale_views'] = [ True if l=="gaussian" else False for l in self.model_opts["likelihoods"] ]
    else:
      self.data_opts['scale_views'] = [ False for l in self.model_opts["likelihoods"] ]

    # Data processing: scale features
    if scale_features:
      assert data_opts['scale_views'] is False, "Scale either entire views or features, not both"
      self.data_opts['scale_features'] = [ True if l=="gaussian" else False for l in self.model_opts["likelihoods"] ]
    else:
      self.data_opts['scale_features'] = [ False for l in self.model_opts["likelihoods"] ]


    # Data processing: mask values
    if maskAtRandom is not None:
      self.data_opts['maskAtRandom'] = data_opts['maskAtRandom']
      assert len(self.data_opts['maskAtRandom'])==M, "Length of MaskAtRandom and input files does not match"
    else:
      self.data_opts['maskAtRandom'] = [0]*M

    if maskNSamples is not None:
      self.data_opts['maskNSamples'] = data_opts['maskNSamples']
      assert len(self.data_opts['maskNSamples'])==M, "Length of MaskAtRandom and input files does not match"
    else:
      self.data_opts['maskNSamples'] = [0]*M

    # Remove incomplete samples?
    self.data_opts['RemoveIncompleteSamples'] = RemoveIncompleteSamples

  # TODO load other necessary things and also handle the multiple samples
  def load_data(self):
    """ Load the data """

    # Load observations
    self.data, self.sample_groups = loadData(self.data_opts)

    # Remove samples with missing views
    if self.data_opts['RemoveIncompleteSamples']:
      self.data = removeIncompleteSamples(self.data)

    # Calculate dimensionalities
    M = self.dimensionalities["M"]

    # TODO fix for multiple groups
    N = self.dimensionalities["N"] = self.data[0].shape[0]
    D = self.dimensionalities["D"] = [self.data[m].shape[1] for m in range(M)]

    # TODO  covariates_files ?
    self.data_opts['covariates'] = None
    self.data_opts['scale_covariates'] = False

    # TODO change that
    self.all_data = {}
    self.all_data['data'] = self.data
    self.all_data['sample_groups'] = self.sample_groups

  def parse_covariates(self):
    """ Parse covariates """

    K = self.dimensionalities["K"]
    M = self.dimensionalities["M"]

    # Define idx for covariate factors
    idx = range(self.data_opts['covariates'].shape[1])

    # Ignore mean and variance in the prior distribution of Z
    # self.model_opts["priorZ"]["mean"][idx] = s.nan
    self.model_opts["priorZ"]["var"][idx] = s.nan

    # Ignore variance in the variational distribution of Z
    # The mean has been initialised to the covariate values
    self.model_opts["initZ"]["var"][idx] = 0.

  def build_and_run(self):
    model_builder = buildBiofam(self.all_data, self.model_opts, self.model_def)
    self.model = model_builder.net
    self.model.setTrainOptions(self.train_opts)

    train_model(self.model, self.train_opts)

    print("Saving model in %s...\n" % self.data_opts['output_file'])
    self.train_opts['schedule'] = '_'.join(self.train_opts['schedule'])
    saveTrainedModel(model=self.model, outfile=self.data_opts['output_file'], train_opts=self.train_opts, model_opts=self.model_opts,
                     view_names=self.data_opts['view_names'], group_names=self.data_opts['group_names'], sample_groups=self.data_opts['sample_groups'],
                     sample_names=self.data_opts['sample_names'], feature_names=self.data_opts['feature_names'])

if __name__ == '__main__':
  a = entry_point()
  infiles = ["../run/test_data//500_0.txt", "../run/test_data//500_1.txt", "../run/test_data//500_2.txt", "../run/test_data//500_2.txt" ]

  views =  ["view_A", "view_A", "view_B", "view_B"]
  groups = ["group_A", "group_B", "group_A", "group_B"]

  lik = ["gaussian", "gaussian"]

  outfile ="/tmp/test.hdf5"

  a.set_data_options(infiles, outfile, views, groups, delimiter=" ", header_cols=False, header_rows=False)
  a.set_train_options(iter=50, tolerance=0.01, dropR2=0.0)
  a.set_model()
  a.set_model_options(factors=10, likelihoods=lik, learn_intercept=False)
  a.set_dataprocessing_options()
  a.load_data()
  a.build_and_run()
