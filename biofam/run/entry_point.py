import argparse
import pandas as pd
import scipy as s
import sys
from time import sleep
from time import time
import pandas as pd

from typing import List, Union, Dict, TypeVar

from biofam.build_model.build_model import *
from biofam.build_model.train_model import train_model

# TODO where do we input the out_file option ??
class entry_point(object):
    def __init__(self):
        self.print_banner()
        self.dimensionalities = {}
        self.model_opts = None
        self.model = None

    def print_banner(self):
        """ Method to print the biofam banner """

        banner = r"""
         _     _        __
        | |__ (_) ___  / _| __ _ _ __ ___
        | '_ \| |/ _ \| |_ / _` | '_ ` _ \
        | |_) | | (_) |  _| (_| | | | | | |
        |_.__/|_|\___/|_|  \__,_|_| |_| |_|

        """

        print(banner)
        # sleep(2)
        sys.stdout.flush()
    
    def set_data_matrix(self, data):
        """ Method to set the data 

        PARAMETERS
        ----------
        data: several options:
        - a dictionary where each key is the view names and the object is a numpy array or a pandas data frame
        - a list where each element is a numpy array or a pandas data frame
        """

        # Sanity check
        if isinstance(data, dict):
          data = list(data.values())
        elif isinstance(data, list):
          pass
        else:
          print("Error: Data not recognised")
          sys.stdout.flush()
          exit()

        for m in range(len(data)):
          if isinstance(data[m], dict):
            data[m] = list(data[m].values())

        assert s.all([ isinstance(data[m][p], s.ndarray) or isinstance(data[m][p], pd.DataFrame) for p in range(len(data[0])) for m in range(len(data)) ]), "Error, input data is not a numpy.ndarray"

        # Verbose message
        for m in range(len(data)):
          for p in range(len(data[0])):
            print("Loaded view %d group %d with %d samples and %d features..." % (m+1, p+1, data[m][p].shape[0], data[m][p].shape[1]))

        # Save dimensionalities
        self.dimensionalities["M"] = len(data)
        self.dimensionalities["P"] = len(data[0])
        self.dimensionalities["N"] = [data[0][p].shape[0] for p in range(len(data[0]))]
        self.dimensionalities["D"] = [data[m][0].shape[1] for m in range(len(data))]

        self.data = data

    def set_data_from_files(self, inFiles, views, groups, header_rows=False, header_cols=False, delimiter=' '):
        """ Load the data """

        # TODO: sanity checks
        # - Check that input file exists
        # - Check that output directory exists: warning if not
        # TODO: Verbosity, print messages

        self.io_opts = {}

        # I/O
        if type(inFiles) is not list:
            inFiles = [inFiles]
        self.io_opts['input_files'] = inFiles
        self.io_opts['delimiter'] = delimiter

        # View names and group names (sample groups)
        if type(views) is not list:
            views = [views]
        if type(groups) is not list:
            groups = [groups]

        self.io_opts['views_names'] = views
        self.io_opts['groups_names'] = groups

        # Headers
        if header_rows is True:
            self.io_opts['rownames'] = 0
        else:
            self.io_opts['rownames'] = None

        if header_cols is True:
            self.io_opts['colnames'] = 0
        else:
            self.io_opts['colnames'] = None

        # Load observations
        self.data_opts.update(self.io_opts)
        # TODO reindex
        self.data, self.samples_groups = loadData(self.data_opts)

        # data frame needs reindexing if no header row
        if not header_rows:
            self.data[0] = self.data[0].reset_index(drop=True)
        # save feature, sample names, sample groups

        self.data_opts['samples_names'] = self.data[0].index
        self.data_opts['features_names'] = [dt.columns.values for dt in self.data]
        # TODO check that we have the right dictionary
        # TODO check that what is used later in the code is ok for this
        self.data_opts['samples_groups'] = self.samples_groups

        # set dimensionalities of the model
        M = self.dimensionalities['M'] = len(set(self.io_opts['views_names']))
        N = self.dimensionalities["N"] = self.data[0].shape[0]
        D = self.dimensionalities["D"] = [self.data[m].shape[1] for m in range(M)]
        self.dimensionalities['P'] = len(set(self.io_opts['groups_names']))

        # NOTE: Usage of covariates is currently not functional
        self.data_opts['covariates'] = None
        self.data_opts['scale_covariates'] = False

    def set_data_df(self, data):
        """Method to define the data

        PARAMETERS
        ----------
        data: pd.DataFrame
            a pandas DataFrame with columns ("sample","sample_group","feature","feature_group","value")
            the order is irrelevant
        """

        # Sanity checks
        assert hasattr(self, 'data_opts'), "Data options not defined"
        assert isinstance(data, pd.DataFrame), "'data' has to be an instance of pd.DataFrame"
        assert 'sample' in data.columns, "'data' has to contain the column 'sample'"
        assert 'sample_group' in data.columns, "'data' has to contain the column 'sample_group'"
        assert 'feature' in data.columns, "'data' has to contain the column 'feature'"
        assert 'feature_group' in data.columns, "'data' has to contain the column 'feature_group'"
        assert 'value' in data.columns, "'data' has to contain the column 'value'"

        # Define feature group names and sample group names
        self.data_opts['views_names'] = data["feature_group"].unique()
        self.data_opts['groups_names'] = data["sample_group"].unique()

        # Define feature and sample names
        self.data_opts['samples_names'] = data["sample"].unique()
        self.data_opts['features_names'] = [ data.loc[data['feature_group'] == m].feature.unique() for m in self.data_opts['views_names'] ]

        # (DEPRECIATED) Create dictionaries with mapping between:
        #  sample names (keys) and samples_groups (values)
        #  feature names (keys) and feature groups (values)
        # TODO CHECK that the order is fine here ...
        # self.data_opts['samples_groups'] = pd.Series(data["sample_group"].values, index=data["sample"].values).to_dict()
        # self.data_opts['feature_groups'] = pd.Series(data.feature_group.values, index=data.feature).to_dict()

        # Define sample groups
        self.data_opts['samples_groups'] = data[['sample', 'sample_group']].drop_duplicates() \
                                            .set_index('sample').loc[self.data_opts['samples_names']] \
                                            .sample_group.tolist()

        # Define dictionary with the dimensionalities
        self.dimensionalities = {}
        self.dimensionalities['D'] = [len(x) for x in self.data_opts['features_names']]
        self.dimensionalities['M'] = len(self.data_opts['views_names'])
        self.dimensionalities['N'] = len(self.data_opts['samples_names'])
        self.dimensionalities['P'] = len(self.data_opts['groups_names'])

        # Convert data frame to list of matrices
        data_matrix = [None]*self.dimensionalities['M']
        for m in range(self.dimensionalities['M']):
            subdata = data.loc[ data['feature_group'] == self.data_opts['views_names'][m] ]
            data_matrix[m] = subdata.pivot(index='sample', columns='feature', values='value')

        # Process the data (i.e center, scale, etc.)
        self.data = process_data(data_matrix, self.data_opts, self.data_opts['samples_groups'])

        # NOTE: Usage of covariates is currently not functional
        self.data_opts['covariates'] = None
        self.data_opts['scale_covariates'] = False

    def set_train_options(self,
        iter=5000, elbofreq=1, startSparsity=0, tolerance=0.01,
        startDrop=1, freqDrop=1, dropR2=0, nostop=False, verbose=False, seed=None,
        schedule=None
        ):
        """ Parse training options """

        # TODO: verbosity, print more messages

        self.train_opts = {}

        # Maximum number of iterations
        self.train_opts['maxiter'] = int(iter)

        # Lower bound computation frequency
        self.train_opts['elbofreq'] = int(elbofreq)

        # Verbosity
        self.train_opts['verbose'] = verbose

        # Minimum Variance explained threshold to drop inactive factors
        self.train_opts['drop'] = { "by_r2":float(dropR2) }
        self.train_opts['start_drop'] = int(startDrop)
        self.train_opts['freq_drop'] = int(freqDrop)
        if (dropR2>0 & verbose==True): print("\nDropping factors with minimum threshold of {0}% variance explained".format(dropR2))

        # Tolerance level for convergence
        self.train_opts['tolerance'] = float(tolerance)

        # Do no stop even when convergence criteria is met
        self.train_opts['forceiter'] = nostop

        # Iteration to activate spike and slab sparsity
        self.train_opts['start_sparsity'] = int(startSparsity)

        # Training schedule
        if schedule is not None:
            self.train_opts['schedule'] = schedule

        # Seed
        if seed is None or seed == 0:  # Seed for the random number generator
            self.train_opts['seed'] = int(round(time() * 1000) % 1e6)
        else:
            self.train_opts['seed'] = seed
        s.random.seed(self.train_opts['seed'])

    def set_model_options(self, factors, likelihoods, sl_z=False, sl_w=False, ard_z=False, ard_w=False, noise_on='features', learnTheta=True, learn_intercept=False):
        """ Set model options """

        # TODO: SANITY CHECKS AND:
        # - learnTheta should be replaced by learn_sparsity


        if self.model_opts is None:
            self.model_opts = {}

        # Define whether to use sample-wise spike and slab prior for Z
        self.model_opts['sl_z'] = sl_z

        # Define whether to use feature-wise spike and slab prior for W
        self.model_opts['sl_w'] = sl_w

        # Define whether to use sample_group and factor-wise ARD prior for Z
        self.model_opts['ard_z'] = ard_z

        # Define whether to use view and factor-wise ARD prior for W
        self.model_opts['ard_w'] = ard_w

        # Define whether to add noise terms on the features or on the samples
        self.model_opts['noise_on'] = noise_on

        # Define initial number of latent factors
        self.dimensionalities["K"] = self.model_opts['factors'] = int(factors)

        # Define likelihoods
        self.model_opts['likelihoods'] = likelihoods
        if isinstance(self.model_opts['likelihoods'], str):
            self.model_opts['likelihoods'] = [self.model_opts['likelihoods']]

        assert len(self.model_opts['likelihoods'])==self.dimensionalities["M"], "Please specify one likelihood for each view"
        assert set(self.model_opts['likelihoods']).issubset(set(["gaussian","bernoulli","poisson"])), "Available likelihoods are 'gaussian','bernoulli' and 'poisson'"

        # Define whether to learn the feature-wise means
        self.model_opts["learn_intercept"] = learn_intercept
        if learn_intercept:
            self.model_opts['factors'] += 1
            self.dimensionalities["K"] += 1

        # Define for which factors and views should we learn the sparsity levels
        if isinstance(learnTheta, bool):
            self.model_opts['sparsity'] = True
            self.model_opts['learnTheta'] = [s.ones(self.dimensionalities["K"]) for m in range(self.dimensionalities["M"])]
        elif isinstance(learnTheta,list):
        	print("Depreciated, '--learn-theta' has to be a boolean")
			# self.model_opts['sparsity'] = True
			# assert len(learnTheta)==M, "--learnTheta has to be a binary vector with length number of views"
			# self.model_opts['learnTheta'] = [ learnTheta[m]*s.ones(K) for m in range(M) ]
        else:
            print("Error, --learn-theta has to be a boolean")
            exit(1)

    def set_data_options(self, lik,
        center_features=False, center_features_per_group=False,
        scale_features=False, scale_views=False,
        maskAtRandom=None, maskNSamples=None,
        features_in_rows=False
        ):

        """ Parse data processing options """

        # RICARD: LIKELIHOOD SHOULD BE IN MODEL_OPTS, NOT IN DATA_OPTIONS
        #       WHY IS SELF_MODEL.OPTS() DEFINED HERE??????

        # TODO: more verbose messages
        # TODO Sanity checks
        self.data_opts = {}
        self.model_opts = {}

        # Define likelihoods
        self.model_opts['likelihoods'] = lik
        if type(self.model_opts['likelihoods']) is not list:
          self.model_opts['likelihoods'] = [self.model_opts['likelihoods']]
        # assert len(self.model_opts['likelihoods'])==M, "Please specify one likelihood for each view"
        assert set(self.model_opts['likelihoods']).issubset(set(["gaussian","bernoulli","poisson"]))
        self.data_opts["likelihoods"] = self.model_opts['likelihoods'] 
        
        M = len(self.model_opts["likelihoods"])
        # assert len(self.data_opts['views_names'])==M, "Length of view names and input files does not match"

        if features_in_rows is True:
            self.data_opts['features_in_rows'] = features_in_rows

        # Center features
        # TO-DO: ITS NOT INTUITIVE TO HARD BOTH CENTER_FEATURES AND CENTER_FEATURES_PER_GROUP, WE NEEED TO FIX THIS
        if center_features_per_group is True:
            self.data_opts['center_features_per_group'] = [ True if l=="gaussian" else False for l in self.model_opts["likelihoods"] ]
            self.data_opts['center_features'] = [ False for l in self.model_opts["likelihoods"] ]
        elif center_features is True:
            self.data_opts['center_features_per_group'] = [ False for l in self.model_opts["likelihoods"] ]
            self.data_opts['center_features'] = [ True if l=="gaussian" else False for l in self.model_opts["likelihoods"] ]
        else:
            # if not self.model_opts["learn_intercept"]: print("\nWarning... you are not centering the data and not learning the mean...\n")
            self.data_opts['center_features'] = [ False for l in self.model_opts["likelihoods"] ]
            self.data_opts['center_features_per_group'] = [ False for l in self.model_opts["likelihoods"] ]


        # Scale views
        if scale_views is True:
            self.data_opts['scale_views'] = [ True if l=="gaussian" else False for l in self.model_opts["likelihoods"] ]
        else:
            self.data_opts['scale_views'] = [ False for l in self.model_opts["likelihoods"] ]

        # Scale features
        if scale_features is True:
            assert self.data_opts['scale_views'][0] is False, "Scale either entire views or features, not both"
            self.data_opts['scale_features'] = [ True if l=="gaussian" else False for l in self.model_opts["likelihoods"] ]
        else:
            self.data_opts['scale_features'] = [ False for l in self.model_opts["likelihoods"] ]


        # Mask values at random
        if maskAtRandom is not None:
            self.data_opts['maskAtRandom'] = data_opts['maskAtRandom']
            assert len(self.data_opts['maskAtRandom'])==M, "Length of MaskAtRandom and input files does not match"
        else:
            self.data_opts['maskAtRandom'] = [0]*M

       # Mask entire samples
        if maskNSamples is not None:
            self.data_opts['maskNSamples'] = data_opts['maskNSamples']
            assert len(self.data_opts['maskNSamples'])==M, "Length of MaskAtRandom and input files does not match"
        else:
            self.data_opts['maskNSamples'] = [0]*M

    def build(self):
        """ Build the model """
        # Sanity checks
        assert hasattr(self, 'model_opts'), "Train options not defined"
        assert hasattr(self, 'dimensionalities'), "Dimensionalities are not defined"

        # Build the BioFAM model
        self.model_builder = buildBiofam(self.data, self.data_opts, self.model_opts, self.dimensionalities)
        self.model = self.model_builder.net

    def run(self):
        """ Run the model """

        # Sanity checks
        assert hasattr(self, 'train_opts'), "Train options not defined"
        assert hasattr(self, 'data_opts'), "Data options not defined"

        # Fetch training schedule (order of updates for the different nodes)
        if 'schedule' in self.train_opts:
            assert all(self.train_opts['schedule'] in self.model.get_nodes().keys()), "Some nodes defined in the training schedule are not present in the model"
            if ~all(self.model.get_nodes().keys() in self.train_opts['schedule']):
                if self.train_opts['verbose']: print("Warning: some nodes are not in the trainign schedule and will not be updated")
        else:
            self.train_opts['schedule'] = self.model_builder.schedule

        self.model.setTrainOptions(self.train_opts)

        # Train the model
        train_model(self.model, self.train_opts)

    def save(self, outfile):
        """ Save the model in an hdf5 file """

        if self.train_opts["verbose"]:
            print("Saving model in %s...\n" % outfile)

        self.train_opts['schedule'] = '_'.join(self.train_opts['schedule'])

        saveTrainedModel(
            model=self.model,
            outfile=outfile,
            train_opts=self.train_opts,
            model_opts=self.model_opts,
            samples_names=self.data_opts['samples_names'],
            features_names=self.data_opts['features_names'],
            views_names=self.data_opts['views_names'],
            groups_names=self.data_opts['groups_names'],
            samples_groups=self.data_opts['samples_groups']
        )