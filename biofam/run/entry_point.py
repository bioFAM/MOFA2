import argparse
import pandas as pd
import scipy as s
import sys
from time import sleep
from time import time
import pandas as pd
import imp

from typing import List, Union, Dict, TypeVar

from biofam.build_model.build_model import *
from biofam.build_model.train_model import train_model

class entry_point(object):
    def __init__(self):
        self.print_banner()
        self.dimensionalities = {}
        self.model_opts = None
        self.model = None

    def print_banner(self):
        """ Method to print the biofam banner """

        banner = """
         _     _        __
        | |__ (_) ___  / _| __ _ _ __ ___
        | '_ \| |/ _ \| |_ / _` | '_ ` _ \
        | |_) | | (_) |  _| (_| | | | | | |
        |_.__/|_|\___/|_|  \__,_|_| |_| |_|

        """

        print(banner)
        # sleep(2)
        sys.stdout.flush()

    def set_data_matrix(self, data, samples_names_dict, features_names_dict):
        """ Method to set the data

        PARAMETERS
        ----------
        data: several options:
        - a dictionary where each key is the view names and the object is a numpy array or a pandas data frame
        - a list where each element is a numpy array or a pandas data frame
        """

        # Sanity check
        if not isinstance(data, list):
            if isinstance(data, dict):
                data = list(data.values())
            else:
                print("Error: Data not recognised"); sys.stdout.flush(); exit()

        # Convert input data to numpy array format
        for m in range(len(data)):
            if isinstance(data[m], dict):
                data[m] = list(data[m].values())
        for p in range(len(data[0])):
            if not isinstance(data[m][p], np.ndarray):
                if isinstance(data[m][p], pd.DataFrame):
                    data[m][p] = data[m][p].values
                else:
                    print("Error, input data is not a numpy.ndarray or a pandas dataframe"); sys.stdout.flush(); exit()

        # Verbose message
        for m in range(len(data)):
            for p in range(len(data[0])):
                print("Loaded view %d group %d with %d samples and %d features..." % (m+1, p+1, data[m][p].shape[0], data[m][p].shape[1]))

        # Save dimensionalities
        self.dimensionalities["M"] = len(data)
        self.dimensionalities["P"] = len(data[0])
        self.dimensionalities["N"] = [data[0][p].shape[0] for p in range(len(data[0]))]
        self.dimensionalities["D"] = [data[m][0].shape[1] for m in range(len(data))]

        # Define feature group names and sample group names
        self.data_opts['views_names']  = [k for k in features_names_dict.keys()]
        self.data_opts['groups_names'] = [k for k in samples_names_dict.keys()]

        # Define feature and sample names
        self.data_opts['samples_names']  = [v for l in samples_names_dict.values() for v in l]
        self.data_opts['features_names'] = [v for v in features_names_dict.values()]

        # Set samples groups
        # FIX THIS
        self.data_opts['samples_groups'] = [list(samples_names_dict.keys())[i] for i in range(len(self.data_opts['groups_names'])) for n in range(len(list(samples_names_dict.values())[i]))]

        # Concatenate groups
        for m in range(len(data)):
            data[m] = np.concatenate(data[m])
        self.dimensionalities["N"] = np.sum(self.dimensionalities["N"])

        self.data = process_data(data, self.data_opts, self.data_opts['samples_groups'])

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
        # taking unique view and group names in the
        ix = np.unique(self.data_opts['views_names'], return_index=True)[1]
        self.data_opts['views_names'] = (np.array(self.data_opts['views_names'])[ix]).tolist()
        ix = np.unique(self.data_opts['groups_names'], return_index=True)[1]
        self.data_opts['groups_names'] = (np.array(self.data_opts['groups_names'])[ix]).tolist()

        # set dimensionalities of the model
        M = self.dimensionalities['M'] = len(set(self.io_opts['views_names']))
        N = self.dimensionalities["N"] = self.data[0].shape[0]
        D = self.dimensionalities["D"] = [self.data[m].shape[1] for m in range(M)]
        self.dimensionalities['P'] = len(set(self.io_opts['groups_names']))

        # NOTE: Usage of covariates is currently not functional
        self.data_opts['covariates'] = None
        self.data_opts['scale_covariates'] = False

        self.data = process_data(self.data, self.data_opts,  self.samples_groups)

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
        self.data_opts['features_names'] = data.groupby(["feature_group"])["feature"].unique()[self.data_opts['views_names']].tolist()
        self.data_opts['samples_names'] = data["sample"].unique()

        # data_matrix = [None]*len(self.data_opts['views_names'])
        # for m in range(len(self.data_opts['views_names'])):
        #     subdata = data.loc[ data['feature_group'] == self.data_opts['views_names'][m] ]
        #     data_matrix[m] = subdata.pivot(index='sample', columns='feature', values='value')
        # data_matrix = data.pivot(index='sample', columns='feature', values='value').split().tolist()

        # Count the number of feature per view
        tmp = data.groupby(['feature_group'])['feature'].nunique()
        nfeatures = tmp.loc[self.data_opts['views_names']]



        # Convert data frame to list of matrices
        data['feature'] = data['feature'] + data['feature_group'] # make sure there are no duplicated feature names before pivoting
        data_matrix = data.pivot(index='sample', columns='feature', values='value')

        # Sort rows and columns according to the sample and feature names
        features_names_tmp = data.groupby(["feature_group"])["feature"].unique()[self.data_opts['views_names']].tolist()
        data_matrix = data_matrix.loc[self.data_opts['samples_names']]
        data_matrix = data_matrix[[y for x in features_names_tmp for y in x]]

        # Split into a list of views, each view being a matrix
        data_matrix = np.split(data_matrix, np.cumsum(nfeatures)[:-1], axis=1)

        # Define feature names
        # self.data_opts['features_names'] = [ y.columns.values.tolist() for y in data_matrix]

        # Define sample groups
        self.data_opts['samples_groups'] = data[['sample', 'sample_group']].drop_duplicates() \
                                            .set_index('sample').loc[self.data_opts['samples_names']] \
                                            .sample_group.tolist()

        # Define dimensionalities
        self.dimensionalities = {}
        self.dimensionalities['M'] = len(self.data_opts['views_names'])
        self.dimensionalities['N'] = len(self.data_opts['samples_names'])
        self.dimensionalities['P'] = len(self.data_opts['groups_names'])
        self.dimensionalities['D'] = [len(x) for x in self.data_opts['features_names']]

        # Process the data (i.e center, scale, etc.)
        data_matrix = process_data(data_matrix, self.data_opts, self.data_opts['samples_groups'])

        # Convert input data to numpy array format
        for m in range(len(data_matrix)):
            if not isinstance(data_matrix[m], np.ndarray):
                if isinstance(data_matrix[m], pd.DataFrame):
                    data_matrix[m] = data_matrix[m].values
                else:
                    print("Error, input data is not a numpy.ndarray or a pandas dataframe"); exit()

        self.data = data_matrix

        # NOTE: Usage of covariates is currently not functional
        self.data_opts['covariates'] = None
        self.data_opts['scale_covariates'] = False

    def set_train_options(self,
        iter=5000, elbofreq=1, startSparsity=0, tolerance=0.01,
        startDrop=1, freqDrop=1, dropR2=0, nostop=False, verbose=False, seed=None,
        schedule=None, gpu_mode=False
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

        # GPU mode
        self.train_opts['gpu_mode'] = gpu_mode
        if gpu_mode:
            try:
                imp.find_module('cupy')
            except ImportError:
                print('For GPU mode, you need to install the CUPY library')
                print ('1 - Make sure that you are running BIOFAM on a machine with an NVIDIA GPU')
                print ('2 - Install CUPY following instructions on https://docs-cupy.chainer.org/en/stable/install.html')
                print ('Alternatively, deselect GPU mode')
                exit(1)

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
        if seed is None:  # Seed for the random number generator
            seed = int(round(time()*1000)%1e6)
        self.train_opts['seed'] = int(seed)
        s.random.seed(self.train_opts['seed'])


    def set_stochasticity_options(self,
                                  tau=1.,
                                  forgetting_rate=1.,
                                  batch_size=.2):

        # snaity checks
        # assert tau > 0, 'tau must be greater thn zero'
        # assert .5 < forgetting_rate <= 1., 'Choose .5 < forgetting_rate <= 1'
        # assert 0. < batch_size <= 1., 'Choose 0. < batch_size <= 1'

        self.train_opts['stochastic'] = True
        self.train_opts['tau'] = tau
        self.train_opts['forgetting_rate'] = forgetting_rate
        self.train_opts['batch_size'] = batch_size

    def set_model_options(self, factors, likelihoods, sl_z=False, sl_w=False, ard_z=False, ard_w=False):
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
        self.model_opts['noise_on'] = "features"

        # Define initial number of latent factors
        self.dimensionalities["K"] = self.model_opts['factors'] = int(factors)

        # Define likelihoods
        self.model_opts['likelihoods'] = likelihoods
        if isinstance(self.model_opts['likelihoods'], str):
            self.model_opts['likelihoods'] = [self.model_opts['likelihoods']]
        assert len(self.model_opts['likelihoods'])==self.dimensionalities["M"], "Please specify one likelihood for each view"
        assert set(self.model_opts['likelihoods']).issubset(set(["gaussian","bernoulli","poisson"])), "Available likelihoods are 'gaussian','bernoulli' and 'poisson'"

        # Define whether to learn the feature-wise means
        self.model_opts["learn_intercept"] = False

        # Define for which factors and views should we learn the sparsity levels
        self.model_opts['sparsity'] = True
        self.model_opts['learnTheta'] = [s.ones(self.dimensionalities["K"]) for m in range(self.dimensionalities["M"])]

    def set_data_options(self, likelihoods,
        center_features=False, center_features_per_group=False,
        scale_features=False, scale_views=False,
        maskAtRandom=None, maskNSamples=None,
        features_in_rows=False
        ):

        """ Parse data processing options """

        # TODO: more verbose messages
        # TODO Sanity checks
        self.data_opts = {}

        if features_in_rows is True:
            self.data_opts['features_in_rows'] = features_in_rows

        # Define likelihoods
        if type(likelihoods) is not list:
          print("You only specified one likelihood, we assume that you only have a single view")
          likelihoods = [likelihoods]
        assert set(likelihoods).issubset(set(["gaussian","bernoulli","poisson"]))
        self.data_opts["likelihoods"] = likelihoods
        M = len(likelihoods)

        # Center features
        # TO-DO: ITS NOT INTUITIVE TO HARD BOTH CENTER_FEATURES AND CENTER_FEATURES_PER_GROUP, WE NEEED TO FIX THIS
        if center_features_per_group is True:
            self.data_opts['center_features_per_group'] = [ True if l=="gaussian" else False for l in likelihoods ]
            self.data_opts['center_features'] = [ False for l in likelihoods ]
        elif center_features is True:
            self.data_opts['center_features_per_group'] = [ False for l in likelihoods ]
            self.data_opts['center_features'] = [ True if l=="gaussian" else False for l in likelihoods ]
        else:
            # if not self.model_opts["learn_intercept"]: print("\nWarning... you are not centering the data and not learning the mean...\n")
            self.data_opts['center_features'] = [ False for l in likelihoods ]
            self.data_opts['center_features_per_group'] = [ False for l in likelihoods ]


        # Scale views
        if scale_views is True:
            self.data_opts['scale_views'] = [ True if l=="gaussian" else False for l in likelihoods ]
        else:
            self.data_opts['scale_views'] = [ False for l in likelihoods ]

        # Scale features
        if scale_features is True:
            assert self.data_opts['scale_views'][0] is False, "Scale either entire views or features, not both"
            self.data_opts['scale_features'] = [ True if l=="gaussian" else False for l in likelihoods ]
        else:
            self.data_opts['scale_features'] = [ False for l in likelihoods ]


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

    def run(self, no_theta=False):
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
          data=self.data,
          train_opts=self.train_opts,
          model_opts=self.model_opts,
            samples_names=self.data_opts['samples_names'],
          features_names=self.data_opts['features_names'],
          views_names=self.data_opts['views_names'],
            groups_names=self.data_opts['groups_names'],
          samples_groups=self.data_opts['samples_groups']
        )

if __name__ == '__main__':

    ent = entry_point()

    infiles = ["../run/test_data/with_nas/500_0.txt", "../run/test_data/with_nas/500_1.txt", "../run/test_data/with_nas/500_2.txt", "../run/test_data/with_nas/500_2.txt" ]
    views =  ["view_A", "view_A", "view_B", "view_B"]
    groups = ["group_A", "group_B", "group_A", "group_B"]

    lik = ["gaussian", "gaussian"]
    out_file = '/tmp/test_biofam.hdf5'
    ent.set_data_options(lik, center_features=False, center_features_per_group=False, scale_features=False, scale_views=False)
    ent.set_data_from_files(infiles, views, groups, delimiter=" ", header_cols=False, header_rows=False)
    ent.set_model_options(ard_z=True, sl_w=True , sl_z=True, ard_w=True, factors=15, likelihoods=lik)
    ent.set_train_options(iter=10, tolerance=1., dropR2=0.0, seed=4, elbofreq=1, verbose=1)
    ent.build()
    ent.run()
    ent.save(out_file)
