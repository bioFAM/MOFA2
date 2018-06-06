import argparse
import pandas as pd
import scipy as s
import sys
from time import sleep
from time import time
import pandas as pd

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

    def set_data(self, data):
        """Method to define the data

        PARAMETERS
        ----------
        data: pd.DataFrame
        	a pandas DataFrame with columns ("sample","sample_group","feature","feature_group","value")
        	the order is irrelevant
        """

        # Sanity checks
        assert isinstance(data, pd.DataFrame), "'data' has to be an instance of pd.DataFrame"
        assert 'sample' in data.columns, "'data' has to contain the column 'sample'"
        assert 'sample_group' in data.columns, "'data' has to contain the column 'sample_group'"
        assert 'feature' in data.columns, "'data' has to contain the column 'feature'"
        assert 'feature_group' in data.columns, "'data' has to contain the column 'feature_group'"
        assert 'value' in data.columns, "'data' has to contain the column 'value'"

        # Defien data options
        self.data_opts = {}

        # Define feature groups and sample groups
        self.data_opts['view_names'] = data["feature_group"].unique()
        self.data_opts['group_names'] = data["sample_group"].unique()

        # Define feature and sample names
        self.data_opts['sample_names'] = data["sample"].unique()
        self.data_opts['feature_names'] = [ data.loc[data['feature_group'] == m].feature.unique() for m in self.data_opts['view_names'] ]

        # Create dictionaries with mapping between:
        #  sample names (keys) and sample_groups (values)
        #  feature names (keys) and feature groups (values)
        # TODO CHECK that the order is fine here ...
        self.data_opts['sample_groups'] = pd.Series(df.sample_group.values, index=df.sample).values()
        # self.data_opts['feature_groups'] = pd.Series(df.feature_group.values, index=df.feature).to_dict()

        # Define dictionary with the dimensionalities
        self.dimensionalities = {}
        self.dimensionalities['D'] = [len(x) for x in self.data_opts['feature_names']]
        self.dimensionalities['M'] = len(self.data_opts['view_names'])
        self.dimensionalities['N'] = len(self.data_opts['sample_names'])
        self.dimensionalities['P'] = len(self.data_opts['group_names'])

        # Convert data frame to nested list of matrices where
        # the first level splits by feature_group (views) and the second level splits by sample_group
        data_matrix = [[None]*self.dimensionalities['P'] for m in range(self.dimensionalities['M'])]
        for m in range(self.dimensionalities['M']):
            for p in range(self.dimensionalities['P']):
                subdata = data.loc[(data['feature_group'] == self.data_opts['view_names'][m]) & (data['sample_group'] == self.data_opts['group_names'][p]) ]
                data_matrix[m][p] = subdata.pivot(index='sample', columns='feature', values='value').values
                # TO-DO: Reorder to match feature_names and sample_names
                # NOT TESTED data_matrix[m][p] = data_matrix[m][p].reindex(self.data_opts['sample_names'])
                # NOT TESTED data_matrix[m][p] = data_matrix[m][p][self.data_opts['feature_names'][m]]

        # TODO check that
        self.data = process_data(data_matrix, self.data_opts, self.data_opts['sample_groups'])
        # NOTE: Usage of covariates is currently not functional
        self.data_opts['covariates'] = None
        self.data_opts['scale_covariates'] = False

    def set_data_from_files(self, inFiles, views, groups, header_rows=False,
                        header_cols=False, delimiter=' '):
        """ Load the data """
                # inFiles, outFile, views, groups,
                # delimiter=" ", header_cols=False, header_rows=False

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
            views = [groups]

        self.io_opts['view_names'] = views
        self.io_opts['group_names'] = groups

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
        self.data, self.sample_groups = loadData(self.data_opts)

        # data frame needs reindexing if no header row
        if not header_rows:
            self.data[0] = self.data[0].reset_index(drop=True)
        # save feature, sample names, sample groups

        self.data_opts['sample_names'] = self.data[0].index
        self.data_opts['feature_names'] = [dt.columns.values for dt in self.data]
        # TODO check that we have the right dictionary
        # TODO check that what is used later in the code is ok for this
        self.data_opts['sample_groups'] = self.sample_groups

        # set dimensionalities of the model
        M = self.dimensionalities['M'] = len(set(self.io_opts['view_names']))
        N = self.dimensionalities["N"] = self.data[0].shape[0]
        D = self.dimensionalities["D"] = [self.data[m].shape[1] for m in range(M)]
        self.dimensionalities['P'] = len(set(self.io_opts['group_names']))

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


    def set_model_options(self,factors, likelihoods,
    	sl_z=False, sl_w=False, ard_z=False, ard_w=False, noise_on='features',
    	learnTheta=True, learn_intercept=False):
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
        if isinstance(self.model_opts['likelihoods'],str):
            self.model_opts['likelihoods'] = [self.model_opts['likelihoods']]

        assert len(self.model_opts['likelihoods'])==self.dimensionalities["M"], "Please specify one likelihood for each view"
        assert set(self.model_opts['likelihoods']).issubset(set(["gaussian","bernoulli","poisson"])), "Available likelihoods are 'gaussian','bernoulli' and 'poisson'"

        # Define whether to learn the feature-wise means
        self.model_opts["learn_intercept"] = learn_intercept
        if learn_intercept:
            self.model_opts['factors'] += 1
            self.dimensionalities["K"] += 1

        # Define for which factors and views should we learn the sparsity levels
        if isinstance(learnTheta,bool):
            self.model_opts['sparsity'] = True
            self.model_opts['learnTheta'] = [s.ones(self.dimensionalities["K"]) for m in range(self.dimensionalities["M"])]
        elif isinstance(learnTheta,list):
        	print("Depreciated, '--learnTheta' has to be a boolean")
			# self.model_opts['sparsity'] = True
			# assert len(learnTheta)==M, "--learnTheta has to be a binary vector with length number of views"
			# self.model_opts['learnTheta'] = [ learnTheta[m]*s.ones(K) for m in range(M) ]
        else:
            print("Error, --learnTheta has to be a boolean")
            exit(1)

        # TODO sort that out
        self.data_opts['features_in_rows'] = False

    def set_data_options(self, lik,
        center_features=False, center_features_per_group=False,
        scale_features=False, scale_views=False,
        maskAtRandom=None, maskNSamples=None
        ):

        """ Parse data processing options """

        # TODO: more verbose messages
        # TODO Sanity checks
        self.data_opts = {}
        self.model_opts = {}

        self.data_opts["likelihoods"] = lik
        self.model_opts["likelihoods"] = lik
        M = len(self.model_opts["likelihoods"])
        # assert len(self.data_opts['view_names'])==M, "Length of view names and input files does not match"

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
            assert data_opts['scale_views'] is False, "Scale either entire views or features, not both"
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

        if no_theta:
            self.train_opts['schedule'].remove('ThetaW')
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
	    	sample_names=self.data_opts['sample_names'],
            feature_names=self.data_opts['feature_names'],
            view_names=self.data_opts['view_names'],
	    	group_names=self.data_opts['group_names'],
            sample_groups=self.data_opts['sample_groups']
        )



    def get_df(self, node):
        assert self.model is not None, 'Model is not built yet'

        nodes = self.model.nodes
        assert node in nodes, "requested node is not in the model"

        sample_names = self.data_opts['sample_names']
        feature_names = self.data_opts['feature_names']
        factor_names = np.array(range(nodes['Z'].K)).astype(str)
        view_names = np.unique(self.data_opts['view_names'])
        sample_groups=self.data_opts['sample_groups']

        # TODO this index_opts could be passed as an argument of a get_df method implemented in each nodes.
        # then each node would know what index option to use to build the dataframe
        # index_opts = {'sample_names': sample_names,
        #               'feature_names': feature_names,
        #               'view_names': view_names,
        #               'sample_groups':sample_groups,
        #               'factor_names':factor_names}
        exp = nodes[node].getExpectation()

        if node == 'W':
            i=0
            all_dfs = []
            for view in view_names:
                e = pd.DataFrame(exp[i], index=feature_names[i], columns=factor_names)
                e.index.name = 'feature'
                e = e.reset_index(drop=True)
                e_melted = pd.melt(e, id_vars=['feature'], var_name='factor', value_name='value')
                e_melted['view'] = view

                all_dfs.append(e_melted)

                i += 1

            res = pd.concat(all_dfs)

        if node == 'Y':
            i=0
            all_dfs = []
            for view in view_names:
                e = pd.DataFrame(exp[i], index=sample_names, columns=feature_names[i])
                e['group']=sample_groups
                e.index.name = 'sample'
                e = e.reset_index(drop=True)
                e_melted = pd.melt(e, id_vars=['group', 'sample'], var_name='feature', value_name='value')
                e_melted['view'] = view

                all_dfs.append(e_melted)

                i += 1

            res = pd.concat(all_dfs)

        if node == 'Z':
            e = pd.DataFrame(exp, index=sample_names, columns=factor_names)
            e['group']=sample_groups
            e.index.name = 'sample'
            e = e.reset_index(drop=True)
            e_melted = pd.melt(e, id_vars=['group', 'sample'], var_name='factor', value_name='value')

            res = e_melted

        return res

    def depreciated_parse_covariates(self):
        """ Parse covariates """
        print("Covariates are not implemented")
        exit()

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


###### IGNORE BELOW ########



class entry_sfa(entry_point):
    def __init__(self):
        super(entry_sfa, self).__init__()

    def set_data_options(self,
      inFiles, outFile, views, groups, x_files, view_has_covariance_prior=None,
      delimiter=" ", header_cols=False, header_rows=False
      ):
        super(entry_sfa, self).set_data_options(inFiles, outFile, views, groups,
          delimiter, header_cols, header_rows)
        # Position files
        self.data_opts['x_files'] = x_files

        # Boolean to say whether the covariance prior is known for each view
        # Only used if covariance prior on W
        if view_has_covariance_prior is None:
            view_has_covariance_prior = [True] * self.dimensionalities['M']
        self.data_opts['view_has_covariance_prior'] = view_has_covariance_prior


    def set_model(self, sl_z=False, sl_w=True, ard_z=False, ard_w=True, noise_on='features', cov_on='samples'):
        # sanity check for the current implementation
        assert (cov_on == samples) or (cov_on == 'features'), 'cov_on argument not understood'

        if cov_on == 'samples':
            assert not sl_z, 'cannot put spike and slab on Z if covariance on samples'
            assert not ard_z, 'cannot put ARD on Z if covariance on samples'

        if cov_on == 'features':
            assert not sl_z, 'cannot put spike and slab on W if covariance on features'
            assert not ard_z, 'cannot put ARD on W if covariance on features'

        # define where you put the covariance prior
        self.model_opts['cov_on'] = cov_on
        super(entry_sfa, self).set_model(sl_z, sl_w, ard_z, ard_w, noise_on)

    def load_data(self):
        super(entry_sfa, self).load_data()
        self.all_data['data_x'] = loadDataX(self.data_opts)

    def build_and_run(self):
        model_builder = buildSpatialBiofam(self.all_data, self.model_opts)

        self.model = model_builder.net
        self.train_opts['schedule'] = model_builder.schedule
        self.model.setTrainOptions(self.train_opts)

        train_model(self.model, self.train_opts)

        print("Saving model in %s...\n" % self.data_opts['output_file'])
        self.train_opts['schedule'] = '_'.join(self.train_opts['schedule'])
        saveTrainedModel(model=self.model, outfile=self.data_opts['output_file'], train_opts=self.train_opts, model_opts=self.model_opts,
                         view_names=self.data_opts['view_names'], group_names=self.data_opts['group_names'], sample_groups=self.all_data['sample_groups'])



if __name__ == '__main__':
    ent = entry_point()

    infiles = ["../run/test_data/with_nas/500_0.txt", "../run/test_data/with_nas/500_1.txt", "../run/test_data/with_nas/500_2.txt", "../run/test_data/with_nas/500_2.txt" ]
    views =  ["view_A", "view_A", "view_B", "view_B"]
    groups = ["group_A", "group_B", "group_A", "group_B"]

    # views =  ["view_0"]
    # groups = ["group_0"]

    lik = ["gaussian", "gaussian"]
    # lik = ["gaussian"]
    #
    # outfile = dir+"test_no_sl.hdf5"
    #
    ent.set_data_options(lik, center_features=True, center_features_per_group=False, scale_features=False, scale_views=True)
    ent.set_data_from_files(infiles, views, groups, delimiter=" ", header_cols=False, header_rows=False)
    ent.set_model_options(ard_z=True, sl_w=True, sl_z=False, ard_w=True, factors=15, likelihoods=lik, learnTheta=False)
    ent.set_train_options(iter=10, tolerance=0.01, dropR2=0.0, seed=1, elbofreq=1)

    ent.build()
    ent.run(no_theta=False)
    # ent.save(outfile)
    #
    # outfile2 = dir+"test_sl.hdf5"
    # ent2 = entry_point()
    # ent2.set_data_options(lik, center_features=True, center_features_per_group=False, scale_features=False, scale_views=True)
    # ent2.set_data_from_files(infiles, views, groups, delimiter=" ", header_cols=False, header_rows=False)
    # ent2.set_model_options(ard_z=False, sl_w=True, sl_z=False, ard_w=True, factors=4, likelihoods=lik, learnTheta=False)
    # ent2.set_train_options(iter=500, tolerance=0.01, dropR2=0.0)
    # ent2.build()
    # ent2.run(no_theta=False)
    # ent2.save(outfile2)
    # ent.get_df('Y')


    # # from biofam.run.entry_point import entry_point
    # file = "/Users/ricard/Downloads/test_biofam/data.txt"
    # lik = ["gaussian", "gaussian"]
    # ent.set_data_options()
    # data = pd.read_csv(file, delimiter="\t")
    # ent = entry_point()
    # ent.set_data(data)
    # ent.set_train_options(iter=10, tolerance=0.01, dropR2=0.0)
    # ent.set_model(sl_z=False, sl_w=False, ard_z=True, ard_w=False, noise_on='features')
    # ent.set_model_options(factors=10, likelihoods=lik)
    # ent.build_and_run()
    # ent.get_df('Y')
