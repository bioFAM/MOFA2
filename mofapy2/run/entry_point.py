import numpy as np
import pandas as pd
import scipy as s
import sys
import os
from time import sleep, time, strftime
from typing import List, Optional, Union
from itertools import chain
from functools import wraps

from mofapy2.core.BayesNet import *
from mofapy2.core import gpu_utils
from mofapy2.build_model.build_model import *
from mofapy2.build_model.save_model import *
from mofapy2.build_model.utils import guess_likelihoods
from mofapy2.build_model.train_model import train_model
# import matplotlib.pyplot as plt



def keyboardinterrupt_saver(func):
    @wraps(func)
    def saver(self, *args, **kwargs):
        try:
            func(self, *args, **kwargs)
        # Internal methods will raise TypeError when interrupted
        except (KeyboardInterrupt, TypeError):
            if self.train_opts["save_interrupted"]:
                print("Attempting to save the model at the current iteration...")
                if self.train_opts["outfile"] is not None and self.train_opts != "":
                    tmp_file = self.train_opts["outfile"]
                    tmp_file = tmp_file.rstrip(".hdf5") + "_interrupted.hdf5"
                else:
                    tmp_file = os.path.join('/tmp', "mofa_{}_interrupted.hdf5".format(strftime('%Y%m%d-%H%M%S')))
                self.save(outfile=tmp_file)
                print("Saved partially trained model in {}. Exiting now.".format(tmp_file))
            else:
                print("Exiting now without saving the partially trained model. To save a partially trained model, set save_interrupted in the training options to true.")
            sys.stdout.flush()
            sys.exit()
    return saver

class entry_point(object):
    def __init__(self):
        self.print_banner()
        self.dimensionalities = {"C":0}
        self.model = None
        self.imputed = False # flag
        self.Zcompleted = False # flag

    def print_banner(self):
        """ Method to print the mofapy2 banner """
        
        banner = """
        #########################################################
        ###           __  __  ____  ______                    ### 
        ###          |  \/  |/ __ \|  ____/\    _             ### 
        ###          | \  / | |  | | |__ /  \ _| |_           ### 
        ###          | |\/| | |  | |  __/ /\ \_   _|          ###
        ###          | |  | | |__| | | / ____ \|_|            ###
        ###          |_|  |_|\____/|_|/_/    \_\              ###
        ###                                                   ### 
        ######################################################### 
       \n 
        """

        print(banner)
        sys.stdout.flush()

    def set_covariates(self, sample_cov, covariates_names=None):
        """"
        sample_cov: a list of matrices per group
                    The dimensions of each matrix must be (samples, covariates)
        			The order of list elements and rows in each matrix must match the structure of data
        covariates_names: XXXXX
        """

        assert self.data is not None, "the data must be specified with set_data before specifying covariates"

        self.smooth_opts = {}    
        self.smooth_opts["GP_factors"] = True

        # Define covariates
        # Define smooth options
        # if not hasattr(self, 'smooth_opts'): 
        #     print("Options for the usage of smooth factors not defined, using default values...\n")
        #     self.set_smooth_options()

        # Sanity check
        if not isinstance(sample_cov, list):
            if isinstance(sample_cov, dict):
                sample_cov = list(sample_cov.values())
            # if providing a single matrix, treat it as G=1
            elif isinstance(sample_cov, pd.DataFrame):
                sample_cov = [sample_cov.values]
            elif isinstance(sample_cov, np.ndarray):
                sample_cov = [sample_cov]
            else:
                print("Error: sample_cov not recognised");
                sys.stdout.flush();
                sys.exit()

        assert len(sample_cov) == self.dimensionalities["G"], "sample_cov needs to be a list of same length as data (same number of groups)"

        for g in range(self.dimensionalities["G"]):
            if not isinstance(sample_cov[g], np.ndarray):
                if isinstance(sample_cov[g], pd.DataFrame):
                    sample_cov[g] = sample_cov[g].values
                else:
                    print("Error, sample_cov is not a numpy.ndarray or a pandas dataframe"); sys.stdout.flush(); sys.exit()
            sample_cov[g] = sample_cov[g].astype(np.float64)

        N = [len(x) for x in self.data_opts['samples_names']]
        if not all([sample_cov[g].shape[0] == N[g] for g in range(self.dimensionalities["G"])]):
            print("Error, number of rows in sample covariates does not match number of samples in input data (N=%d vs. N=%d)" % ([sample_cov[g].shape[0] for g in range(G)], N))
            sys.stdout.flush(); sys.exit()

        # concatenate groups in sample_cov and standardize sample_cov to avoid scale differences
        sample_cov = np.concatenate(sample_cov, axis = 0)

        # Define dimensionality
        self.dimensionalities["C"] = sample_cov.shape[1]

        # If sample_cov loaded successfully, print verbose message
        if not sample_cov is None:
            print("Loaded %d covariate(s) for each sample..." % (sample_cov.shape[1]))
            print("\n")

        # Define covariate names
        if covariates_names is None:
            print("Covariates names not provided, using default naming convention:")
            print("- covariate1, ..., covariateC\n")
            self.smooth_opts['covariates_names'] = ["covariate%d" % (c) for c in range(self.dimensionalities["C"])]
        else:
            if isinstance(covariates_names,str):
                covariates_names = [covariates_names]
            assert isinstance(covariates_names, list), "covariates_names must be a string or a list"
            assert len(covariates_names)==self.dimensionalities["C"], "covariates_names must be of length equivalent to the number of covariates"
            self.smooth_opts['covariates_names'] = covariates_names

        self.sample_cov = sample_cov

        return None

    def set_data_matrix(self, data, likelihoods=None, views_names=None, groups_names=None, samples_names=None, features_names=None):
        """ Method to input the data in a wide matrix

        PARAMETERS
        ----------
        data: a nested list, first dimension for views, second dimension for groups.
              The dimensions of each matrix must be (samples,features)
        """

        if not hasattr(self, 'data_opts'): 
            # print("Data options not defined, using default values...\n")
            self.set_data_options()

        # Sanity check
        if not isinstance(data, list):
            if isinstance(data, dict):
                data = list(data.values())
            # if providing a single matrix, treat it as G=1 and M=1
            elif isinstance(data, pd.DataFrame):
                data = [[data.values]]
            elif isinstance(data, np.ndarray):
                data = [[data]]
            else:
                print("Error: Data not recognised"); sys.stdout.flush(); sys.exit()
        if len(data)==0:
            print("Error: Data is empty"); sys.stdout.flush(); sys.exit()

        # Convert input data to numpy array float64 format
        for m in range(len(data)):
            if isinstance(data[m], dict):
                data[m] = list(data[m].values())
            for p in range(len(data[m])):
                if not isinstance(data[m][p], np.ndarray):
                    if isinstance(data[m][p], pd.DataFrame):
                        data[m][p] = data[m][p].values
                    else:
                        print("Error, input data is not a numpy.ndarray or a pandas dataframe"); sys.stdout.flush(); sys.exit()
                data[m][p] = data[m][p].astype(np.float64)

        # Save dimensionalities
        M = self.dimensionalities["M"] = len(data)
        G = self.dimensionalities["G"] = len(data[0])
        N = self.dimensionalities["N"] = [data[0][p].shape[0] for p in range(len(data[0]))]
        D = self.dimensionalities["D"] = [data[m][0].shape[1] for m in range(len(data))]

        # Define views names
        if views_names is None:
            print("View names not provided, using default naming convention:")
            print("- view1, view2, ..., viewM\n")
            self.data_opts['views_names'] = [ "view" + str(m) for m in range(M) ]
        else:
            assert len(views_names)==self.dimensionalities["M"], "Length of views names is not the same as the number of views"
            self.data_opts['views_names']  = views_names

        # Define features names
        if features_names is None:
            print("Features names not provided, using default naming convention:")
            print("- feature1_view1, featureD_viewM\n")
            self.data_opts['features_names'] = [ ["feature%d_view%d" % (d,m) for d in range(D[m])] for m in range(M) ]
        else:
            assert len(features_names)==self.dimensionalities["M"], "views_names must a nested list with length equivalent to the number of views"
            self.data_opts['features_names'] = features_names

        # Define groups names
        if groups_names is None:
            print("Groups names not provided, using default naming convention:")
            print("- group1, group2, ..., groupG\n")
            self.data_opts['groups_names'] = [ "group" + str(g) for g in range(G) ]
        else:
            assert len(groups_names)==self.dimensionalities["G"], "Length of groups names is not the same as the number of groups"
            self.data_opts['groups_names'] = groups_names

        # Define samples names
        if samples_names is None:
            print("Samples names not provided, using default naming convention:")
            print("- sample1_group1, sample2_group1, sample1_group2, ..., sampleN_groupG\n")
            self.data_opts['samples_names'] = [ ["sample%d_group%d" % (n,g) for n in range(N[g])] for g in range(G) ]
        else:
            assert len(samples_names)==self.dimensionalities["G"], "samples_names must a nested list with length equivalent to the number of groups"
            self.data_opts['samples_names'] = samples_names

        # Check for duplicated entries
        assert len(self.data_opts['groups_names']) == len(set(self.data_opts['groups_names'])), "Duplicated groups names"
        assert len(self.data_opts['views_names']) == len(set(self.data_opts['views_names'])), "Duplicated views names"

        tmp = list(chain(*self.data_opts['samples_names']))
        assert len(tmp) == len(set(tmp)), "Duplicated entries found in samples_names. Make sure that samples names are not duplicated across different groups"
        
        tmp = list(chain(*self.data_opts['features_names']))
        assert len(tmp) == len(set(tmp)), "Duplicated entries found in features_names. Make sure that feature names are not duplicated across different views"

        # Set samples groups (list with dimensionality N where each row is the corresponding group name)
        # self.data_opts['samples_groups'] = [list(samples_names_dict.keys())[i] for i in range(len(self.data_opts['groups_names'])) for n in range(len(list(samples_names_dict.values())[i]))]

        # WHY THIS DOES NOT WORK?? asd = [ [self.data_opts['groups_names'][g]]*N[g] for g in range(G) ]
        self.data_opts['samples_groups'] = []
        for g in range(G):
            self.data_opts['samples_groups'].append( [self.data_opts['groups_names'][g]]*N[g] )
        self.data_opts['samples_groups'] = np.concatenate(self.data_opts['samples_groups'])

        # If everything successful, print verbose message
        for m in range(M):
            for g in range(G):
                print("Successfully loaded view='%s' group='%s' with N=%d samples and D=%d features..." % (self.data_opts['views_names'][m],self.data_opts['groups_names'][g], data[m][g].shape[0], data[m][g].shape[1]))
        print("\n")

        # Store intercepts
        self.intercepts = [ [ np.nanmean(data[m][g],axis=0) for g in range(G)] for m in range(M) ]

        # Concatenate groups in data
        for m in range(len(data)):
            data[m] = np.concatenate(data[m])
            if type(data[m]) != np.ndarray:
                data[m] = np.array(data[m])
        self.dimensionalities["N"] = np.sum(self.dimensionalities["N"])

        # Define likelihoods
        if likelihoods is None:
            likelihoods = guess_likelihoods(data)
        elif isinstance(likelihoods, str):
            likelihoods = [likelihoods]
        assert len(likelihoods)==self.dimensionalities["M"], "Please specify one likelihood for each view"
        assert set(likelihoods).issubset(set(["gaussian","bernoulli","poisson"])), "Available likelihoods are 'gaussian','bernoulli', 'poisson'"
        self.likelihoods = likelihoods

        # Process the data (center, scaling, etc.)
        self.data = process_data(data, likelihoods, self.data_opts, self.data_opts['samples_groups'])

    def set_train_options(self,
        iter=1000, startELBO=1, freqELBO=1, startSparsity=100, tolerance=None, convergence_mode="medium",
        startDrop=20, freqDrop=10, dropR2=None, nostop=False, verbose=False, quiet=False, seed=None,
        schedule=None, gpu_mode=False, Y_ELBO_TauTrick=True, save_parameters=False, weight_views = False,
        outfile=None, save_interrupted=False):
        """ Set training options """

        # Sanity checks
        assert hasattr(self, 'model_opts'), "Model options have to be defined before training options"

        self.train_opts = {}

        # Maximum number of iterations
        self.train_opts['maxiter'] = int(iter)

        # Lower bound computation frequency
        if freqELBO is None or freqELBO==0: freqELBO=iter+1
        self.train_opts['freqELBO'] = int(freqELBO)
        if startELBO==0: startELBO=1
        if startELBO==None: startELBO=iter+1
        self.train_opts['start_elbo'] = int(startELBO)

        # Verbosity
        self.train_opts['verbose'] = bool(verbose)
        self.train_opts['quiet'] = bool(quiet)

        # GPU mode
        if gpu_mode:
            try:
                import cupy as cp
                print("\nGPU mode is activated\n")
            except ImportError:
                print("\nGPU mode is activated, but GPU not found... switching to CPU mode")
                print('For GPU mode, you need:')
                print('1 - Make sure that you are running MOFA+ on a machine with an NVIDIA GPU')
                print('2 - Install CUPY following instructions on https://docs-cupy.chainer.org/en/stable/install.html\n')
                sys.stdout.flush()
                gpu_mode = False
        self.train_opts['gpu_mode'] = gpu_mode

        # Minimum Variance explained threshold to drop inactive factors
        if dropR2 is False: dropR2 = None
        if dropR2 is not None:
            dropR2 = float(dropR2)
            if dropR2 < 0: dropR2 = None
        self.train_opts['drop'] = { "min_r2":dropR2 }
        self.train_opts['start_drop'] = int(startDrop)
        self.train_opts['freq_drop'] = int(freqDrop)
        if ((dropR2 is not None) & (verbose is True)): print("\nDropping factors with minimum threshold of {0}% variance explained\n".format(dropR2))

        if ((dropR2 is not None) & (self.dimensionalities["N"]>1e4)):
            print("Warning: actively dropping factors during model training can be slow in large data sets.")
            print("Consider training the model with set drop_factor_threshold = -1 and prune them a posteriori")


        # Tolerance level for convergence
        if tolerance is not None:
            print("Warning: tolerance argument is depreciated, use the 'convergence_mode' argument instead")
            self.train_opts['tolerance'] = float(tolerance)

        # Convergence mode
        self.train_opts['convergence_mode'] = str(convergence_mode)
        if verbose: print("Convergence mode: %s\n" % convergence_mode)

        # Do no stop even when convergence criteria is met
        self.train_opts['forceiter'] = nostop

        # Iteration to activate spike and slab sparsity
        self.train_opts['start_sparsity'] = int(startSparsity)

        # By default, stochastic training is not activated
        self.train_opts['stochastic'] = False

        # Training schedule
        if schedule is None:
            schedule = ['Y', 'W', 'Z', 'Tau']

            # Insert ThetaW after W if Spike and Slab prior on W
            if self.model_opts['spikeslab_weights']:
                ix = schedule.index("W"); schedule.insert(ix+1, 'ThetaW')

            # Insert ThetaZ after Z if Spike and Slab prior on Z
            if self.model_opts['spikeslab_factors']:
                ix = schedule.index("Z"); schedule.insert(ix+1, 'ThetaZ')

            # Insert AlphaW after W if ARD prior on W
            if self.model_opts['ard_weights']:
                ix = schedule.index("W"); schedule.insert(ix+1, 'AlphaW')

            # Insert AlphaZ after Z if ARD prior on Z
            if self.model_opts['ard_factors']:
                ix = schedule.index("Z"); schedule.insert(ix+1, 'AlphaZ')

        else:
            assert set(["Y","W","Z","Tau"]) <= set(schedule)
            if self.model_opts['ard_factors']: assert "AlphaZ" in schedule
            if self.model_opts['ard_weights']: assert "AlphaW" in schedule
            if self.model_opts['spikeslab_factors']: assert "ThetaZ" in schedule
            if self.model_opts['spikeslab_weights']: assert "ThetaW" in schedule
            # if self.smooth_opts['GP_factors']: assert "Sigma" in schedule
            # if self.smooth_opts['sparseGP']: assert "U" in schedule

        self.train_opts['schedule'] = schedule

        # Seed
        if seed is None:  # Seed for the random number generator
            seed = int(round(time()*1000)%1e6)
        self.train_opts['seed'] = int(seed)
        # s.random.seed(self.train_opts['seed'])

        # Use TauTrick to speed up ELBO computation?
        self.train_opts['Y_ELBO_TauTrick'] = Y_ELBO_TauTrick

        # Save variational parameters?
        self.train_opts['save_parameters'] = save_parameters

        # Weight the views to avoid imbalance problems?
        self.train_opts['weight_views'] = weight_views

        # Output file name
        self.train_opts['outfile'] = outfile

        # If to save the partially trained model when the training is interrupted
        self.train_opts['save_interrupted'] = save_interrupted

    def set_smooth_options(self, 
        scale_cov = False,
        start_opt=20, n_grid=20, opt_freq=10, GP_factors = True, multivariate_Z = False,
        warping = False, warping_freq = 20, warping_ref = 0, warping_open_begin = True, warping_open_end = True,
        model_groups = False):

        # Sanity checks
        assert hasattr(self, 'smooth_opts'), "Please run set_covariates() before set_smooth_options()"
        assert hasattr(self, 'model_opts'), "Model options not defined. Please run set_model_opts() before set_smooth_options()"
        assert hasattr(self, 'train_opts'), "Training options not defined. Please run set_train_opts() before set_smooth_options()"
        assert self.sample_cov is not None, "Before setting smooth options, you need to define the covariates with set_covariates"

        self.smooth_opts['GP_factors'] = True

        # Define whether to scale covariates to unit variance
        self.smooth_opts['scale_cov'] = scale_cov
        # if (scale_cov): print("Scaling covariates to unit variance...\n")
        if self.smooth_opts['scale_cov']:
            self.sample_cov = (self.sample_cov - self.sample_cov.mean(axis=0)) / self.sample_cov.std(axis=0)

        # Define at which iteration to start optimizing the lengthscales, at which frequency and how many grid points
        start_opt = max(0, start_opt)
        self.smooth_opts['start_opt'] = int(start_opt)
        self.smooth_opts['n_grid'] = int(n_grid)
        self.smooth_opts['opt_freq'] = int(opt_freq)

        # inactivate group-wise ARD when using the SMOFA framework
        if self.model_opts['ard_factors'] is True:
            print("SMOFA framework is activated (GP_factors = True). This is not compatible with ARD prior on factors. Setting ard_factors to False...\n")
            self.model_opts['ard_factors'] = False

        # inactivate spike-slab on the factors when using the SMOFA framework
        if self.model_opts['spikeslab_factors'] is True:
            print("SMOFA framework is activated (GP_factors = True). This is not compatible with Spike-and-Slab prior on factors. Setting spikeslab_factors to False...\n")
            self.model_opts['spikeslab_factors'] = False


        # By default, no sparse GPs are used
        self.smooth_opts['sparseGP'] = False
        self.smooth_opts['idx_inducing'] = None


        # Define warping
        self.smooth_opts['warping'] = bool(warping)
        self.smooth_opts['warping_freq'] = int(warping_freq)
        self.smooth_opts['warping_ref'] = int(warping_ref)
        self.smooth_opts['warping_open_begin'] = bool(warping_open_begin)
        self.smooth_opts['warping_open_end'] = bool(warping_open_end)

        if self.smooth_opts['warping'] is True:
            if self.dimensionalities["C"] > 1:
                self.smooth_opts['warping'] = False
                print("Warping only implemented for one dimensional covariates, setting to False.")

        # Insert Sigma into training schedule if GP prior on Z
        if self.smooth_opts['GP_factors']:
            self.train_opts['schedule'].insert(len(self.train_opts['schedule']), "Sigma")

        # Define whether to model a group covariance structure
        self.smooth_opts['model_groups'] = model_groups
        self.smooth_opts['use_gpytorch'] = False # experimental, this could be passes as a model_option but to keep options uncluttered set to False

        # print("- Gaussian process prior on the factors: %s \n" % str(GP_factors))

    def set_stochastic_options(self, learning_rate=1., forgetting_rate=0., batch_size=1., start_stochastic=1):

        # Sanity checks
        if self.smooth_opts['GP_factors']:
            print("Stochastic inference is not possible when using covariates - train_opts['stochastic'] is set to False - consider using the option 'sparseGP' instead.")
            self.train_opts['stochastic'] = False
            return None
        assert hasattr(self, 'train_opts'), "Train options not defined"
        assert 0 < learning_rate <= 1, 'Learning rate must range from 0 and 1'
        # assert 0 < forgetting_rate <= 1, 'Forgetting rate must range from 0 and 1'
        assert 0 < batch_size <= 1, 'Batch size must range from 0 to 1'
        assert start_stochastic >= 1, 'start_stochastic must be >= 1'

        if self.train_opts['drop']["min_r2"] is not None:
            print("Dropping factors is currently disabled with stochastic inference...")
            self.train_opts['drop']["min_r2"] = None

        # Edit schedule: Z should come first (after Y) in the training schedule
        # (THIS IS DONE IN THE BAYESNET CLASS)
        # self.train_opts['schedule'].pop( self.train_opts['schedule'].index("Z") )
        # self.train_opts['schedule'].insert(1,"Z")

        self.train_opts['stochastic'] = True
        # self.train_opts['Y_ELBO_TauTrick'] = False # TauTrick speed up only works in non-stochastic mode
        self.train_opts['learning_rate'] = learning_rate
        self.train_opts['forgetting_rate'] = forgetting_rate
        self.train_opts['start_stochastic'] = start_stochastic
        self.train_opts['batch_size'] = batch_size

        self.train_opts['drop']["min_r2"] = None

    def set_sparseGP_options(self, n_inducing = None, idx_inducing = None, seed_inducing = None):
        """ Set options for sparse GPs (only when using smooth factors)"""

        # Sanity check
        assert hasattr(self, 'smooth_opts'), "Please run set_covariates() and set_smooth_options() before set_sparseGP_options()"
        assert len(self.smooth_opts.keys())>5, "Smooth options not defined. Please run set_smooth_options() before set_sparseGP_options()"
        assert hasattr(self, 'model_opts'), "Model options not defined. Please run set_model_options() before set_sparseGP_options()"
        assert hasattr(self, 'train_opts'), "Training options not defined. Please run set_train_options() before set_sparseGP_options()"
        assert self.sample_cov is not None, "Before setting sparse GP options, you need to define the covariates with set_covariates()"

        if not self.smooth_opts['GP_factors']:
            print("Sparse Gaussian Processes can only be used when having covariates and GP_factors set to True, setting to False")
            self.smooth_opts['sparseGP'] = False
            return None

        # Set the number of inducing points
        if n_inducing is None:
            n_inducing = max(0.2 * self.dimensionalities["N"], 100) # note: groups are already concatenated, N is total number of samples
        else:
            assert isinstance(n_inducing,int), "n_inducing has to be an integer"
            if n_inducing > self.dimensionalities["N"]:
                print("Number of inducing points is higher than original number of samples - using non-sparse GP inference")
                self.smooth_opts['sparseGP'] = False
                return None

        self.smooth_opts['sparseGP'] = True

        # Set the identity of the inducing points
        if idx_inducing is None:
            missing_sample_per_view = np.ones((self.dimensionalities["N"], self.dimensionalities["M"]))
            for m in range(len(self.data)):
                missing_sample_per_view[:,m] = np.isnan(self.data[m]).all(axis = 1)
            nonmissing_samples = np.where(missing_sample_per_view.sum(axis=1) != self.dimensionalities["M"])[0]
            N_nonmissing = len(nonmissing_samples)
            n_inducing = min(n_inducing, N_nonmissing)
            init_inducing_random = False # not used, could be passed as option
            if init_inducing_random:
                if not seed_inducing is None:
                    s.random.seed(int(seed_inducing))
                idx_inducing = np.random.choice(self.dimensionalities["N"], n_inducing, replace = False)
                idx_inducing.sort()
            else:
                N = self.dimensionalities["N"]
                loc = self.sample_cov.sum(axis = 1)
                groups = self.data_opts['samples_groups']
                nonmissing_samples_tiesshuffled = nonmissing_samples[np.lexsort((np.random.random(N_nonmissing), loc[nonmissing_samples]))] # shuffle ties randomly (e.g. between groups)
                grid_ix = np.floor(np.arange(0, N_nonmissing, step=N_nonmissing / n_inducing)).astype('int')
                if grid_ix[-1] == N_nonmissing: # avoid out of bound
                    grid_ix = grid_ix[:-1]
                idx_inducing = nonmissing_samples_tiesshuffled[grid_ix]

                # Show inducing points
                # import matplotlib.pyplot as plt
                # plt.figure(10)
                # color_labels = np.unique(self.data_opts['samples_groups'])
                # rgb_values = ['blue', 'red', 'green']
                # color_map = dict(zip(color_labels, rgb_values))
                # colors = [color_map[x] for x in self.data_opts['samples_groups']]
                # plt.scatter(loc, [x in idx_inducing for x in range(N)], c = colors)
                # plt.pause(1)

        # Insert U in schedule
        ix = self.train_opts['schedule'].index("Z")
        self.train_opts['schedule'].insert(ix, 'U')

        self.smooth_opts['n_inducing'] = n_inducing
        self.smooth_opts['idx_inducing'] = idx_inducing

    def set_model_options(self, factors=10, spikeslab_factors=False, spikeslab_weights=True, ard_factors=False, ard_weights=True):
        """ Set model options """

        self.model_opts = {}

        # Define whether to use sample-wise spike and slab prior for Z
        self.model_opts['spikeslab_factors'] = spikeslab_factors

        # Define whether to use feature-wise spike and slab prior for W
        self.model_opts['spikeslab_weights'] = spikeslab_weights

        # Define whether to use group and factor-wise ARD prior for Z
        # if ((self.dimensionalities["G"]>1) & (ard_factors==False)): 
        #     print("WARNING: 'ard_factors' should be set to True in model_options if using multiple groups\n")
        # always use the ARD prior on Z with more than one group except when using GPs
        if self.dimensionalities["G"]>1: ard_factors = True
        self.model_opts['ard_factors'] = ard_factors

        # Define whether to use view and factor-wise ARD prior for W
        # if ((self.dimensionalities["M"]>1) & (ard_weights==False)): 
        #     print("WARNING: 'ard_weights' should be set to True in model_options if using multiple views\n")
        if self.dimensionalities["M"]>1: ard_weights = True
        self.model_opts['ard_weights'] = ard_weights

        # Define initial number of latent factors
        self.dimensionalities["K"] = self.model_opts['factors'] = int(factors)


        # Define likelihoods
        self.model_opts['likelihoods'] = self.likelihoods

        print("Model options:")
        print("- Automatic Relevance Determination prior on the factors: %s" % str(ard_factors))
        print("- Automatic Relevance Determination prior on the weights: %s" % str(ard_weights))
        print("- Spike-and-slab prior on the factors: %s" % str(spikeslab_factors))
        print("- Spike-and-slab prior on the weights: %s" % str(spikeslab_weights))
        print("Likelihoods:")
        for m in range(self.dimensionalities["M"]):
          print("- View %d (%s): %s" % (m,self.data_opts["views_names"][m],self.likelihoods[m]) )
        print("\n")

    def set_data_options(self, scale_views=False, scale_groups = False, center_groups = True):
        """ Set data processing options """

        if not hasattr(self, 'data_opts'): self.data_opts = {}
        self.data_opts['center_groups'] = center_groups

        # Scale views to unit variance
        self.data_opts['scale_views'] = scale_views
        if (scale_views): print("Scaling views to unit variance...\n")
        
        # Scale groups to unit variance
        self.data_opts['scale_groups'] = scale_groups
        if (scale_groups): print("Scaling groups to unit variance...\n")

    def build(self):
        """ Build the model """

        # Sanity checks
        assert hasattr(self, 'train_opts'), "Training options not defined"
        assert hasattr(self, 'model_opts'), "Model options not defined"
        assert hasattr(self, 'dimensionalities'), "Dimensionalities are not defined"
        if hasattr(self, 'smooth_opts'):
            assert len(self.smooth_opts)>2, "Smooth covariates applied but smooth options not defined. Please define set_smooth_options() before build()"

        if np.any(np.array(self.dimensionalities["D"])<15):
            print("\nWarning: some view(s) have less than 15 features, MOFA won't be able to learn meaningful factors for these view(s)...\n")
        
        _, counts = np.unique(self.data_opts["samples_groups"], axis=0, return_counts=True)
        if np.any(counts<15):
            print("\nWarning: some group(s) have less than 15 samples, MOFA won't be able to learn meaningful factors for these group(s)...\n")

        # Build the nodes
        if hasattr(self, 'smooth_opts'):
            tmp = build_mofa_smooth(self.data, self.sample_cov, self.dimensionalities, self.data_opts, self.model_opts, self.train_opts, self.smooth_opts)
        else:
            tmp = buildBiofam(self.data, self.dimensionalities, self.data_opts, self.model_opts, self.train_opts)
        tmp.main()

        # Create BayesNet class
        if self.train_opts['stochastic']:
            self.model = StochasticBayesNet(self.dimensionalities, tmp.get_nodes())
        else:
            self.model = BayesNet(self.dimensionalities, tmp.get_nodes())

    @keyboardinterrupt_saver
    def run(self):
        """ Run the model """

        # Sanity checks
        assert hasattr(self, 'model_opts'), "Model options not defined"
        assert hasattr(self, 'train_opts'), "Train options not defined"
        assert hasattr(self, 'data_opts'), "Data options not defined"

        # Fetch training schedule (order of updates for the different nodes)
        if 'schedule' in self.train_opts:
            assert set(self.train_opts['schedule']) == set(list(self.model.getNodes().keys())), "Some nodes defined in the training schedule are not present in the model, or viceversa"
        else:
            self.train_opts['schedule'] = self.model_builder.schedule

        # Set training options
        self.model.setTrainOptions(self.train_opts)

        # Train the model
        train_model(self.model)

    def mask_outliers(self):

        Z = self.model.nodes['Z'].getExpectation()
        zscore_cutoff = 3 * 1.96   # z-score cutoff
        value_cutoff = 1       # max factor value

        for g in range(len(self.data_opts['groups_names'])):
            idx = np.where(np.array(self.data_opts["samples_groups"]) == self.data_opts['groups_names'][g])[0]
            Ztmp = Z[idx,:] # is this by reference? apparently not

            # calculate outlier score
            z_score = np.absolute((Ztmp-Ztmp.mean(axis=0))) / np.std(Ztmp, axis=0)

            # mask outliers with np.nan
            Ztmp[(z_score>zscore_cutoff) & (np.absolute(Z[idx,:])>value_cutoff)] = np.nan
            Z[idx,:] = Ztmp

    def predict_factor(self, new_covariates = None, uncertainty=True, groups = "all"):
        """
        Predict factor values at new covariate values or in missing groups
        new_covariates: new covariate values to predict at, should be of same format as covariates.
                        If None all present covariate values are used
        uncertainty: provide uncertainty measures for the prediction
        groups: Groups for which to predict, default is "all"
        """

        assert self.smooth_opts['GP_factors'], "Using predict_factors requires the use of GP_factors," \
                                              " maybe you want to use impute instead?"

        if new_covariates is None:
            new_covariates = self.sample_cov

        # get group-covariate combinations included in the model
        old_groups = self.model.nodes['Sigma'].groupsidx
        old_covariates = self.model.nodes['Sigma'].sample_cov_transformed
        all_covariates = np.unique(np.vstack([new_covariates, old_covariates]), axis =0)
        N = all_covariates.shape[0]
        M = new_covariates.shape[0]

        # get the necessary expectations
        Z = self.model.nodes['Z'].getExpectations()['E']
        K = Z.shape[1]
        Sigma_terms = self.model.nodes['Sigma'].getExpectations()
        Sigma = Sigma_terms['cov']
        Sigma_inv = Sigma_terms['inv']
        GP_param = self.model.nodes['Sigma'].getParameters()
        if groups == "all":
            G = len(self.model.nodes['Sigma'].groups)
            groups = np.arange(G)
        else:
            G = len(groups)
        if not self.model_opts['model_groups']:
            Kg = np.ones([K, G, G])
        else:
            Kg = self.model.nodes['Sigma'].getParameters()['Kg']
            Kg = Kg[:, groups, :][:, :, groups]

        # which rows/columns in Sigma_new correspond to original Sigma and which to new test covariates?
        old = np.hstack([old_groups[:,None], old_covariates])
        all = np.hstack([np.repeat(groups, all_covariates.shape[0])[:,None], np.tile(all_covariates.transpose(), G).transpose()])
        new = np.hstack([np.repeat(groups, new_covariates.shape[0])[:,None], np.tile(new_covariates.transpose(), G).transpose()])
        oldix = np.concatenate([np.where(np.all(np.equal(old[j,:], all), axis = 1))[0] for j in range(old.shape[0])])
        newidx = np.concatenate([np.where(np.all(np.equal(new[j,:], all), axis = 1))[0] for j in range(new.shape[0])])

        Z_new_mean = np.zeros([M*G, K])
        if uncertainty:
            Z_new_var = np.zeros([M*G, K])
        for k in range(K):
            Kc_new = self.model.nodes['Sigma'].Kc.eval_at_newpoints_k(all_covariates,k)
            K_new_k = GP_param['scale'][k] * np.kron(Kg[k,:,:], Kc_new)
            Sigma_new_k = K_new_k + (1 - GP_param['scale'][k]) * np.eye(N * G)
            Z_new_mean[:,k] = gpu_utils.dot(K_new_k[newidx, :][:, oldix], gpu_utils.dot(Sigma_inv[k,:,:], Z[:,k]))
            if uncertainty:
                # marginal variances p(z, z*|c, c*, y) =  p(z*|z, y,c,c*) *p(z|y,c,c*)  = p(z*|z, c*) * p(z|y,c)
                # Var[z*] = Var([E(z*|z)]) + Var(E(z*|z))
                # if not self.model_opts['mv_Znode']:
                #     Z2 = np.diag(self.model.nodes['Z'].getExpectations()['E2'][:,k] - self.model.nodes['Z'].getExpectations()['E'][:,k]**2)
                # else:
                Z2 = self.model.nodes['Z'].getExpectations()['cov'][k,:,:]
                # TODO: avoid calculating full matrix, only require diagonal
                Z_new_var[:, k] = np.diag(Sigma_new_k[newidx, :][:, newidx] - \
                                          gpu_utils.dot(K_new_k[newidx, :][:, oldix],
                                                        gpu_utils.dot(Sigma_inv[k, :, :],
                                                                      K_new_k[oldix, :][:, newidx])) + \
                                                            gpu_utils.dot(K_new_k[newidx, :][:, oldix], gpu_utils.dot(Sigma_inv[k, :, :],
                                                                           gpu_utils.dot(Z2,gpu_utils.dot(Sigma_inv[k, :, :],K_new_k[oldix,:][:,newidx])))))

        if uncertainty:
            assert np.all(Z_new_var >= 0), "Something went wrong in the prediction: variances of the predictive distribution are negative"
            self.Zpredictions = {"mean": Z_new_mean, "variance": Z_new_var, 'values' : new_covariates, 'groups' : groups}
        else:
            self.Zpredictions = {"mean": Z_new_mean, "variance": None, 'values' : new_covariates, 'groups' : groups}

        self.Zcompleted = True

    def impute(self, uncertainty=True, mask_outliers = True):
        """
        impute missing values with or without uncertainty estimates
        """

        # detect and mask outliers that could skew the results
        if mask_outliers:
            self.mask_outliers()

        # get the necessary expectations
        W = [w['E'] for w in self.model.nodes['W'].getExpectations()]
        Z = self.model.nodes['Z'].getExpectations()['E']

        # Predict the mean
        pred_mean = [Z.dot(w.T)for w in W]

        # for non-gaussian likelihoods, convert from pseudodata space to observation space
        for m in range(len(pred_mean)):
            if self.model_opts["likelihoods"][m]=="bernoulli":
                pred_mean[m] = np.round(np.exp(pred_mean[m])/(1+np.exp(pred_mean[m])))
            elif self.model_opts["likelihoods"][m]=="poisson":
                pred_mean[m] = np.round(np.log(1.+np.exp(pred_mean[m])))

        # Predict the variance
        if uncertainty:
            W2 = [w['E2'] for w in self.model.nodes['W'].getExpectations()]
            Z2 = self.model.nodes['Z'].getExpectations()['E2']
            Tau = [tau['E'] for tau in self.model.nodes['Tau'].getExpectations()]
            pred_var = [Z2.dot(W2[v].T) - (Z**2.).dot(W[v].T**2.) + 1./Tau[v] for v in range(len(W))]
            self.imputed_data = { "mean":pred_mean, "variance":pred_var }
        else:
            self.imputed_data = { "mean":pred_mean, "variance":None }

        # Only impute missing data
        for m in range(len(W)):
            mask = self.model.nodes['Y'].getNodes()[m].getMask()
            self.imputed_data["mean"][m][~mask] = self.data[m][~mask]
            self.imputed_data["variance"][m][~mask] = np.nan

        self.imputed = True # change flag

    def save(self, outfile=None, save_data=True, save_parameters=False, expectations=None):
        """ Save the model in an hdf5 file """

        # Sanity checks
        assert hasattr(self, 'data'), "Data has to be defined before training the model"
        assert hasattr(self, 'model'), "No trained model found"

        # Use outfile from training options if an outfile to override it is not provided
        if outfile is None or outfile == "":
            if self.train_opts["outfile"] is None or self.train_opts["outfile"] == "":
                outfile = os.path.join('/tmp', "mofa_{}.hdf5".format(strftime('%Y%m%d-%H%M%S')))
                print("No output file name provided as a training options or to the save method. Saving to {} .".format(
                    outfile
                    ))
            else:
                outfile = self.train_opts["outfile"]

        if os.path.isfile(outfile):
            print("Warning: Output file {} already exists, it will be replaced".format(outfile))

        # Create output directory
        if not os.path.isdir(os.path.dirname(outfile)) and (os.path.dirname(outfile) != ''):
            print("Output directory does not exist, creating it...")
            os.makedirs(os.path.dirname(outfile))
        print("Saving model in %s..." % outfile)
        sys.stdout.flush()

        # Save the model
        tmp = saveModel(
          model = self.model,
          outfile = outfile,
          data = self.data,
          intercepts = self.intercepts,
          samples_groups = self.data_opts['samples_groups'],
          sample_cov = self.sample_cov,
          train_opts = self.train_opts,
          model_opts = self.model_opts,
          samples_names = self.data_opts['samples_names'],
          features_names = self.data_opts['features_names'],
          views_names = self.data_opts['views_names'],
          groups_names = self.data_opts['groups_names'],
          covariates_names=self.smooth_opts['covariates_names'],
          samples_metadata = self.data_opts["samples_metadata"] if "samples_metadata" in self.data_opts else None,
          features_metadata = self.data_opts["features_metadata"] if "features_metadata" in self.data_opts else None,
          compression_level = 9
        )

        # Save sample and feature names
        tmp.saveNames()

        # Save metadata
        tmp.saveMetaData()

        # Save expectations

        # If all likelihoods are gaussian there is no need to save the expectations of Y, just saving the data is enough
        # TO-DO: THERE IS STH WRONG WITH THIS, CHECK WITH NON-GAUSS LIK
        # if all([i=="gaussian" for i in self.model_opts["likelihoods"]]):
        #     tmp.saveExpectations(nodes=["W","Z"])
        # else:
        #     tmp.saveExpectations(nodes=["W","Z"])

        if expectations is None:
            # Default is to save only W and Z nodes
            expectations = ["W", "Z", "Sigma"]
        
        tmp.saveExpectations(nodes=expectations)

        # Save parameters
        if save_parameters:
            tmp.saveParameters(nodes=["W","Z"])

        # Save model options
        tmp.saveModelOptions()

        # Save training options
        tmp.saveTrainOptions()

        # Save training statistics
        tmp.saveTrainingStats()

        # Save variance explained values
        tmp.saveVarianceExplained()

        # Save data
        if save_data: 
            tmp.saveData()

        # Save imputed data
        if self.imputed:
            tmp.saveImputedData(self.imputed_data["mean"], self.imputed_data["variance"])

        # Save predictions
        if self.Zcompleted:
            tmp.saveZpredictions(self.Zpredictions["mean"], self.Zpredictions["variance"], self.Zpredictions['values'], self.Zpredictions['groups'])


def mofa(adata, groups_label: bool = None, use_raw: bool = False, use_layer: bool = None, 
         features_subset: Optional[str] = None,
         likelihood: Optional[Union[str, List[str]]] = None, n_factors: int = 10,
         scale_views: bool = False, scale_groups: bool = False,
         ard_weights: bool = True, ard_factors: bool = True,
         spikeslab_weights: bool = True, spikeslab_factors: bool = False,
         n_iterations: int = 1000, convergence_mode: str = "fast",
         gpu_mode: bool = False, Y_ELBO_TauTrick: bool = True, 
         save_parameters: bool = False, save_data: bool = True, save_metadata: bool = True,
         seed: int = 1, outfile: Optional[str] = None,
         expectations: Optional[List[str]] = None,
         save_interrupted: bool = False,
         verbose: bool = False, quiet: bool = True, copy: bool = False):
    """
    Helper function to init and build the model in a single call
    from annotation data object

    PARAMETERS
    ----------
    adata: an AnnotationData object
    groups_label (optional): a column name in adata.obs for grouping the samples
    use_raw (optional): use raw slot of AnnData as input values
    use_layer (optional): use a specific layer of AnnData as input values (supersedes use_raw option)
    features_subset (optional): .var column with a boolean value to select genes (e.g. "highly_variable"), None by default
    likelihood (optional): likelihood to use, default is guessed from the data
    n_factors (optional): number of factors to train the model with
    scale_views (optional): scale views to unit variance
    scale_groups (optional): scale groups to unit variance
    ard_weights (optional): use view-wise sparsity
    ard_factors (optional): use group-wise sparsity
    spikeslab_weights (optional): use feature-wise sparsity (e.g. gene-wise)
    spikeslab_factors (optional): use sample-wise sparsity (e.g. cell-wise)
    n_iterations (optional): upper limit on the number of iterations
    convergence_mode (optional): fast, medium, or slow convergence mode
    gpu_mode (optional): if to use GPU mode
    Y_ELBO_TauTrick (optional): if to use ELBO Tau trick to speed up computations
    save_parameters (optional): if to save training parameters
    save_data (optional): if to save training data
    save_metadata (optional): if to load metadata from the AnnData object (.obs and .var tables) and save it, False by default
    seed (optional): random seed
    outfile (optional): path to HDF5 file to store the model
    expectations (optional): which nodes should be used to save expectations for (will save only W and Z by default);
    possible expectations names include Y, W, Z, Tau, AlphaZ, AlphaW, ThetaW, ThetaZ
    outfile (optional): output file name
    save_interrupted (optional): if to save partially trained model when the training is interrupted
    verbose (optional): print verbose information during traing
    quiet (optional): silence messages during training procedure
    copy (optional): return a copy of AnnData instead of writing to the provided object
    """

    ent = entry_point()

    lik = [likelihood] if likelihood is not None else None

    ent.set_data_options(scale_views=scale_views, scale_groups=scale_groups)
    ent.set_data_from_anndata(adata, groups_label=groups_label, use_raw=use_raw, use_layer=use_layer,
                              likelihoods=lik, features_subset=features_subset, save_metadata=save_metadata)
    ent.set_model_options(ard_factors=ard_factors, ard_weights=ard_weights, 
                          spikeslab_weights=spikeslab_weights, spikeslab_factors=spikeslab_factors, 
                          factors=n_factors)
    ent.set_train_options(iter=n_iterations, convergence_mode=convergence_mode, 
                          gpu_mode=gpu_mode, Y_ELBO_TauTrick=Y_ELBO_TauTrick,
                          seed=seed, verbose=verbose, quiet=quiet, outfile=outfile, save_interrupted=save_interrupted)

    ent.build()
    ent.run()

    ent.save(outfile, save_data=save_data, save_parameters=save_parameters, expectations=expectations)

    try:
        import h5py
    except ImportError:
        h5py = None

    if h5py:
        if outfile is not None:
            f = h5py.File(outfile)
            if copy:
                adata = adata.copy()
            adata.obsm['X_mofa'] = np.concatenate([v[:,:] for k, v in f['expectations']['Z'].items()], axis=1).T
            if features_subset is None:
                # Loadings can be saved only if all the features were used in training
                adata.varm['LFs'] = np.concatenate([v[:,:] for k, v in f['expectations']['W'].items()], axis=1).T
            if copy:
                return adata
            else:
                if features_subset is None:
                    print("Saved MOFA embeddings in adata.obsm['X_mofa'] slot and their loadings in adata.varm['LFs'].")
                else:
                    print("Saved MOFA embeddings in adata.obsm['X_mofa'] slot.")
    else:
        print("Can not add embeddings and loadings to AnnData object since h5py is not installed.")
