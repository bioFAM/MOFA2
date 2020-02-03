from __future__ import division
import numpy as np
import scipy as s
import pandas as pd
import numpy.ma as ma
import os
import h5py
from mofapy2.core.nodes import *
from mofapy2.core.nodes import *

# To keep same order of views and groups in the hdf5 file
# h5py.get_config().track_order = True

class saveModel():
    def __init__(self, model, outfile, data, intercepts, samples_groups, train_opts, model_opts, features_names, views_names, samples_names, groups_names, compression_level=9):

        # Check that the model is trained
        assert model.trained, "Model is not trained"        
        self.model = model

        # Initialise hdf5 file
        self.hdf5 = h5py.File(outfile,'w')
        self.compression_level = compression_level

        # Initialise training data
        self.data = data

        # Define masks
        self.mask = [ x.mask for x in model.getNodes()["Y"].getNodes() ]

        # Initialise samples groups
        assert len(samples_groups) == data[0].shape[0], "length of samples groups does not match the number of samples in the data"
        self.samples_groups = samples_groups

        # Initialise intercepts
        self.intercepts = intercepts

        # Initialise options
        self.train_opts = train_opts
        self.model_opts = model_opts

        # Initialise dimension names
        self.views_names = views_names
        self.samples_names = samples_names
        self.features_names = features_names
        self.groups_names = groups_names


    def saveNames(self):
        """ Method to save sample and feature names"""

        # Save group names
        groups_grp = self.hdf5.create_group("groups")
        groups_grp.create_dataset("groups", data=np.array(self.groups_names, dtype='S50'))        

        # Save views names
        views_grp = self.hdf5.create_group("views")
        views_grp.create_dataset("views", data=np.array(self.views_names, dtype='S50'))

        # Save samples names
        samples_grp = self.hdf5.create_group("samples")
        for g in range(len(self.groups_names)):
            samples_grp.create_dataset(self.groups_names[g], data=np.array(self.samples_names[g], dtype='S50'))

        # Save feature names
        features_grp = self.hdf5.create_group("features")
        for m in range(len(self.data)):
            features_grp.create_dataset(self.views_names[m], data=np.array(self.features_names[m], dtype='S50'))

    def saveData(self):
        """ Method to save the training data"""
        
        # Create HDF5 groups
        data_grp = self.hdf5.create_group("data")
        intercept_grp = self.hdf5.create_group("intercepts")

        for m in range(len(self.data)):
            data_subgrp = data_grp.create_group(self.views_names[m])
            intercept_subgrp = intercept_grp.create_group(self.views_names[m])
            for g in range(len(self.groups_names)):

                # Subset group
                samples_idx = np.where(np.array(self.samples_groups) == self.groups_names[g])[0]
                tmp = self.data[m][samples_idx,:]

                # Mask missing values
                tmp[self.mask[m][samples_idx,:]] = np.nan
                
                # Create hdf5 data set for data
                data_subgrp.create_dataset(self.groups_names[g], data=tmp, compression="gzip", compression_opts=self.compression_level)
                
                # Create hdf5 data set for intercepts
                intercept_subgrp.create_dataset(self.groups_names[g], data=self.intercepts[m][g])

    def saveImputedData(self, mean, variance):
        """ Method to save the training data"""
        
        # Create HDF5 groups
        data_grp = self.hdf5.create_group("imputed_data")

        # Save mean
        for m in range(len(mean)):
            view_subgrp = data_grp.create_group(self.views_names[m])
            for g in range(len(self.groups_names)):

                # Subset group
                samples_idx = np.where(np.array(self.samples_groups) == self.groups_names[g])[0]

                # Create HDF5 subgroup
                group_subgrp = view_subgrp.create_group(self.groups_names[g])

                # Create hdf5 data sets for the mean and the variance
                group_subgrp.create_dataset("mean", data=mean[m][samples_idx,:], compression="gzip", compression_opts=self.compression_level)
                if variance is not None:
                    group_subgrp.create_dataset("variance", data=variance[m][samples_idx,:], compression="gzip", compression_opts=self.compression_level)
                
    def saveExpectations(self, nodes="all"):

        # Get nodes from the model
        nodes_dic = self.model.getNodes()
        if type(nodes) is str:
            nodes = list(nodes_dic.keys()) if nodes=="all" else [nodes]
        elif type(nodes) is list or type(nodes) is tuple:
            assert set(nodes).issubset(["Z","W","Y","Tau","AlphaW","AlphaZ","ThetaZ","ThetaW"]), "Unrecognised nodes"
        nodes_dic = {x: nodes_dic[x] for x in nodes if x in nodes_dic}

        # Define nodes which special characteristics 
        # (note that this is ugly and is not proper class-oriented programming)
        multigroup_nodes = ["Y","Tau","Z"]

        # Create HDF5 group
        grp = self.hdf5.create_group("expectations")

        # Iterate over nodes
        for n in nodes_dic:
            # Create subgroup for the node
            node_subgrp = grp.create_group(n)

            # Collect node expectation
            exp = nodes_dic[n].getExpectation()

            # Multi-view nodes
            if isinstance(nodes_dic[n],Multiview_Node):
                for m in range(nodes_dic[n].M):

                    # Multi-groups nodes (Tau and Y)
                    if n in multigroup_nodes:

                        # Create subgroup for the view
                        view_subgrp = node_subgrp.create_group(self.views_names[m])
                        
                        for g in self.groups_names:

                            # Add missing values to Tau and Y nodes
                            exp[m][self.mask[m]] = np.nan

                            # create hdf5 data set for the expectation
                            samp_indices = np.where(np.array(self.samples_groups) == g)[0]

                            view_subgrp.create_dataset(g, data=exp[m][samp_indices,:], compression="gzip", compression_opts=self.compression_level)

                    # Single-groups nodes (W)
                    else:
                        node_subgrp.create_dataset(self.views_names[m], data=exp[m].T, compression="gzip", compression_opts=self.compression_level)

            # Single-view nodes
            else:

                # Multi-group nodes
                if n in multigroup_nodes:
                    for g in self.groups_names:
                        samp_indices = np.where(np.array(self.samples_groups) == g)[0]
                        node_subgrp.create_dataset(g, data=exp[samp_indices,:].T, compression="gzip", compression_opts=self.compression_level)

                # Single-group nodes
                else:
                    node_subgrp.create_dataset("E", data=exp.T, compression="gzip", compression_opts=self.compression_level)

        pass

    def saveParameters(self, nodes="all"):

        # Get nodes from the model
        nodes_dic = self.model.getNodes()
        if type(nodes) is str:
            nodes = list(nodes_dic.keys()) if nodes=="all" else [nodes]
        elif type(nodes) is list or type(nodes) is tuple:
            assert set(nodes).issubset(["Z","W","Tau","AlphaW","AlphaZ","ThetaZ","ThetaW"]), "Unrecognised nodes"
        nodes_dic = {x: nodes_dic[x] for x in nodes if x in nodes_dic}

        # Define nodes which special characteristics 
        # (note that this is ugly and is not proper class-oriented programming)
        multigroup_nodes = ["Y","Tau","Z"]

        # Create HDF5 group
        grp = self.hdf5.create_group("parameters")

        # Iterate over nodes
        for n in nodes_dic:
            # Create subgroup for the node
            node_subgrp = grp.create_group(n)

            # Collect node parameters
            par = nodes_dic[n].getParameters()

            # Multi-view nodes
            if isinstance(nodes_dic[n],Multiview_Node):
                for m in range(nodes_dic[n].M):

                    # Create subgroup for the view
                    view_subgrp = node_subgrp.create_group(self.views_names[m])

                    # Multi-groups nodes
                    if n in multigroup_nodes:

                        for g in self.groups_names:
                            grp_subgrp = view_subgrp.create_group(g)

                            # create hdf5 data set for the parameter
                            samp_indices = np.where(np.array(self.samples_groups) == g)[0]

                            for k in par[m].keys():
                                tmp = par[m][k][samp_indices,:]
                                grp_subgrp.create_dataset(k, data=tmp, compression="gzip", compression_opts=self.compression_level)

                    # Single-groups nodes
                    else:
                        for k in par[m].keys():
                            if k not in ["mean_B0","var_B0"]:
                                tmp = par[m][k].T
                                view_subgrp.create_dataset(k, data=tmp, compression="gzip", compression_opts=self.compression_level)

            # Single-view nodes
            else:

                # Multi-group nodes
                if n in multigroup_nodes:
                    for g in self.groups_names:
                        grp_subgrp = node_subgrp.create_group(g)
                        samp_indices = np.where(np.array(self.samples_groups) == g)[0]

                        for k in par.keys():
                            tmp = par[k][samp_indices,:].T
                            grp_subgrp.create_dataset(k, data=tmp, compression="gzip", compression_opts=self.compression_level)

                # Single-group nodes
                else:
                    for k in par.keys():
                        node_subgrp.create_dataset(k, data=par[k].T, compression="gzip", compression_opts=self.compression_level)

        pass

    def saveModelOptions(self):

        # Subset model options
        options_to_save = ["likelihoods", "spikeslab_factors", "spikeslab_weights", "ard_factors", "ard_weights"]
        opts = dict((k, np.asarray(self.model_opts[k]).astype('S')) for k in options_to_save)

        # Sort values by alphabetical order of views
        # order = np.argsort(self.views_names)
        # opts["likelihoods"] = opts["likelihoods"][order]

        # Create HDF5 group
        grp = self.hdf5.create_group('model_options')

        # Create HDF5 data sets
        for k, v in opts.items():
            grp.create_dataset(k, data=v)
        grp[k].attrs['names'] = np.asarray(list(opts.keys())).astype('S')

    def saveTrainOptions(self):
        """ Method to save the training options """

        # TO-DO:
        # Currently only numeric options can be saved due to compatibility problems between hdf5 and strings
        # For more information see: https://github.com/h5py/h5py/pull/1032 or https://github.com/h5py/h5py/issues/289

        # Subset training options
        opts = dict((k, self.train_opts[k]) for k in ["maxiter", "freqELBO", "start_elbo", "gpu_mode", "stochastic", "seed"])

        # Replace dictionaries (not supported in hdf5) by lists 
        # opts = self.train_opts
        for k,v in opts.copy().items():
            if type(v)==dict:
                for k1,v1 in v.items():
                    opts[str(k)+"_"+str(k1)] = v1
                opts.pop(k)

        # Remove strings from training options
        # self.train_opts['schedule'] = '_'.join(self.train_opts['schedule'])
        # if 'schedule' in opts.keys():
        #     del opts['schedule']
        # if 'convergence_mode' in opts.keys():
        #     del opts['convergence_mode']

        # Remove some training options
        # del opts['quiet']; del opts['start_drop']; del opts['freq_drop']; del opts['forceiter']; del opts['start_sparsity']; del opts['Y_ELBO_TauTrick']

        # Create data set: only numeric options 
        self.hdf5.create_dataset("training_opts".encode('utf8'), data=np.array(list(opts.values()), dtype=np.float))
        self.hdf5['training_opts'].attrs['names'] = np.asarray(list(opts.keys())).astype('S')

    def saveVarianceExplained(self):

        # Sort values by alphabetical order of views
        # order = np.argsort(self.views_names)
        # # order = [ i[0] for i in sorted(enumerate(self.views_names), key=lambda x:x[1]) ]

        # Store variance explained per factor in each view and group
        grp = self.hdf5.create_group("variance_explained")

        subgrp = grp.create_group("r2_per_factor")
        r2 = self.model.calculate_variance_explained()
        for g in range(len(self.groups_names)):
            # subgrp.create_dataset(self.groups_names[g], data=r2[g][order], compression="gzip",
            subgrp.create_dataset(self.groups_names[g], data=r2[g], compression="gzip",
                               compression_opts=self.compression_level)

        # Store total variance explained for each view and group (using all factors)
        subgrp = grp.create_group("r2_total")
        r2 = self.model.calculate_variance_explained(total=True)
        for g in range(len(self.groups_names)):
            # subgrp.create_dataset(self.groups_names[g], data=r2[g][order], compression="gzip",
            subgrp.create_dataset(self.groups_names[g], data=r2[g], compression="gzip",
                               compression_opts=self.compression_level)

    def saveTrainingStats(self):
        """ Method to save the training statistics """

        # Get training statistics
        stats = self.model.getTrainingStats()

        # Create HDF5 group
        stats_grp = self.hdf5.create_group("training_stats")

        stats_grp.create_dataset("number_factors", data=stats["number_factors"])
        stats_grp.create_dataset("time", data=stats["time"])
        stats_grp.create_dataset("elbo", data=stats["elbo"])
        # stats_grp.create_dataset("elbo_terms", data=stats["elbo_terms"].T)
        # stats_grp['elbo_terms'].attrs['colnames'] = [a.encode('utf8') for a in stats["elbo_terms"].columns.values]
