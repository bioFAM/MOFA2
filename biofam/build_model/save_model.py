from __future__ import division
import numpy as np
import scipy as s
import pandas as pd
import numpy.ma as ma
import os
import h5py

from biofam.core.nodes import *


class saveModel():
    def __init__(self, model, outfile, data, samples_groups, train_opts, model_opts, features_names, views_names, samples_names):

        # Check that the model is trained
        assert model.trained, "Model is not trained"        
        self.model = model

        # Initialise hdf5 file
        self.hdf5 = h5py.File(outfile,'w')

        # Initialise training data
        self.data = data

        # Initialise samples groups
        assert len(samples_groups) == data[0].shape[0], "length of samples groups does not match the number of samples in the data"
        self.samples_groups = samples_groups

        # Initialise options
        self.train_opts = train_opts
        self.model_opts = model_opts

        # Initialise dimension names
        self.views_names = views_names
        self.samples_names = samples_names
        self.features_names = features_names
        self.groups_names = set(samples_groups)

    def saveData(self):
        """ Method to save the training data"""

        # Create HDF5 groups
        data_grp = self.hdf5.create_group("data")
        features_grp = self.hdf5.create_group("features")
        samples_grp  = self.hdf5.create_group("samples")

        # Save samples names
        for g in self.groups_names:
            samples_idx = np.where(np.array(self.samples_groups) == g)[0]
            samples_names = [s for s in [self.samples_names[e] for e in samples_idx]]
            # samples_grp.create_dataset(g, data=[str(s).encode('utf8') for s in [self.samples_names[e] for e in samples_idx]])
            samples_grp.create_dataset(g, data=np.array(samples_names, dtype='S50'))

        # Save feature names
        for m in range(len(self.data)):
            # features_grp.create_dataset(self.views_names[m], data=[str(x).encode('utf8') for x in self.features_names[m]])
            features_grp.create_dataset(self.views_names[m], data=np.array(self.features_names[m], dtype='S50'))

        # Save data
        for m in range(len(self.data)):
            view_subgrp = data_grp.create_group(self.views_names[m])
            for g in self.groups_names:
                samples_idx = np.where(np.array(self.samples_groups) == g)[0]
                view_subgrp.create_dataset(g, data=self.data[m][samples_idx,:])

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

                    # Multi-groups nodes
                    if n in multigroup_nodes:
                        for g in self.groups_names:

                            # Create subgroup for the view
                            view_subgrp = node_subgrp.create_group(self.views_names[m])

                            # create hdf5 data set for the expectation
                            samp_indices = np.where(np.array(self.samples_groups) == g)[0]
                            view_subgrp.create_dataset(g, data=exp[m][samp_indices,:])

                    # Single-groups nodes
                    else:
                        node_subgrp.create_dataset(self.views_names[m], data=exp[m].T)

            # Single-view nodes
            else:

                # Multi-group nodes
                if n in multigroup_nodes:
                    for g in self.groups_names:
                        samp_indices = np.where(np.array(self.samples_groups) == g)[0]
                        node_subgrp.create_dataset(g, data=exp[samp_indices,:].T)

                # Single-group nodes
                else:
                    node_subgrp.create_dataset("E", data=exp.T)

        pass

    def saveParameters(self, nodes="all"):
        print("saveParameters is depreciated")
        pass

    def saveModelOptions(self):
        # TO-DO: 
        # (1) RENAME SL_Z and SL_W to SPARsity_z AND SPARSITY_w

        # Subset model options
        opts = dict((k, self.model_opts[k]) for k in ["likelihoods", "sl_z", "sl_w"])

        # Create HDF5 group
        grp = self.hdf5.create_group('model_options')

        # Create HDF5 data sets
        for k, v in opts.items():
            grp.create_dataset(k, data=np.asarray(v).astype('S'))
        grp[k].attrs['names'] = np.asarray(list(opts.keys())).astype('S')

    def saveTrainOptions(self):
        """ Method to save the training options """

        # TO-DO:
        # Currently only numeric options can be saved due to compatibility problems between hdf5 and strings
        # For more information see: https://github.com/h5py/h5py/pull/1032 or https://github.com/h5py/h5py/issues/289


        # Replace dictionaries (not supported in hdf5) by lists 
        opts = self.train_opts
        for k,v in opts.copy().items():
            if type(v)==dict:
                for k1,v1 in v.items():
                    opts[str(k)+"_"+str(k1)] = v1
                opts.pop(k)

        # Remove schedule from training options
        # self.train_opts['schedule'] = '_'.join(self.train_opts['schedule'])
        if 'schedule' in opts.keys():
            del opts['schedule']

        # Create data set: only numeric options 
        self.hdf5.create_dataset("training_opts".encode('utf8'), data=np.array(list(opts.values()), dtype=np.float))
        self.hdf5['training_opts'].attrs['names'] = np.asarray(list(opts.keys())).astype('S')

    def saveTrainingStats(self):
        """ Method to save the training statistics """

        # Get training statistics
        # stats = self.model.getTrainingStats()

        # Create HDF5 group
        # stats_grp = self.hdf5.create_group("training_stats")

        # stats_grp.create_dataset("activeK", data=stats["activeK"])
        # stats_grp.create_dataset("elbo", data=stats["elbo"])
        # stats_grp.create_dataset("elbo_terms", data=stats["elbo_terms"].T)
        # stats_grp['elbo_terms'].attrs['colnames'] = [a.encode('utf8') for a in stats["elbo_terms"].columns.values]

        print("Depreciated")
        pass
