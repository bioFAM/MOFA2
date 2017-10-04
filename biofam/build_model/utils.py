from __future__ import division
from time import sleep

import numpy as np
import pandas as pd
import numpy.ma as ma
import os
import h5py


def removeIncompleteSamples(data):
    """ Method to remove samples with missing views

    PARAMETERS
    ----------
    data: list
    """
    print("Removing incomplete samples...")

    M = len(data)
    N = data[0].shape[0]
    samples_to_remove = []
    for n in range(N):
        for m in range(M):
            if pd.isnull(data[m].iloc[n][0]):
                samples_to_remove.append(n)
                break

    print("A total of " + str(len(samples_to_remove)) + " sample(s) have at least a missing view and will be removed")

    data_filt = [None]*M
    samples_to_keep = np.setdiff1d(range(N),samples_to_remove)
    for m in range(M):
        data_filt[m] = data[m].iloc[samples_to_keep]

    return data_filt

def maskData(data, data_opts):
    """ Method to mask values of the data,
    It is mainly to test missing values and to evaluate imputation

    PARAMETERS
    ----------
    data_opts: dic
    """
    print("Masking data with the following options:")
    print("at random:")
    print(data_opts['maskAtRandom'])
    print("full cases:")
    print(data_opts['maskNSamples'])

    for m in range(len(data)):

        # Mask values at random
        D = data[m].shape[1]
        N = data[m].shape[0]
        p2Mask = data_opts['maskAtRandom'][m]
        if p2Mask != 0:
            idxMask = np.zeros(N*D)
            idxMask[:int(round(N*D*p2Mask))] = 1
            np.random.shuffle(idxMask)
            idxMask = np.reshape(idxMask, [N, D])
            data[m] = data[m].mask(idxMask==1)

        # Mask samples in a complete view
        Nsamples2Mask = data_opts['maskNSamples'][m]
        if Nsamples2Mask != 0:
            idxMask = np.random.choice(N, size=Nsamples2Mask, replace = False)
            # idxMask = np.arange(Nsamples2Mask)
            # print idxMask
            tmp = data[m].copy()
            tmp.ix[idxMask,:] = pd.np.nan
            data[m] = tmp

    return data

def loadData(data_opts, verbose=True):
    """ Method to load the data

    PARAMETERS
    ----------
    data_opts: dic
    verbose: boolean
    """

    print ("\n")
    print ("#"*18)
    print ("## Loading data ##")
    print ("#"*18)
    print ("\n")
    sleep(1)

    M = len(data_opts['input_files'])

    Y =  [None]*M
    for m in range(M):

        # Read file
        file = data_opts['input_files'][m]
        Y[m] = pd.read_csv(file, delimiter=data_opts["delimiter"], header=data_opts["colnames"], index_col=data_opts["rownames"]).astype(pd.np.float32)

        # Y[m] = pd.read_csv(file, delimiter=data_opts["delimiter"])
        print("Loaded %s with dim (%d,%d)..." % (file, Y[m].shape[0], Y[m].shape[1]))

        # Removing features with no variance
        var = Y[m].std(axis=0)
        if np.any(var==0.):
            print("Warning: %d features(s) have zero variance, removing them..." % (var==0.).sum())
            Y[m].drop(Y[m].columns[np.where(var==0.)], axis=1, inplace=True)

        # Center the features
        if data_opts['center_features'][m]:
            print("Centering features for view " + str(m) + "...")
            Y[m] = (Y[m] - Y[m].mean(axis=0))

        # Scale the views to unit variance
        if data_opts['scale_views'][m]:
            print("Scaling view " + str(m) + " to unit variance...")
            Y[m] = Y[m] / np.nanstd(Y[m].as_matrix())

        # Scale the features to unit variance
        if data_opts['scale_features'][m]:
            print("Scaling features for view " + str(m) + " to unit variance...")
            Y[m] = Y[m] / np.std(Y[m], axis=0, )

        print("\n")

    return Y

def corr(A,B):
    """ Method to efficiently compute correlation coefficients between two matrices

    PARMETERS
    ---------
    A: np array
    B: np array
    """

    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:,None]
    B_mB = B - B.mean(1)[:,None]

    # Sum of squares across rows
    ssA = (A_mA**2).sum(1);
    ssB = (B_mB**2).sum(1);

    # Finally get corr coeff
    return np.dot(A_mA,B_mB.T)/np.sqrt(np.dot(ssA[:,None],ssB[None]))

def saveParameters(model, hdf5, view_names=None):
    """ Method to save the parameters of the model in an hdf5 file

    PARAMETERS
    ----------
    model: a BayesNet instance
    hdf5:
    view_names
    """

    # Get nodes from the model
    nodes = model.getNodes()

    # Create groups
    param_grp = hdf5.create_group("parameters")

    # Iterate over nodes
    for node in nodes:

        # Collect node parameters
        parameters = nodes[node].getParameters()

        # Create node subgroup
        node_subgrp = param_grp.create_group(node)

        # Multi-view nodes
        if type(parameters) == list:
            # Loop through the views
            for m in range(len(parameters)):
                if view_names is not None:
                    tmp = view_names[m]
                else:
                    tmp = "%d" % m
                # Create subsubgroup for the view
                view_subgrp = node_subgrp.create_group(tmp)
                # Loop through the parameters of the view
                if parameters[m] is not None:
                    # Variational nodes
                    if type(parameters[m]) == dict:
                        for param_name in parameters[m].keys():
                            if parameters[m][param_name] is not None:
                                view_subgrp.create_dataset(param_name, data=parameters[m][param_name].T)
                    # Non-variational nodes (no distributions)
                    elif type(parameters[m]) == np.ndarray:
                           view_subgrp.create_dataset("value", data=parameters[m].T)

        # Single-view nodes
        else:
            for param_name in parameters.keys():
                node_subgrp.create_dataset("%s" % (param_name), data=parameters[param_name].T)
    pass

def saveExpectations(model, hdf5, view_names=None, only_first_moments=True):
    """ Method to save the expectations of the model in an hdf5 file

    PARAMETERS
    ----------
    model: a BayesNet instance
    hdf5:
    view_names
    only_first_moments
    """
    # Get nodes from the model
    nodes = model.getNodes()

    exp_grp = hdf5.create_group("expectations")

    # Iterate over nodes
    for node in nodes:

        # Collect node expectations
        expectations = nodes[node].getExpectations()

        # Create subgroup for the node
        node_subgrp = exp_grp.create_group(node)

        # Multi-view nodes
        if type(expectations) == list:

            # Iterate over views
            for m in range(len(expectations)):
                if view_names is not None:
                    tmp = view_names[m]
                else:
                    tmp = "%d" % m

                # Create subsubgroup for the view
                view_subgrp = node_subgrp.create_group(tmp)

                # Loop through the expectations
                if only_first_moments:
                    if node == "SW":
                        expectations[m] = {'E':expectations[m]["E"], 'ES':expectations[m]["ES"], 'EW':expectations[m]["EW"]}
                    else:
                        expectations[m] = {'E':expectations[m]["E"]}

                if expectations[m] is not None:
                    for exp_name in expectations[m].keys():
                        if type(expectations[m][exp_name]) == ma.core.MaskedArray:
                            tmp = ma.filled(expectations[m][exp_name], fill_value=np.nan)
                            view_subgrp.create_dataset(exp_name, data=tmp.T)
                        else:
                            view_subgrp.create_dataset(exp_name, data=expectations[m][exp_name].T)

        # Single-view nodes
        else:
            if only_first_moments: expectations = {'E':expectations["E"]}
            for exp_name in expectations.keys():
                node_subgrp.create_dataset("%s" % (exp_name), data=expectations[exp_name].T)

def saveTrainingStats(model, hdf5):
    """ Method to save the training statistics in an hdf5 file

    PARAMETERS
    ----------
    model: a BayesNet instance
    hdf5:
    """
    stats = model.getTrainingStats()
    stats_grp = hdf5.create_group("training_stats")
    stats_grp.create_dataset("activeK", data=stats["activeK"])
    stats_grp.create_dataset("elbo", data=stats["elbo"])
    stats_grp.create_dataset("elbo_terms", data=stats["elbo_terms"].T)
    stats_grp['elbo_terms'].attrs['colnames'] = [a.encode('utf8') for a in stats["elbo_terms"].columns.values]

def saveTrainingOpts(opts, hdf5):
    """ Method to save the training options in an hdf5 file

    PARAMETERS
    ----------
    opts:
    hdf5:
    """
    # Remove dictionaries from the options
    for k,v in opts.copy().items():
        if type(v)==dict:
            for k1,v1 in v.items():
                opts[str(k)+"_"+str(k1)] = v1
            opts.pop(k)

    # Create HDF5 data set
    hdf5.create_dataset("training_opts", data=np.array(list(opts.values()), dtype=np.float))
    hdf5['training_opts'].attrs['names'] = np.asarray(list(opts.keys())).astype('S')

def saveModelOpts(opts, hdf5):
    """ Method to save the model options in an hdf5 file

    PARAMETERS
    ----------
    opts:
    hdf5:
    """
    opts_interest = ["learnIntercept","schedule","likelihood"]
    opts = dict((k, opts[k]) for k in opts_interest)
    grp = hdf5.create_group('model_opts')
    for k,v in opts.items():
        grp.create_dataset(k, data=np.asarray(v).astype('S'))
    grp[k].attrs['names'] = np.asarray(list(opts.keys())).astype('S')

def saveTrainingData(model, hdf5, view_names=None, sample_names=None, feature_names=None):
    """ Method to save the training data in an hdf5 file

    PARAMETERS
    ----------
    model: a BayesNet instance
    hdf5:
    view_names
    sample_names
    feature_names
    """
    data = model.getTrainingData()
    data_grp = hdf5.create_group("data")
    featuredata_grp = hdf5.create_group("features")
    hdf5.create_dataset("samples", data=sample_names)
    for m in range(len(data)):
        view = view_names[m] if view_names is not None else str(m)
        data_grp.create_dataset(view, data=data[m].data.T)
        if feature_names is not None:
            featuredata_grp.create_dataset(view, data=feature_names[m])

def saveModel(model, outfile, train_opts, model_opts, view_names=None, sample_names=None, feature_names=None):
    """ Method to save the model in an hdf5 file

    PARAMETERS
    ----------
    """
    assert model.trained == True, "Model is not trained yet"
    assert len(np.unique(view_names)) == len(view_names), 'View names must be unique'
    assert len(np.unique(sample_names)) == len(sample_names), 'Sample names must be unique'

    # For some reason h5py orders the datasets alphabetically, so we have to modify the likelihood accordingly
    idx = sorted(range(len(view_names)), key=lambda k: view_names[k])
    tmp = [model_opts["likelihood"][idx[m]] for m in range(len(model_opts["likelihood"]))]
    model_opts["likelihood"] = tmp

    hdf5 = h5py.File(outfile,'w')
    saveExpectations(model,hdf5,view_names)
    saveParameters(model,hdf5,view_names)
    saveTrainingStats(model,hdf5)
    saveTrainingOpts(train_opts,hdf5)
    saveModelOpts(model_opts,hdf5)
    saveTrainingData(model, hdf5, view_names, sample_names, feature_names)
    hdf5.close()
