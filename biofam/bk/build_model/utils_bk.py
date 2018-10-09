def loadDataX(data_opts, transpose = False):
    """ Method to load the data of the samples positions and assigned clusters
    """

    print ("\n")
    print ("#"*18)
    print ("## Loading samples positions data ##")
    print ("#"*18)
    print ("\n")

    M = len(data_opts["views_names"])

    if data_opts['x_files'] is None:

        X = None
        sigma_clust = None

    else:

        if transpose:

            assert M == len(data_opts['x_files']), "Length of view names and samples positions input files does not match"

            X = [None] * M
            sigma_clust = [None] * M
            for m in range(M):
                file = data_opts['x_files'][m]
                if file != "None":
                    try:
                        X[m] = np.loadtxt(file, delimiter=',')
                    except:
                        X[m] = np.loadtxt(file, delimiter=' ')

            if data_opts['permute_samples'] == 1:
                for m in range(M):
                    perm = np.random.permutation(D[m])
                    X[m] = X[m][perm, :]

            # load sigma cluster if among arguments
            if data_opts['sigma_clusters_file'] is not None:
                assert M == len(data_opts['sigma_clusters_file']), "Length of view names and samples clusters input files does not match"
                for m in range(M):
                    file = data_opts['sigma_clusters_file'][m]
                    if file != "None":
                        sigma_clust[m] = np.loadtxt(file)

            #if [np.all(X_m == None) for X_m in X] == [True]*M:
            if [X_m is None for X_m in X] == [True] * M: #simplified expression
                X = None

        else:

            assert 1 == len(data_opts['x_files']), "Length of view names and samples positions input files does not match"

            file = data_opts['x_files'][0]
            if file != "None":
                try:
                    X = np.loadtxt(file, delimiter=',')
                except:
                    X = np.loadtxt(file, delimiter=' ')

                if data_opts['permute_samples'] == 1:
                    perm = np.random.permutation(N)
                    X = X[perm, :]

                # load sigma cluster if among arguments
                if data_opts['sigma_clusters_file'] is not None:
                    assert 1 == len(data_opts['sigma_clusters_file']), "Length of view names and samples clusters input files does not match"
                    sigma_clust = np.loadtxt(data_opts['sigma_clusters_file'][0])
                else:
                    sigma_clust = None

            else:
                X = None
                sigma_clust = None

    return X, sigma_clust

def saveDataTxt(model, outDir, views_names=None, samples_names=None, features_names=None, shared_features=False):
    """ Method to save the training data in text files

    PARAMETERS
    ----------
    model: a BayesNet instance
    outDir
    views_names
    samples_names
    features_names
    """
    data = model.getTrainingData()
    for m in range(len(data)):
        view = views_names[m] if views_names is not None else 'view_' + str(m)
        file_name = outDir + '/' + view
        to_save = pd.DataFrame(data[m].data)
        if features_names is not None: #shared features or not -> to distinguish
            if shared_features:
                to_save.columns = features_names
            else:
                to_save.columns = features_names[m]
        if samples_names is not None:
            if shared_features:
                to_save.index = samples_names[m]
            else:
                to_save.index = samples_names
        to_save.to_csv(file_name)

def saveDataXTxt(model, outDir, transpose, views_names=None, samples_names=None):
    """ Method to save the X_Files data in text files (X_Files : positions of the samples)

    PARAMETERS
    ----------
    model: a BayesNet instance
    outDir
    views_names
    samples_names
    """
    if transpose:
        dataX = [SigmaAlphaW_m["X"] for SigmaAlphaW_m in model.getNodes()["SigmaAlphaW"].getParameters()]
        for m in range(len(dataX)):
            view = views_names[m] if views_names is not None else 'view_' + str(m)
            file_name = outDir + '/' + "X_file_" + view
            to_save = pd.DataFrame(dataX[m])
            #to_save.columns = ["x1", "x2"]
            if samples_names is not None:
                to_save.index = samples_names
            to_save.to_csv(file_name, index=False, header=False)
    else:
        dataX = model.getNodes()["SigmaZ"].getParameters()["X"]
        file_name = outDir + '/' + "X_file"
        to_save = pd.DataFrame(dataX)
        #to_save.columns = ["x1","x2"]
        if samples_names is not None:
            to_save.index = samples_names
        to_save.to_csv(file_name, index=False, header=False)

def overwriteExpectations(net):
    """
    methods to overwrite the expectations of the Q distributions with sampled
    values in cases where we don't train the network but do only simulations

    This enables saving the values more easily
    """
    for node in net.nodes.keys():
        if isinstance(net.nodes[node], Multiview_Node):
            overwriteExpectationsMV(net.nodes[node])
        if isinstance(net.nodes[node], Unobserved_Variational_Node):
            net.nodes[node].Q.expectations["E"] = net.nodes[node].samp
        if isinstance(net.nodes[node], Constant_Variational_Node):
            net.nodes[node].value = net.nodes[node].samp
        if node=='Sigma':
            net.nodes[node].ix = net.nodes[node].samp

def overwriteExpectationsMV(MV):
    """
    specific overwrite functions for multiview nodes
    """
    for node in MV.nodes:
        if isinstance(node, Unobserved_Variational_Node):
            node.Q.expectations["E"] = node.samp
        if isinstance(node, Constant_Variational_Node):
            node.value = node.samp
        if isinstance(node, Y_Node):
            node.value = node.samp
            node.mask()
        if isinstance(node, PseudoY):
            node.value = node.samp
            node.mask()

def saveSimulatedModel(model, outfile, train_opts, model_opts, views_names=None, groups_names=None, samples_names=None, features_names=None):
    """ Method to save the model in an hdf5 file

    PARAMETERS
    ----------
    """
    assert model.simulated, "Model is not simulated yet"

    if views_names is not None:
        assert len(np.unique(views_names)) == len(views_names), 'View names must be unique'

        # For some reason h5py orders the datasets alphabetically, so we have to modify the likelihood accordingly
        idx = sorted(range(len(views_names)), key=lambda k: views_names[k])
        tmp = [model_opts["likelihoods"][idx[m]] for m in range(len(model_opts["likelihoods"]))]
        model_opts["likelihoods"] = tmp

    # if features_names is not None:
    #         assert len(np.unique(features_names)) == len(features_names), 'Feature names must be unique'
    if samples_names is not None:
        assert len(np.unique(samples_names)) == len(samples_names), 'Sample names must be unique'

    overwriteExpectations(model)
    if 'outDir' in model_opts:
        saveDataTxt(model, model_opts['output_dir'], views_names, groups_names, samples_names, features_names)

        if model_opts['sample_x']:
            #saveDataXTxt(model, model_opts['outDir'], model_opts["transpose_sparsity"], views_names=views_names) #, samples_names=samples_names)
            saveDataXTxt(model, model_opts['output_dir'], views_names=views_names) #, samples_names=samples_names)

    hdf5 = h5py.File(outfile,'w')
    saveExpectations(model, hdf5, views_names, groups_names)
    saveParameters(model,hdf5, views_names, groups_names, samples_groups)
    saveModelOpts(model_opts, hdf5)
    saveTrainingData(model, hdf5, views_names, groups_names, samples_names, features_names)

    hdf5.close()
