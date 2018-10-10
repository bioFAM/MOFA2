class buildSpatialBiofam(buildBiofam):

    def __init__(self, data_opts, model_opts):
        assert 'data_x' in data_opts, 'positions not found in data options'
        super(buildSpatialBiofam, self).__init__(data_opts, model_opts)

    def build_all(self):
        # build the nodes which are missing in build_all
        self.build_Sigma()
        super(buildSimulationBiofam, self).build_all()

    def build_Sigma(self):
        if self.model_opts['cov_on'] == 'samples':
            # TODO add a if statement to check if there is a sigma_clust argument to see if blockSigma is needed
            if self.data_opts['dataClust'] is None:
                self.init_model.initSigmaZ(self.data_opts['data_x'])
            else:
                self.init_model.initSigmaBlockZ(self.data_opts['data_x'], clust=self.data_opts['dataClust'])
        else:
            params = [None] * M
            for m in range(M):
                if self.data_opts['view_has_covariance_prior'][m]:
                    params[m]={'X':self.data_opts['data_x'][m],
                    'sigma_clust':self.data_opts['dataClust'],'n_diag':0}
                else:
                    params[m]={'pa': 1e-14, 'pb':1e-14, 'qa':1., 'qb':1.}
            self.init_model.initMixedSigmaAlphaW(view_has_covariance_prior,params)

    def createMarkovBlankets(self):
        super(buildSpatialBiofam, self).createMarkovBlankets()
        nodes = self.init_model.getNodes()

        # create the markov blanket for the sigma nodes
        if self.model_opts['cov_on'] == 'samples':
            nodes["Sigma"].addMarkovBlanket(Z=nodes["Z"])
            nodes["Z"].addMarkovBlanket(Sigma=nodes["Sigma"])

        if self.model_opts['cov_on'] == 'features':
            nodes["Sigma"].addMarkovBlanket(W=nodes["W"])
            nodes["W"].addMarkovBlanket(Sigma=nodes["Sigma"])

    def createSchedule(self):
        super(buildSpatialBiofam, self).createSchedule()

        # add the sigma node at the right position in the schedule
        if self.model_opts['cov_on'] == 'samples':
            ix = self.find_node(self.schedule, 'Z')[0][0]
            np.insert(self.schedule, ix + 1, 'Sigma')

        if self.model_opts['cov_on'] == 'features':
            ix = self.find_node(self.schedule, 'W')[0][0]
            np.insert(self.schedule, ix + 1, 'Sigma')


class buildSimulationBiofam(buildBiofam):
    def __init__(self, model_opts):
        M = model_opts['M']
        N = model_opts['N']
        D = model_opts['D']

        self.init_model = initNewModel(dim, data, model_opts["likelihoods"])

    def build_all():
        # TODO add somewhere
        # if 'spatialFact' in model_opts:
        #     n_diag = (1-model_opts['spatialFact']) * K
        # else:
        #     n_diag = K
        # if model_opts["sample_X"]:
        #     if model_opts["transpose_sparsity"]:
        #        dataX = [s.random.normal(0, 1, [D[m], 2]) for m in range(M)]
        #        dataClust = [None] * M
        #        view_has_covariance_prior = [True] * M
        #     else:
        #        dataX = s.random.normal(0, 1, [N, 2])
        #        dataClust  = None
        # data = [np.ones([N, D[m]]) * np.nan for m in range(M)]
        # if model_opts["transpose_sparsity"]:
        #
        #    if dataX is None :
        #        view_has_covariance_prior = [None] * M
        #    else:
        #        view_has_covariance_prior = [dataX[m] is not None for m in range(M)]
        #
        #    #for test purpose
        #    if (dataX is None) and (model_opts["sample_X"]):
        #        dataX = [s.random.normal(0, 1, [D[m], 2]) for m in range(M)]
        #        dataClust = [None] * M
        #        view_has_covariance_prior = [True] * M
        #
        # #for test purpose
        # else:
        #    if (dataX is None) and (model_opts["sample_X"]):
        #        dataX = s.random.normal(0, 1, [N, 2])
        #        dataClust = None

        super(buildSimulationBiofam, self).build_all()
        self.build_Sigma()
        # build other stuff

    def build_Sigma():
        self.init_model.initSigma()
        self.init_model['Sigma'].addMarkovBlanket(self.init_model['Z'])

    def build_Tau():
        # for simulations we should change the parameters
        pass
    def build_ThetaZ():
        # TODO reimplement that so we can simulate different levels of sparsity
        if 'sparsity' in model_opts:
            assert 0. <= model_opts['sparsity'] <= 1., 'sparsty level must be between 0 and 1'
            priorTheta_a = 10.
            priorTheta_b = ((1 - model_opts['sparsity']) / model_opts['sparsity']) * priorTheta_a
        elif 'noise' in model_opts:
            priorTheta_a = 10.
            priorTheta_b = 10.

    def build_Tau(self):
        #TODO enable different levels of Noise    if 'noise' in model_opts:  # simulation
        pa = 20.; pb = pa * model_opts['noise']
        pa = 1e-14; pb = 1e-14

    def build_AlphaW(self):
        # also requires different prior hyperparameters
        pass
