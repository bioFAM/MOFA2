
  # TODO have some config file for that as an input maybe ?
  def define_priors(self):
    """ Define priors of the model"""

    N = self.dimensionalities["N"]
    K = self.dimensionalities["K"]
    M = self.dimensionalities["M"]
    D = self.dimensionalities["D"]

    # Latent Variables
    self.model_opts["priorZ"] = { 'mean':s.zeros((N,K)) }
    self.model_opts["priorZ"]['var'] = s.ones((K,))*1.

    # Weights
    self.model_opts["priorSW"] = { 'Theta':[s.nan]*M, 'mean_S0':[s.nan]*M, 'var_S0':[s.nan]*M, 'mean_S1':[s.nan]*M, 'var_S1':[s.nan]*M } # Not required
    self.model_opts["priorAlphaW"] = { 'a':[s.ones(K)*1e-14]*M, 'b':[s.ones(K)*1e-14]*M }

    # Theta
    self.model_opts["priorTheta"] = { 'a':[s.ones(K,) for m in range(M)], 'b':[s.ones(K,) for m in range(M)] }
    for m in range(M):
      for k in range(K):
        if self.model_opts['learnTheta'][m][k]==0:
          self.model_opts["priorTheta"]["a"][m][k] = s.nan
          self.model_opts["priorTheta"]["b"][m][k] = s.nan

    # Tau
    self.model_opts["priorTau"] = { 'a':[s.ones(D[m])*1e-14 for m in range(M)], 'b':[s.ones(D[m])*1e-14 for m in range(M)] }

  def initialise_variational(self, initTheta=1.):
    """ Initialise variational distributions of the model"""

    N = self.dimensionalities["N"]
    K = self.dimensionalities["K"]
    M = self.dimensionalities["M"]
    D = self.dimensionalities["D"]

    # Latent variables
    self.model_opts["initZ"] = { 'mean':"random", 'var':s.ones((K,)), 'E':None, 'E2':None }

    # Tau
    self.model_opts["initTau"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[s.ones(D[m])*100 for m in range(M)] }

    # ARD of weights
    self.model_opts["initAlphaW"] = { 'a':[s.nan]*M, 'b':[s.nan]*M, 'E':[s.ones(K)*1. for m in range(M)] }

    # Theta
    self.model_opts["initTheta"] = { 'a':[s.ones(K,) for m in range(M)], 'b':[s.ones(K,) for m in range(M)], 'E':[s.nan*s.zeros((D[m],K)) for m in range(M)] }
    if type(initTheta) is float:
      self.model_opts['initTheta']['E'] = [s.ones((D[m],K))*initTheta for m in range(M)]
    # elif type(self.model_opts["initTheta"]) == list:
    #   assert len(self.model_opts["initTheta"]) == M, "--initTheta has to be a binary vector with length number of views"
    #   self.model_opts['initTheta']['E']= [ self.model_opts["initTheta"][m]*s.ones((D[m],K)) for m in range(M) ]
    else:
       print("Error: 'initTheta' must be a float")
       exit()

    for m in range(M):
      for k in range(K):
        if self.model_opts['learnTheta'][m][k]==0.:
          self.model_opts["initTheta"]["a"][m][k] = s.nan
          self.model_opts["initTheta"]["b"][m][k] = s.nan

    # Weights
    self.model_opts["initSW"] = {
      'Theta':[ self.model_opts['initTheta']['E'][m] for m in range(M)],
      'mean_S0':[s.zeros((D[m],K)) for m in range(M)],
      'var_S0':[s.nan*s.ones((D[m],K)) for m in range(M)],
      'mean_S1':[s.zeros((D[m],K)) for m in range(M)],
      # 'mean_S1':[stats.norm.rvs(loc=0, scale=1, size=(D[m],K)) for m in range(M)],
      'var_S1':[s.ones((D[m],K)) for m in range(M)],
      'ES':[None]*M, 'EW_S0':[None]*M, 'EW_S1':[None]*M # It will be calculated from the parameters
    }

        # Load covariates (NOT IMPLEMENTED)
        # if self.data_opts['covariatesFile'] is not None:

        #   self.data_opts['covariates'] = pd.read_csv(self.data_opts['covariatesFile'], delimiter=" ", header=None).as_matrix()
        #   print("Loaded covariates from " + self.data_opts['covariatesFile'] + "with shape " + str(self.data_opts['covariates'].shape) + "...")

        #   # Scale covariates
        #   self.data_opts['scale_covariates'] = self.data_opts['scale_covariates']
        #   if len(self.data_opts['scale_covariates']) == 1 and self.data_opts['covariates'].shape[1] > 1:
        #     self.data_opts['scale_covariates'] = self.data_opts['scalecovariates'][0] * s.ones(self.data_opts['covariates'].shape[1])
        #   elif type(self.data_opts['scale_covariates'])==list:
        #     assert len(self.data_opts['scale_covariates']) == self.data_opts['covariates'].shape[1], "'scale_covariates' has to be the same length as the number of covariates"
        #   self.data_opts['scale_covariates'] = [ bool(x) for x in self.data_opts['scale_covariates'] ]

        #   # Add the number of covariates to the total number of factors
        #   self.model_opts['factors'] += self.data_opts['covariates'].shape[1]

        #   # Parse covariates
        #   self.parse_covariates()

        # else:
        #   self.data_opts['scale_covariates'] = False
        #   self.data_opts['covariates'] = None

  def parse_intercept(self):
    """ Parse intercept factor """

    K = self.dimensionalities["K"]
    M = self.dimensionalities["M"]
    N = self.dimensionalities["N"]

    # If we want to learn the intercept, we add a constant covariate of 1s
    if self.model_opts['learn_intercept']:
      if self.data_opts['covariates'] is not None:
        self.data_opts['covariates'] = s.insert(self.data_opts['covariates'], obj=0, values=1, axis=1)
        self.data_opts['scale_covariates'].insert(0,False)
      else:
        self.data_opts['covariates'] = s.ones((N,1))
        self.data_opts['scale_covariates'] = [False]

      # Parse intercept
      # self.model_opts['factors'] += 1
      # self.dimensionalities["K"] += 1

      # Remove sparsity from the Intercept factor
      # TO-DO: CHECK THAT THE MODEL IS ALREADY NOT SPARSE
      # TO-DO: RECHECK THIS, ITS UGLY
      # stop if not self.model_opts["learn_intercept"] == TRUE

      for m in range(M):

        # Weights
        if self.model_opts["likelihoods"][m]=="gaussian":
          self.model_opts["initSW"]["mean_S1"][m][:,0] = self.data[m].mean(axis=0)
          self.model_opts["initSW"]["var_S1"][m][:,0] = 1e-5

        # Theta
        self.model_opts['learnTheta'][m][0] = 0.
        self.model_opts["initSW"]["Theta"][m][:,0] = 1.
        self.model_opts["priorTheta"]['a'][m][0] = s.nan
        self.model_opts["priorTheta"]['b'][m][0] = s.nan
        self.model_opts["initTheta"]["a"][m][0] = s.nan
        self.model_opts["initTheta"]["b"][m][0] = s.nan
        self.model_opts["initTheta"]["E"][m][:,0] = 1.
