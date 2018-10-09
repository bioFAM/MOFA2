def initAlphaZ_noGroups(self, pa=1e-14, pb=1e-14, qa=1., qb=1., qE=None, qlnE=None):
    """Method to initialise the ARD prior on Z per factor

    PARAMETERS
    ----------
    pa: float
        'a' parameter of the prior distribution
    pb :float
        'b' parameter of the prior distribution
    qb: float
        initialisation of the 'b' parameter of the variational distribution
    qE: float
        initial expectation of the variational distribution
    qlnE: float
        initial log expectation of the variational distribution
    """
    self.nodes["AlphaZ"] = AlphaZ_Node_k(dim=(self.K,), pa=pa, pb=pb, qa=qa, qb=qb, qE=qE, qlnE=qlnE)


def initSigmaZ(self, X, n_diag=0):
    """Method to initialise the covariance prior structure on Z

    (TO-DO)
    PARAMETERS
    ----------
    X:
    n_diag:
    """
    dim = (self.K,)
    self.Sigma = SigmaGrid_Node(dim, X, n_diag=n_diag)
    self.nodes["SigmaZ"] = self.Sigma

def initSigmaZ_Block(self, X, clust, n_diag=0):
    """Method to initialise the covariance prior structure on Z, for clusters assigned to samples

    (TO-DO)
    PARAMETERS
    ----------
    X:
    n_diag:
    """
    dim = (self.K,)
    self.Sigma = BlockSigmaGrid_Node(dim, X, clust, n_diag=n_diag)
    self.nodes["SigmaZ"] = self.Sigma


def initSigmaW_Mixed(self,view_has_covariance_prior,params_views):
    '''Method to initialise the covariance prior structure or the Gamma variance on W for each view
    For each view, we choose between a covariance prior structure or a Gamma variance.

    PARAMETERS :
    -----------
    view_has_covariance_prior : list with for each view a boolean which takes value True if a covariance prior
    structure is choosen (or False if Gamma variance choosen)

    params_views : list with for each view the dict {pa ; pb ; qa ; qb; qE} of the Gamma variance
    or the dict {X, clust, n_diag} of the covariance prior structure
    '''

    AlphaSigmaNodes = [None] * self.M

    for m in range(self.M):

        params = params_views[m]
        dim = (self.K,)

        if view_has_covariance_prior[m]:
            # TODO add a if statement to check if there is a sigma_clust argument to see if blockSigma is needed
            if params['sigma_clust'] is None:
                AlphaSigmaNodes[m]=SigmaGrid_Node(dim,params['X'], n_diag = params['n_diag'])
            else:
                AlphaSigmaNodes[m]=SigmaBlockW_k(dim,params['X'], clust = params['sigma_clust'], n_diag = params['n_diag'])

        else:
            AlphaSigmaNodes[m] = AlphaW_Node_mk(dim,pa = params['pa'], pb = params['pb'], qa = params['qa'], qb = params['qb'])

    self.nodes["SigmaAlphaW"] = Multiview_Mixed_Node(self.M, *AlphaSigmaNodes)

def initThetaZ_Learn(self, pa=1., pb=1., qa=1., qb=1., qE=None):
    """Method to initialise the sparsity parameter of the spike and slab factors

    PARAMETERS
    ----------
     pa: float
        'a' parameter of the prior distribution
     pb :float
        'b' parameter of the prior distribution
     qb: float
        initialisation of the 'b' parameter of the variational distribution
     qE: float
        initial expectation of the variational distribution
    """
    self.nodes["ThetaZ"] = ThetaZ_Node_k(dim=(self.K,), pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)

def initThetaZ_Mixed(self, idx, groups = None, pa=1., pb=1., qa=1., qb=1., qE=1.):
    """Method to initialise the sparsity parameter of the spike and slab factors
    In contrast with initThetaLearn, where a sparsity parameter is learnt by each feature and factor, and initThetaConst, where the sparsity is not learnt,
    in initThetaMixed the sparsity parameter is learnt by a subset of factors and features

    PARAMETERS
    ----------
     pa: float
        'a' parameter of the prior distribution
     pb :float
        'b' parameter of the prior distribution
     qb: float
        initialisation of the 'b' parameter of the variational distribution
     qE: (...)
    idx:list with binary matrices with dim (N,K)

    """

    # Do some sanity checks on the arguments
    if isinstance(qE, (int, float)):
        qE = s.ones((self.N, self.K)) * qE
    elif isinstance(qE, s.ndarray):
        assert qE.shape == (self.N, self.K), "Wrong dimensionality of Theta"

    # Initialise constant node
    Kconst = idx == 0
    if Kconst.sum() == 0:
        ConstThetaNode = None
    else:
        if qE is None:
            print("Wrong initialisation for Theta");
            exit(1)
        else:
            ConstThetaNode = ThetaZ_Constant_Node_k(dim=(self.N, s.sum(Kconst),), value=qE[:, Kconst], N_cells=1)
            self.nodes["ThetaZ"] = ConstThetaNode

    # Initialise non-constant node
    Klearn = idx == 1
    if Klearn.sum() == 0:
        LearnThetaNode = None
    else:
        # FOR NOW WE JUST TAKE THE FIRST ROW BECAUSE IT IS EXPANDED, THIS IS UGLY
        if qE is None:
            LearnThetaNode = ThetaZ_Node_k(dim=(s.sum(Klearn),), pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        else:
            LearnThetaNode = ThetaZ_Node_k(dim=(s.sum(Klearn),), pa=pa, pb=pb, qa=qa, qb=qb, qE=qE[0, Klearn])
        self.nodes["ThetaZ"]=LearnThetaNode

    # Initialise mixed node
    if (ConstThetaNode is not None) and (LearnThetaNode is not None):
        self.nodes["ThetaZ"] = Mixed_ThetaZ_Nodes_k(LearnTheta=LearnThetaNode, ConstTheta=ConstThetaNode, idx=idx)

# TODO fix that to account for groups when learning the theta
def initThetaZ_Mixed_groups(self, idx, groups, pa=1., pb=1., qa=1., qb=1., qE=1.):
    """Method to initialise the sparsity parameter of the spike and slab factors
    In contrast with initThetaLearn, where a sparsity parameter is learnt by each feature and factor, and initThetaConst, where the sparsity is not learnt,
    in initThetaMixed the sparsity parameter is learnt by a subset of factors and features

    PARAMETERS
    ----------
     pa: float
        'a' parameter of the prior distribution
     pb :float
        'b' parameter of the prior distribution
     qb: float
        initialisation of the 'b' parameter of the variational distribution
     qE: (...)
    idx:list with binary matrices with dim (N,K)

    """

    # Do some sanity checks on the arguments
    if isinstance(qE, (int, float)):
        qE = s.ones((self.N, self.K)) * qE
    elif isinstance(qE, s.ndarray):
        assert qE.shape == (self.N, self.K), "Wrong dimensionality of Theta"

    # Initialise constant node
    Kconst = idx == 0
    if Kconst.sum() == 0:
        ConstThetaNode = None
    else:
        if qE is None:
            print("Wrong initialisation for Theta");
            exit(1)
        else:
            ConstThetaNode = ThetaZ_Constant_Node_k(dim=(self.N, s.sum(Kconst),), value=qE[:, Kconst], N_cells=1)
            self.nodes["ThetaZ"] = ConstThetaNode

    # Initialise non-constant node
    Klearn = idx == 1
    if Klearn.sum() == 0:
        LearnThetaNode = None
    else:
        # FOR NOW WE JUST TAKE THE FIRST ROW BECAUSE IT IS EXPANDED, THIS IS UGLY
        if qE is None:
            LearnThetaNode = ThetaZ_Node_groups(groups, K=Klearn, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        else:
            # TODO not sure that works here
            LearnThetaNode = ThetaZ_Node_groups(groups, K=Klearn, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE[0, Klearn])
        self.nodes["ThetaZ"]=LearnThetaNode

    # Initialise mixed node
    if (ConstThetaNode is not None) and (LearnThetaNode is not None):
        self.nodes["ThetaZ"] = Mixed_ThetaZ_Nodes_k(LearnTheta=LearnThetaNode, ConstTheta=ConstThetaNode, idx=idx)

def initThetaZ_Const(self, value=1.):
    """Method to initialise a constant sparsity parameter of the spike and slab factors

    PARAMETERS
    ----------
     value: ndarray
        constant value from 0 to 1 to initialise the node, 0 corresponds to complete sparsity (all weights are zero) and 1 corresponds to no sparsity
    """

    self.nodes["ThetaZ"] = ThetaZ_Constant_Node_k(dim=(self.N, self.K,), value=s.ones((self.N, self.K)) * value,
                                                  N_cells=1.)

def initThetaW_Mixed(self, idx, pa=1., pb=1., qa=1., qb=1., qE=1.):
    """Method to initialise the sparsity parameter of the spike and slab weights
    In contrast with initThetaLearn, where a sparsity parameter is learnt by each feature and factor, and initThetaConst, where the sparsity is not learnt,
    in initThetaMixed the sparsity parameter is learnt by a subset of factors and features

    PARAMETERS
    ----------
     pa: float
        'a' parameter of the prior distribution
     pb :float
        'b' parameter of the prior distribution
     qb: float
        initialisation of the 'b' parameter of the variational distribution
     qE: (...)
    idx:list with binary matrices with dim (D[m],K)

    """

    # Do some sanity checks on the arguments
    if isinstance(qE,list):
        assert len(qE) == self.M, "Wrong dimensionality"
        for m in range(self.M):
            if isinstance(qE[m],(int,float)):
                qE[m] = s.ones((self.D[m],self.K)) * qE[m]
            elif isinstance(qE[m],s.ndarray):
                assert qE[m].shape == (self.D[m],self.K), "Wrong dimensionality of Theta"
            else:
                print("Wrong initialisation for Theta"); exit(1)

    elif isinstance(qE,s.ndarray):
        assert qE.shape == (self.D[m],self.K), "Wrong dimensionality of Theta"
        tmp = [ qE for m in range(self.M)]
        qE = tmp # IS THIS REQUIRED????


    elif isinstance(qE,(int,float)):
        tmp = [ s.ones((self.D[m],self.K)) * qE for m in range(self.M)]
        qE = tmp # IS THIS REQUIRED????

    Theta_list = [None] * self.M
    for m in range(self.M):

        # Initialise constant node
        Kconst = idx[m]==0
        if Kconst.sum() == 0:
            ConstThetaNode = None
        else:
            if qE is None:
                print("Wrong initialisation for Theta");
                exit(1)
            else:
                ConstThetaNode = ThetaW_Constant_Node_mk(dim=(self.D[m],s.sum(Kconst),), value=qE[m][:,Kconst], N_cells=1)
                Theta_list[m] = ConstThetaNode

        # Initialise non-constant node
        Klearn = idx[m]==1
        if Klearn.sum() == 0:
            LearnThetaNode = None
        else:
            if qE is None:
                LearnThetaNode = ThetaW_Node_mk(dim=(s.sum(Klearn),), pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
            else:
                # FOR NOW WE JUST TAKE THE FIRST ROW BECAUSE IT IS EXPANDED, THIS IS UGLY
                LearnThetaNode = ThetaW_Node_mk(dim=(s.sum(Klearn),), pa=pa, pb=pb, qa=qa, qb=qb, qE=qE[m][0,Klearn])
            Theta_list[m] = LearnThetaNode

        # Initialise mixed node
        if (ConstThetaNode is not None) and (LearnThetaNode is not None):
            Theta_list[m] = Mixed_ThetaW_Nodes_mk(LearnTheta=LearnThetaNode, ConstTheta=ConstThetaNode, idx=idx[m])

    self.Theta = Multiview_Mixed_Node(self.M, *Theta_list)
    self.nodes["ThetaW"] = self.Theta

def initThetaW_Const(self, value=1.):
    """Method to initialise a constant Theta of the spike and slab on W

    PARAMETERS
    ----------
    value: float ranging from 0 to 1, where:
        0 corresponds to complete sparsity (all weights are zero)
        1 corresponds to no sparsity (all weight allowed to be non-zero)
    """
    Theta_list = [None] * self.M
    for m in range(self.M):
        Theta_list[m] = ThetaW_Constant_Node_mk(dim=(self.D[m],self.K,), value=s.ones((self.D[m],self.K))*value, N_cells=1.)
    self.Theta = Multiview_Constant_Node(self.M, *Theta_list)
    self.nodes["ThetaW"] = self.Theta
