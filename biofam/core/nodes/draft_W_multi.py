#Below, changes for W multivariate gaussian in W_nodes.py, init_model.py, build_model.py
#We would have also to change updates for Y and tau nodes (since E(W1.W2)!=E(W1)E(W2) when breaking the squared sums... sum out of the
#other factors), and maybe also AlphaW and ThetaW...

#1. in W nodes.py

class W_Node(MultivariateGaussian_Unobserved_Variational_Node):
   def __init__(self, dim, pmean, pcov, qmean, qcov, qE=None, qE2=None):
       MultivariateGaussian_Unobserved_Variational_Node.__init__(self, dim=dim,  pmean=pmean, pcov=pcov, qmean=qmean, qcov=qcov, qE=qE)
       self.precompute()

   def precompute(self):
       self.D = self.dim[0]
       self.K = self.dim[1]
       self.factors_axis = 1

   def updateParameters(self):
       Z = self.markov_blanket["TZ"].getExpectation()
       ZZ = self.markov_blanket["TZ"].getExpectations()["EXXT"]
       alpha = self.markov_blanket["AlphaW"].getExpectation()

       #tau = (self.markov_blanket["Tau"].getExpectation())[:, None, None]
       tau = (self.markov_blanket["Tau"].getExpectation())
       Y = self.markov_blanket["Y"].getExpectation()

       #Qcov = s.linalg.inv(tau * s.repeat(ZZ[None, :, :], self.D, 0) + s.diag(alpha))
       Qcov = np.zeros((self.D,self.K,self.K))

       for d in range(self.D):
           Qcov[d,:,:] = s.linalg.inv(tau[d] * np.sum(ZZ,axis=0) + s.diag(alpha))
       tmp1 = s.repeat(s.repeat(tau[:,None,None],self.K,axis=1),self.K,axis=2) * Qcov
       tmp2 = np.dot(Y.T, Z) #tmp2 = ma.dot(Y.T, Z).data
       Qmean = (tmp1[:, :, :] * tmp2[:, None, :]).sum(axis=2)

       # Save updated parameters of the Q distribution
       self.Q.setParameters(mean=Qmean, cov=Qcov)

   def calculateELBO(self):

       # Collect parameters and expectations
       alpha = self.markov_blanket["AlphaW"].getExpectations()["E"]
       logalpha = self.markov_blanket["AlphaW"].getExpectations()["lnE"]
       Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()

       lb_p = self.D*s.sum(logalpha) - s.sum(Qexp['EXXT'] * s.diag(alpha)[None,:,:])
       lb_q = -self.D*self.K - logdet(Qpar['cov']).sum()

       return (lb_p - lb_q)/2

#2. in init.py

    def initW(self, pmean, pcov, qmean, qcov, qE=None, qE2=None):
        """Method to initialise the weights
        PARAMETERS
        ----------
        pmean: initial value of the mean of the prior distribution
        pcov: initial value of the covariance of the prior distribution
        qmean: initial value of the mean of the variational distribution
        qcov: initial value of the covariance of the variational distribution
        qE: initial value of the expectation of the variational distribution
        qE2: initial value of the second moment of the variational distribution
        """
        W_list = [None] * self.M
        for m in range(self.M):

            # mean
            pmean = s.ones((self.D[m], self.K)) * pmean

            # covariance
            if isinstance(pcov, (int, float)):
                tmp = s.zeros((self.D[m],self.K,self.K))
                for d in range(self.D[m]):
                    tmp[d,:,:]=np.eye(self.K)*pcov
                pcov=tmp

            ## Initialise variational distribution (Q) ##

            # mean
            if qmean=="random":
                # Random initialisation
                qmean = stats.norm.rvs(loc=0, scale=1, size=(self.D[m], self.K))

            elif isinstance(qmean, s.ndarray):
                assert qmean.shape == (self.D[m], self.K), "Wrong shape for the expectation of the Q distribution of Z"

            elif isinstance(qmean, (int, float)):
                qmean = s.ones((self.D[m],self.K)) * qmean

            else:
                print("Wrong initialisation for W")
                exit()

            # covariance

            if isinstance(qcov, (int, float)):
                tmp = s.zeros((self.D[m],self.K,self.K))
                for d in range(self.D[m]):
                    tmp[d,:,:]=np.eye(self.K)*qcov
                qcov=tmp

            # variance
            # qvar = s.ones((self.N, self.D[m])) * qvar

            W_list[m] = W_Node(dim=(self.D[m], self.K),pmean=pmean,pcov=pcov,qmean=qmean, qcov=qcov)

        self.nodes["W"] = Multiview_Variational_Node(self.M, *W_list)

#3. in build_model.py

#Initialise weights
    if model_opts['transpose']:
        pmean = 0.
        pcov = 1. #to change
        qmean = "random"
        qcov = 1. #to change
        init.initW(pmean=pmean,pcov=pcov,qmean=qmean, qcov=qcov)


#below, not related directly
'''
in W nodes
def getExpectation(self, dist="Q"):
    if dist == "Q":
        expectations = self.Q.getExpectations()
    elif dist == "P":
        expectations = self.P.getExpectations()
    return np.transpose(expectations["E"]) #since we used the transpose version of multivariate_gaussian

def getExpectations(self, dist="Q"):
    if dist == "Q":
        tmp = self.Q.getExpectations()
    elif dist == "P":
        tmp = self.P.getExpectations()
    expectations={}
    for (k,v) in tmp.items():
        print(k,v.shape)
        expectations[k]=np.transpose(v) #since we used the transpose version of multivariate_gaussian
    return expectations
'''
'''
#in univariate_gaussian
#TO : remove this loop
EXXT = np.zeros((self.dim[0], self.dim[1], self.dim[1]))
for i in range(self.dim[0]):
    EXXT[i, :, :] = np.diag(self.params['var'][i, :]) + np.dot(self.expectations['E'][i, :].T, self.expectations['E'][i, :])
self.expectations['EXXT'] = EXXT
'''
'''
class W_Node(UnivariateGaussian_Unobserved_Variational_Node):
    def __init__(self, dim, pmean, pvar, qmean, qvar, qE=None, qE2=None):
        super(W_Node,self).__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE, qE2=qE2)
        self.precompute()

    def precompute(self):
        # Precompute terms to speed up computation
        self.D = self.dim[0]
        self.factors_axis = 1

    def updateParameters(self):

        # Collect expectations from the markov blanket
        Y = deepcopy(self.markov_blanket["Y"].getExpectation())
        TZtmp = self.markov_blanket["TZ"].getExpectations()
        tau = deepcopy(self.markov_blanket["Tau"].getExpectation())

        if "AlphaW" in self.markov_blanket:
            Alpha = self.markov_blanket['AlphaW'].getExpectation()
            Alpha = s.repeat(Alpha[None,:], self.D, axis=0)
        else:
            Alpha = 1./self.P.getParameters()["var"]

        # Check dimensionality of Tau and expand if necessary (for Jaakola's bound only)
        N = np.shape(Y[0])[0]
        for m in range(len(Y)):
            if tau[m].shape != Y[m].shape:
                tau[m] = s.repeat(tau[m].copy()[None,:], N, axis=0)

        # Collect parameters from the P and Q distributions of this node
        Q = self.Q.getParameters().copy()
        Qmean, Qvar = Q['mean'], Q['var']

        M = len(Y)
        K = np.shape(TZtmp["E"])[1]
        for k in range(K):
            foo = s.zeros((N,))
            bar = s.zeros((N,))
            for m in range(M):
                foo += np.dot(tau[m],TZtmp["EBNN"][:,k])
                bar += np.dot(tau[m]*(Y[m] - s.dot( Qmean[:,s.arange(self.dim[1])!=k] , TZtmp["E"][:,s.arange(self.dim[1])!=k].T )), TZtmp["E"][:,k])
            Qvar[:,k] = 1./(Alpha[:,k]+foo)
            Qmean[:,k] = Qvar[:,k] * (  Alpha[:,k]*Mu[:,k] + bar )

        # Save updated parameters of the Q distribution
        self.Q.setParameters(mean=Qmean, var=Qvar)

    def calculateELBO(self):
        # Collect parameters and expectations of current node
        Qpar,Qexp = self.Q.getParameters(), self.Q.getExpectations()
        Qmean, Qvar = Qpar['mean'], Qpar['var']
        QE, QE2 = Qexp['E'],Qexp['E2']
        PE, PE2 = self.P.getParameters()["mean"], s.zeros((self.N,self.dim[1]))

        if "AlphaZ" in self.markov_blanket:
            Alpha = self.markov_blanket['AlphaZ'].getExpectations().copy() # Notice that this Alpha is the ARD prior on Z, not on W.
            Alpha["E"] = s.repeat(Alpha["E"][None,:], self.N, axis=0)
            Alpha["lnE"] = s.repeat(Alpha["lnE"][None,:], self.N, axis=0)
        else:
            Alpha = { 'E':1./self.P.getParameters()["var"], 'lnE':s.log(1./self.P.getParameters()["var"]) }

        # compute term from the exponential in the Gaussian
        tmp1 = 0.5*QE2 - PE*QE + 0.5*PE2
        tmp1 = -(tmp1 * Alpha['E']).sum()

        # compute term from the precision factor in front of the Gaussian
        tmp2 = 0.5*Alpha["lnE"].sum()

        lb_p = tmp1 + tmp2
        lb_q = -(s.log(Qvar).sum() + self.N)/2.

        return lb_p-lb_q
'''