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