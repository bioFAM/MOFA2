class MuZ_Node(UnivariateGaussian_Unobserved_Variational_Node):
    """ """

    def __init__(self, pmean, pvar, qmean, qvar, clusters, n_Z, cluster_dic=None, qE=None, qE2=None):
        # For now clusters have to be integers from 0 to n_clusters
        # compute dim from numbers of clusters (n_clusters * Z)
        self.clusters = clusters
        self.N = len(self.clusters)
        self.n_clusters = len(np.unique(clusters))
        dim = (self.n_clusters, n_Z)
        self.factors_axis = 1
        super(Cluster_Node_Gaussian, self).__init__(dim=dim, pmean=pmean, pvar=pvar, qmean=qmean, qvar=qvar, qE=qE,
                                                    qE2=qE2)

    def getExpectations(self):
        # reshape the values to N_samples * N_factors and return
        QExp = self.Q.getExpectations()
        expanded_expectation = QExp['E'][self.clusters, :]
        expanded_E2 = QExp['E2'][self.clusters, :]
        # do we need to expand the variance as well -> not used I think
        return {'E': expanded_expectation, 'E2': expanded_E2}

    def updateParameters(self):
        Ppar = self.P.getParameters()
        Z = self.markov_blanket['Z'].Q.getExpectation()

        if "AlphaZ" in self.markov_blanket:
            Alpha = self.markov_blanket[
                'AlphaZ'].getExpectation().copy()  # Notice that this Alpha is the ARD prior on Z, not on W.
        elif "SigmaZ" in self.markov_blanket:
            Sigma = self.markov_blanket['SigmaZ'].getExpectation().copy()
        else:
            Sigma = self.markov_blanket['Z'].P.getParameters()["cov"]

        Qmean, Qvar = self.Q.getParameters()['mean'], self.Q.getParameters()['var']

        # update of the variance

        for c in range(self.n_clusters):
            mask = (self.clusters == c)
            if "AlphaZ" in self.markov_blanket:
                tmp = (Alpha[mask, :]).sum(axis=0)
            else:
                tmp = np.matrix.trace(Sigma[:, mask, mask])
            Qvar[c, :] = tmp
        Qvar += 1. / Ppar['var']
        Qvar = 1. / Qvar

        # update of the mean

        if "AlphaZ" in self.markov_blanket:
            tmp = Z * Alpha
        else:
            # TODO : check if we should ask the mask l. 462
            tmp = np.zeros(self.dim)
            for k in range(self.dim[1]):
                tmp[:, k] = np.dot(Sigma[k, :, :] - np.diag(Sigma[k, :, :]), Z[:, k])

        for c in range(self.n_clusters):
            mask = (self.clusters == c)
            tmp = (tmp[mask, :]).sum(axis=0)
            Qmean[c, :] = tmp
        Qmean = Qmean + Ppar['mean'] / Ppar['var']
        Qmean *= Qvar

        self.Q.setParameters(mean=Qmean, var=Qvar)

    def calculateELBO(self):
        PParam = self.P.getParameters()
        PVar, Pmean = PParam['var'], PParam['mean']

        QExp = self.Q.getExpectations()
        QE2, QE = QExp['E2'], QExp['E']

        Qvar = self.Q.getParameters()['var']

        # Cluster terms corresponding to covariates should not intervene
        # filtering the covariates out
        latent_variables = self.markov_blanket['Z'].getLvIndex()
        PVar, Pmean = PVar[:, latent_variables], Pmean[:, latent_variables]
        QE2, QE = QE2[:, latent_variables], QE[:, latent_variables]
        Qvar = Qvar[:, latent_variables]

        # minus cross entropy
        tmp = -(0.5 * s.log(PVar)).sum()
        tmp2 = - ((0.5 / PVar) * (QE2 - 2. * QE * Pmean + Pmean ** 2.)).sum()

        # entropy of Q
        tmp3 = 0.5 * (s.log(Qvar)).sum()
        tmp3 += 0.5 * self.dim[0] * len(latent_variables)

        return tmp + tmp2 + tmp3



# Z node with multivariate prior
# class Z_Node(UnivariateGaussian_Unobserved_Variational_Node_with_MultivariateGaussian_Prior):
#     def __init__(self, dim, pmean, pcov, qmean, qvar, qE=None, qE2=None, idx_covariates=None, precompute_pcovinv=True):
#         super().__init__(dim=dim, pmean=pmean, pcov=pcov, qmean=qmean, qvar=qvar, axis_cov=0, qE=qE, qE2=qE2)
#
#         self.precompute_pcovinv = precompute_pcovinv
#
#         self.N = self.dim[0]
#         self.K = self.dim[1]
#
#         # Define indices for covariates
#         if idx_covariates is not None:
#             self.covariates[idx_covariates] = True
#
#     def precompute(self):
#         # Precompute terms to speed up computation
#         self.covariates = np.zeros(self.dim[1], dtype=bool)
#         self.factors_axis = 1
#
#         if self.precompute_pcovinv:
#             p_cov = self.P.params["cov"]
#
#             self.p_cov_inv = []
#             self.p_cov_inv_diag = []
#
#             for k in range(self.K):
#                 if p_cov[k].__class__.__name__ == 'dia_matrix':
#                     diag_inv = 1 / p_cov[k].data
#                     diag_inv = diag_inv.flatten()
#                     inv = np.diag(diag_inv)
#                 elif p_cov[k].__class__.__name__ == 'ndarray':
#                     diagonal = np.diagonal(p_cov[k])
#                     if np.all(np.diag(diagonal) == p_cov[k]):
#                         diag_inv = 1. / diagonal
#                         inv = np.diag(diag_inv)
#                     else:
#                         inv = np.linalg.inv(p_cov[k])
#                         diag_inv = np.diagonal(inv)
#                 else:
#                     #TODO : deal with sparse non diagonal input matrices as pcov
#                     print("Not implemented yet")
#                     exit()
#
#                 self.p_cov_inv.append(inv)
#                 self.p_cov_inv_diag.append(diag_inv)
#
#         else:
#             self.p_cov_inv = None
#             self.p_cov_inv_diag = None
#
#     def getLvIndex(self):
#         # Method to return the index of the latent variables (without covariates)
#         latent_variables = np.array(range(self.dim[1]))
#         if any(self.covariates):
#             # latent_variables = np.delete(latent_variables, latent_variables[self.covariates])
#             latent_variables = latent_variables[~self.covariates]
#         return latent_variables
#
#     def removeFactors(self, idx, axis=1):
#         super(Z_Node, self).removeFactors(idx, axis)
#
#         if self.p_cov_inv is not None:
#             for i in idx :
#                 del self.p_cov_inv[i]
#                 del self.p_cov_inv_diag[i]
#
#     def updateParameters(self):
#         # Collect expectations from the markov blanket
#         Y = self.markov_blanket["Y"].getExpectation()
#         SWtmp = self.markov_blanket["W"].getExpectations()
#         tau = self.markov_blanket["Tau"].getExpectation()
#         latent_variables = self.getLvIndex()  # excluding covariates from the list of latent variables
#         mask = [ma.getmask(Y[m]) for m in range(len(Y))]
#
#         # Collect parameters from the prior or expectations from the markov blanket
#         if "MuZ" in self.markov_blanket:
#             Mu = self.markov_blanket['MuZ'].getExpectation()
#         else:
#             Mu = self.P.getParameters()["mean"]
#
#         if "AlphaZ" in self.markov_blanket:
#             Alpha = self.markov_blanket['AlphaZ'].getExpectation(expand=True)
#
#         else:
#             if "SigmaZ" in self.markov_blanket:
#                 Sigma = self.markov_blanket['SigmaZ'].getExpectations()
#                 p_cov_inv = Sigma['inv']
#                 p_cov_inv_diag = Sigma['inv_diag']
#             else:
#                 p_cov_inv = self.p_cov_inv
#                 p_cov_inv_diag = self.p_cov_inv_diag
#
#         # Check dimensionality of Tau and expand if necessary (for Jaakola's bound only)
#         for m in range(len(Y)):
#             # Mask tau
#             # tau[m] = ma.masked_where(ma.getmask(Y[m]), tau[m]) # important to keep this out of the loop to mask non-gaussian tau
#             tau[m][mask[m]] = 0.
#             # Mask Y
#             Y[m] = Y[m].data
#             Y[m][mask[m]] = 0.
#
#         # Collect parameters from the P and Q distributions of this node
#         Q = self.Q.getParameters()
#         Qmean, Qvar = Q['mean'], Q['var']
#
#         M = len(Y)
#         for k in latent_variables:
#             foo = s.zeros((self.N,))
#             bar = s.zeros((self.N,))
#             for m in range(M):
#                 foo += np.dot(tau[m], SWtmp[m]["E2"][:, k])
#
#                 bar_tmp1 = SWtmp[m]["E"][:,k]
#
#                 # NOTE slow bit but hard to optimise
#                 # bar_tmp2 = - fast_dot(Qmean[:, s.arange(self.dim[1]) != k], SWtmp[m]["E"][:, s.arange(self.dim[1]) != k].T)
#                 bar_tmp2 = - s.dot(Qmean[:, s.arange(self.dim[1]) != k], SWtmp[m]["E"][:, s.arange(self.dim[1]) != k].T)
#                 bar_tmp2 += Y[m]
#                 bar_tmp2 *= tau[m]
#                 ##############################
#
#                 bar += np.dot(bar_tmp2, bar_tmp1)
#
#             if "AlphaZ" in self.markov_blanket:
#                 Qvar[:, k] = 1. / (Alpha[:, k] + foo)
#                 Qmean[:, k] = Qvar[:, k] * (bar + Alpha[:, k] * Mu[:, k])
#
#             else:
#                 Qvar[:, k] = 1. / (foo + p_cov_inv_diag[k])
#
#                 if self.P.params["cov"][k].__class__.__name__ == 'dia_matrix':
#                     Qmean[:, k] = Qvar[:, k] * bar
#                 else:
#                     tmp = p_cov_inv[k] - p_cov_inv_diag[k] * s.eye(self.N)
#                     for n in range(self.N):
#                         Qmean[n, k] = Qvar[n, k] * (bar[n] + np.dot(tmp[n, :], Mu[:, k] - Qmean[:, k]))
#
#         # Save updated parameters of the Q distribution
#         self.Q.setParameters(mean=Qmean, var=Qvar)
#
#     # TODO, problem here is that we need to make sure k is in the latent variables first
#     def calculateELBO_k(self, k):
#         '''Compute the ELBO for factor k in absence of Alpha node in the markov blanket of Z'''
#         Qpar, Qexp = self.Q.getParameters(), self.Q.getExpectations()
#         Qmean, Qvar = Qpar['mean'], Qpar['var']
#         QE, QE2 = Qexp['E'], Qexp['E2']
#
#         if "SigmaZ" in self.markov_blanket:
#             Sigma = self.markov_blanket['SigmaZ'].getExpectations()
#             p_cov = Sigma['cov']
#             p_cov_inv = Sigma['inv']
#             p_cov_inv_diag = Sigma['inv_diag']
#         else:
#             p_cov = self.P.params['cov']
#             p_cov_inv = self.p_cov_inv
#             p_cov_inv_diag = self.p_cov_inv_diag
#
#         # compute cross entropy term
#         tmp1 = 0
#         if p_cov[k].__class__.__name__ == 'ndarray':
#             mat_tmp = p_cov_inv[k] - p_cov_inv_diag[k] * s.eye(self.N)
#             tmp1 += QE[:, k].transpose().dot(mat_tmp).dot(QE[:, k])
#         tmp1 += p_cov_inv_diag[k].dot(QE2[:, k])
#         tmp1 = -.5 * tmp1
#         # tmp1 = 0.5*QE2 - PE*QE + 0.5*PE2
#         # tmp1 = -(tmp1 * Alpha['E']).sum()
#
#         # compute term from the precision factor in front of the Gaussian
#         tmp2 = 0  # constant here
#         # if self.n_iter> 4:
#         if p_cov[k].__class__.__name__ == 'dia_matrix':
#             tmp2 += np.sum(np.log(p_cov[k].data.flatten()))
#         elif p_cov[k].__class__.__name__ == 'ndarray':
#             tmp2 += np.linalg.slogdet(p_cov[k])[1]
#         else:
#             print("Not implemented yet")
#             exit()
#         tmp2 *= (-.5)
#         # tmp2 = 0.5*Alpha["lnE"].sum()
#
#         lb_p = tmp1 + tmp2
#         # lb_q = -(s.log(Qvar).sum() + self.N*self.dim[1])/2. # I THINK THIS IS WRONG BECAUSE SELF.DIM[1] ICNLUDES COVARIATES
#         lb_q = -.5 * s.log(Qvar[:, k]).sum()
#
#         return lb_p - lb_q
#
#     def calculateELBO(self):
#         if not ("AlphaZ" in self.markov_blanket):
#             latent_variables = self.getLvIndex()
#             elbo = 0
#             for k in latent_variables:
#                 elbo += self.calculateELBO_k(k)
#
#             elbo += .5 * self.N * len(latent_variables)
#
#             return elbo
#
#         else:
#             # Collect parameters and expectations of current node
#             Qpar, Qexp = self.Q.getParameters(), self.Q.getExpectations()
#             Qmean, Qvar = Qpar['mean'], Qpar['var']
#             QE, QE2 = Qexp['E'], Qexp['E2']
#
#             if "MuZ" in self.markov_blanket:
#                 PE, PE2 = self.markov_blanket['MuZ'].getExpectations()['E'], \
#                           self.markov_blanket['MuZ'].getExpectations()['E2']
#             else:
#                 PE, PE2 = self.P.getParameters()["mean"], s.zeros((self.N, self.dim[1]))
#
#             Alpha = self.markov_blanket[
#                 'AlphaZ'].getExpectations(expand=True).copy()  # Notice that this Alpha is the ARD prior on Z, not on W.
#
#             # This ELBO term contains only cross entropy between Q and P,and entropy of Q. So the covariates should not intervene at all
#             latent_variables = self.getLvIndex()
#             Alpha["E"], Alpha["lnE"] = Alpha["E"][:, latent_variables], Alpha["lnE"][:, latent_variables]
#             Qmean, Qvar = Qmean[:, latent_variables], Qvar[:, latent_variables]
#             PE, PE2 = PE[:, latent_variables], PE2[:, latent_variables]
#             QE, QE2 = QE[:, latent_variables], QE2[:, latent_variables]
#
#             # compute term from the exponential in the Gaussian
#             tmp1 = 0.5 * QE2 - PE * QE + 0.5 * PE2
#             tmp1 = -(tmp1 * Alpha['E']).sum()
#
#             # compute term from the precision factor in front of the Gaussian
#             tmp2 = 0.5 * Alpha["lnE"].sum()
#
#             lb_p = tmp1 + tmp2
#             # lb_q = -(s.log(Qvar).sum() + self.N*self.dim[1])/2. # I THINK THIS IS WRONG BECAUSE SELF.DIM[1] ICNLUDES COVARIATES
#             lb_q = -(s.log(Qvar).sum() + self.N * len(latent_variables)) / 2.
#
#             return lb_p - lb_q
#
#     def sample(self, dist='P'):
#         if "MuZ" in self.markov_blanket:
#             p_mean = self.markov_blanket['MuZ'].sample()
#         else:
#             p_mean = self.P.params['mean']
#
#         if "AlphaZ" in self.markov_blanket:
#             alpha = self.markov_blanket['AlphaZ'].sample()
#             p_var = s.square(1. / alpha)
#             #p_cov = s.diag(p_var)
#             p_cov = [p_var[k] * np.eye(self.N) for k in range(self.K)]
#         else:
#             if "SigmaZ" in self.markov_blanket:
#                 p_cov = self.markov_blanket['SigmaZ'].sample()
#             else:
#                 p_cov = self.P.params['cov']
#
#         # simulating
#
#         samp_tmp = []
#         for i in range(self.dim[1]):
#             # samp_tmp.append(s.random.multivariate_normal(p_mean[:, i], p_cov[i])) #does it yield the correct result for sparse input matrices ?
#             if p_cov[i].__class__.__name__ == 'dia_matrix':
#                 samp_tmp.append(s.random.multivariate_normal(p_mean[:, i], p_cov[i].toarray()))  # inefficient
#             elif p_cov[i].__class__.__name__ == 'ndarray':
#                 samp_tmp.append(s.random.multivariate_normal(p_mean[:, i], p_cov[i]))
#             else:
#                 print("Not implemented yet")
#                 exit()
#
#         # self.samp = s.array([tmp/tmp.std() for tmp in samp_tmp]).transpose()
#
#         self.samp = s.array([tmp - tmp.mean() for tmp in samp_tmp]).transpose()
#         #self.samp = np.array(samp_tmp).T
#
#         return self.samp
