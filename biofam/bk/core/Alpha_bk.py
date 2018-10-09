class AlphaZ_Node_k(Gamma_Unobserved_Variational_Node):
    def __init__(self, dim, pa, pb, qa, qb, qE=None, qlnE=None):
        # Gamma_Unobserved_Variational_Node.__init__(self, dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE)
        super().__init__(dim=dim, pa=pa, pb=pb, qa=qa, qb=qb, qE=qE, qlnE=qlnE)

    def precompute(self, options=None):
        # self.lbconst = self.K * ( self.P.a*s.log(self.P.b) - special.gammaln(self.P.a) )
        # self.lbconst = s.sum( self.P.params['a']*s.log(self.P.params['b']) - special.gammaln(self.P.params['a']) )
        self.factors_axis = 0   # TODO check that !!!!!

    def getExpectations(self, expand=False):
        QExp = self.Q.getExpectations()
        if expand:
            N = self.markov_blanket['Z'].N
            expanded_E = s.repeat(QExp['E'][None, :], N, axis=0)
            expanded_lnE = s.repeat(QExp['lnE'][None, :], N, axis=0)
            return {'E': expanded_E, 'lnE': expanded_lnE}
        else:
            return QExp

    def getExpectation(self, expand=False):
        QExp = self.getExpectations(expand)
        return QExp['E']

    def updateParameters(self):

        # TODO: muZ node not accounted for here
        # Collect expectations from other node
        tmp = self.markov_blanket["Z"].getExpectations()
        if 'ENN' in tmp:
            EZZ = tmp["ENN"]
        else:
            EZZ = tmp["E2"]

        # Collect parameters from the P and Q distributions of this node
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb = P['a'], P['b']

        # Perform update
        Qa = Pa + 0.5*EZZ.shape[0]
        Qb = Pb + 0.5*EZZ.sum(axis=0)

        # Save updated parameters of the Q distribution
        self.Q.setParameters(a=Qa, b=Qb)

    def calculateELBO(self):
        # Collect parameters and expectations
        P,Q = self.P.getParameters(), self.Q.getParameters()
        Pa, Pb, Qa, Qb = P['a'], P['b'], Q['a'], Q['b']
        QE, QlnE = self.Q.getExpectations()['E'], self.Q.getExpectations()['lnE']

        # Do the calculations
        lb_p = (Pa*s.log(Pb)).sum() - special.gammaln(Pa).sum() + ((Pa-1.)*QlnE).sum() - (Pb*QE).sum()
        lb_q = (Qa*s.log(Qb)).sum() - special.gammaln(Qa).sum() + ((Qa-1.)*QlnE).sum() - (Qb*QE).sum()

        return lb_p - lb_q
