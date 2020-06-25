import scipy as s
import numpy as np
import scipy.stats as stats
import sys
import scipy.spatial as SS
import math

def simulate_data(N=200, seed=1234567, views = ["0", "1", "2", "3"], D = [500, 200, 500, 200], noise_level = 1,
                  K = 4, G = 1, lscales = [0.2, 0.8, 0.0, 0.0], sample_cov = "equidistant"):
    """
    Function to simulate test data for MOFA (without ARD or spike-and-slab on factors)
    """
    # simulate some test data
    np.random.seed(seed)
    M = len(views)
    N = int(N)
    if not sample_cov is None:
        if sample_cov == "equidistant":
            sample_cov = np.linspace(0,1,N)
            # sample_cov = np.array(range(N))
            sample_cov = sample_cov.reshape(N, 1)
        else:
            assert sample_cov.shape[0] == N, "sample_cov and N does not match"
            if len(np.repeat(np.arange(0,100,1),2).shape) == 1:
                sample_cov = sample_cov.reshape(N, 1)

        distC = SS.distance.pdist(sample_cov, 'euclidean')**2.
        distC = SS.distance.squareform(distC)

    else:
        lscales = [0]* K

    # simulate Sigma
    Sigma =[]
    for k in range(K):
        if lscales[k] > 0:
            Sigma.append(np.exp(-distC / (2 * lscales[k] ** 2)))
        elif lscales[k] == 0:
            Sigma.append(np.eye(N))
        else:
            sys.exit("All lengthscales need to be positive")

    # simulate same factor values for factors with non-zero lengthscales across groups, other can differ (avoids expanding the covaraicen matrix across groups)
    Z = []
    for g in range(G):
        Zks = []
        for k in range(K):
            sig = Sigma[k]
            if lscales[k] == 0:
                Zks.append(np.random.multivariate_normal(np.zeros(N), sig, 1))
            else:
                if g == 0:
                    Zks.append(np.random.multivariate_normal(np.zeros(N), sig, 1))
                else:
                    Zks.append(Z[0][k])
        Z.append(Zks)
    Z = [np.vstack(zlist).transpose() for zlist in Z]

    W = []

    # simulate alpha and theta, each factor should be active in at least one view
    inactive = 1000
    active = 1
    theta = 0.5 * np.ones([M, K])
    alpha_tmp = [s.ones(M)*inactive]*K
    for k in range(K):
        while s.all(alpha_tmp[k]==inactive):
            alpha_tmp[k] = s.random.choice([active,inactive], size=M, replace=True)
    alpha = [ s.array(alpha_tmp)[:,m] for m in range(M) ]

    for m in range(M):
        W.append(np.column_stack(
            [np.random.normal(0, np.sqrt(1/alpha[m][k]), D[m]) * np.random.binomial(1, theta[m][k], D[m]) for k in range(K)]))

    # simulate heteroscedastic noise
    noise = []
    for m in range(M):
        tau_m = stats.uniform.rvs(loc=1, scale=3, size=D[m]) * noise_level
        noise.append(np.random.multivariate_normal(np.zeros(D[m]), np.eye(D[m]) * 1 / np.lib.scimath.sqrt(tau_m), N))

    data = []
    for m in range(M):
        tmp = []
        for g in range(G):
            tmp.append(Z[g].dot(W[m].transpose()) + noise[m])
        data.append(tmp)

    # store as list of groups
    if not sample_cov is None:
        sample_cov = [sample_cov] * G
    return {'data': data, 'W': W, 'Z': Z, 'noise': noise, 'sample_cov': sample_cov, 'Sigma': Sigma,
            'views': views, 'lscales': lscales, 'N': N}

# mask values
def mask_data(sim, perc = 0.2):
    data = sim['data']
    N = sim['N']
    M = len(sim['views'])
    masked_samples = [np.random.choice(N, math.floor(N * perc), replace = False) for m in range(M)]

    for m in range(len(data)):
        data[m][0][masked_samples[m],:] = s.nan

    return data
