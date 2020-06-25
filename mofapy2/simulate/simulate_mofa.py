import scipy as s
import numpy as np
import scipy.stats as stats
import sys
import scipy.spatial as SS
import math

def simulate_data(N=200, seed=1234567, views = ["0", "1", "2", "3"], D = [500, 200, 500, 200], noise_level = 1,
                  K = 4, lscales = [0.2, 0.8, 0.0, 0.0], sample_cov = None):
    """
    Function to simulate test data for MOFA (with one group)
    """
    # simulate some test data
    np.random.seed(seed)
    M = len(views)
    N = int(N)
    if sample_cov is None:
        sample_cov = np.linspace(0,1,N)
        # sample_cov = np.array(range(N))
        sample_cov = sample_cov.reshape(N, 1)
    else:
        assert sample_cov.shape[0] == N, "sample_cov and N does not match"
        if len(np.repeat(np.arange(0,100,1),2).shape) == 1:
            sample_cov = sample_cov.reshape(N, 1)

    distC = SS.distance.pdist(sample_cov, 'euclidean')**2.
    distC = SS.distance.squareform(distC)

    Sigma =[]
    for k in range(K):
        if lscales[k] > 0:
            Sigma.append(np.exp(-distC / (2 * lscales[k] ** 2)))
        elif lscales[k] == 0:
            Sigma.append(np.eye(N))
        else:
            sys.exit("All lengthscales need to be positive")

    # Sigma = [np.exp(-np.array([[(i - j) ** 2 for i in sample_cov] for j in sample_cov]) / (2 * ls ** 2)) for
    #          ls in lscales]

    Z = np.vstack([np.random.multivariate_normal(np.zeros(N), sig, 1) for sig in
                   Sigma]).transpose()  # mix of spatial and non-spatial factors

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
        data.append(Z.dot(W[m].transpose()) + noise[m])
    data = [[data[m]] for m in range(M)]

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
