import scipy as s
import numpy as np
import scipy.stats as stats
import sys
import scipy.spatial as SS
import math
import seaborn as sns
import matplotlib.pyplot as plt

def simulate_data(N=200, seed=1234567, views = ["0", "1", "2", "3"], D = [500, 200, 500, 200], noise_level = 1,
                  K = 4, G = 1, lscales = [0.2, 0.8, 0.0, 0.0], sample_cov = "equidistant", scales = [1, 0.8, 0, 0],
                  model_groups = False, plot = False):
    """
    Function to simulate test data for MOFA (without ARD or spike-and-slab on factors)

    N: Number of time points/ samples per group
    seed: seed to use for simulation
    views: list of view names
    K: number of factors
    G: Number of groups
    D: list of number of features per view (same length as views)
    noise_level: variance of the residuals (1/tau);
                 per feature it is multiplied by a uniform random number in [0.5, 1.5] to model differences in features' noise
    scales, lscales: hyperparameters of the GP per factor (length as given by K)
    sample_cov: sample_covariates to use (can be of shape N X C) or "equidistant" or None
    model_groups: If False group covariance is assumed to be 1 across all groups, else it is filled with 0,1s at random
    plot: If True, simulation results are plotted
    """

    # simulate some test data
    np.random.seed(seed)
    M = len(views)
    N = int(N)
    groupidx = np.repeat(range(G), N)
    if not sample_cov is None:
        if sample_cov == "equidistant":
            sample_cov = np.linspace(0,1,N)
            sample_cov = sample_cov.reshape(N, 1)
        else:
            assert sample_cov.shape[0] == N, "sample_cov and N does not match"
            if len(np.repeat(np.arange(0,100,1),2).shape) == 1:
                sample_cov = sample_cov.reshape(N, 1)

        # expand across groups
        if not model_groups:
            sample_cov = np.tile(sample_cov, (G,1))
        distC = SS.distance.pdist(sample_cov, 'euclidean')**2.
        distC = SS.distance.squareform(distC)

    else:
        lscales = [0]* K

    if model_groups:
        Gmat = np.random.binomial(1, 0.5, G * K * G).reshape(K, G ,G)
        for k in range(K):
            Gmat[k,:,:] = np.dot(Gmat[k,:,:], Gmat[k,:,:].transpose())
            np.fill_diagonal(Gmat[k,:,:], 1)

    # simulate Sigma
    Sigma =[]
    for k in range(K):
        if lscales[k] > 0:
            Kmat = scales[k] * np.exp(-distC / (2 * lscales[k] ** 2))
            if model_groups:
                Kmat = np.kron(Gmat[k,:,:], Kmat)
            Sigma.append( Kmat + (1-scales[k]) * np.eye(N*G))
        elif lscales[k] == 0:
            Kmat = scales[k] * (distC == 0).astype(float)
            if model_groups:
                Kmat = np.kron(Kmat,  Gmat[k,:,:])
            Sigma.append(Kmat + (1-scales[k]) * np.eye(N*G))
            # Sigma.append(np.eye(N*G))
        else:
            sys.exit("All lengthscales need to be non-negative")

    # plot covariance structure
    if plot:
        fig, axs = plt.subplots(1, K, sharex=True, sharey=True)
        for k in range(K):
            sns.heatmap(Sigma[k], ax =axs[k])

    # simulate factor values
    Zks = []
    for k in range(K):
        sig = Sigma[k]
        Zks.append(np.random.multivariate_normal(np.zeros(N * G), sig, 1))
    Zks = np.vstack(Zks).transpose()

    Z = []
    for g in range(G):
        Z.append(Zks[groupidx == g,])

    # simulate alpha and theta, each factor should be active in at least one view
    inactive = 1000
    active = 1
    theta = 0.5 * np.ones([M, K])
    alpha_tmp = [s.ones(M) * inactive]*K
    for k in range(K):
        while s.all(alpha_tmp[k]==inactive):
            alpha_tmp[k] = s.random.choice([active,inactive], size=M, replace=True)
    alpha = [ s.array(alpha_tmp)[:,m] for m in range(M) ]

    # simulate weights
    W = []
    for m in range(M):
        W.append(np.column_stack(
            [np.random.normal(0, np.sqrt(1/alpha[m][k]), D[m]) * np.random.binomial(1, theta[m][k], D[m]) for k in range(K)]))

    # simulate heteroscedastic noise
    noise = []
    for m in range(M):
        tau_m = stats.uniform.rvs(loc=0.5, scale=1, size=D[m]) * 1/noise_level # uniform between 0.5 and 1.5 scaled by noise level
        noise.append(np.random.multivariate_normal(np.zeros(D[m]), np.eye(D[m]) * 1 / np.lib.scimath.sqrt(tau_m), N))

    # generate data
    data = []
    for m in range(M):
        tmp = []
        for g in range(G):
            tmp.append(Z[g].dot(W[m].transpose()) + noise[m])
        data.append(tmp)

    # store as list of groups
    if not sample_cov is None:
        if not model_groups:
            sample_cov = [sample_cov[groupidx == g] for g in range(G)]
        else:
            sample_cov = [sample_cov] * G

    if not model_groups:
        return {'data': data, 'W': W, 'Z': Z, 'noise': noise, 'sample_cov': sample_cov, 'Sigma': Sigma,
                'views': views, 'lscales': lscales, 'N': N}
    else:
        return {'data': data, 'W': W, 'Z': Z, 'noise': noise, 'sample_cov': sample_cov, 'Sigma': Sigma,
                'views': views, 'lscales': lscales, 'N': N, 'Gmat' : Gmat }



def mask_samples(sim, perc = 0.2, perc_all_views = 0):
    """
    Function to mask values at randomly sampled time points in each group and view.

    Param:
    perc: this fraction of time points are drawn in each view and group independently and all feature values at this time point set to NaN
    perc_all_views: this fraction of time points are drawn in each group independently and all feature values in all views are set to NaN
    """

    data = sim['data']
    N = sim['N']
    M = len(sim['views'])
    G = len(sim['data'][0])
    masked_samples = [[np.random.choice(N, math.floor(N * perc), replace = False) for g in range(G)] for m in range(M)]
    for m in range(M):
        for g in range(G):
            data[m][g][masked_samples[m][g],:] = s.nan

    if perc_all_views > 0:
        masked_samples_all_views = [np.random.choice(N, math.floor(N * perc_all_views), replace = False) for g in range(G)]
        for m in range(len(data)):
            for g in range(G):
                data[m][g][masked_samples_all_views[g], :] = s.nan


    return data
