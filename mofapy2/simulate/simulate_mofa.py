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
                  group_modelling = False, plot = False):
    """
    Function to simulate test data for MOFA (without ARD or spike-and-slab on factors)
    noise_level: variance of the residuals (1/tau); for each feature this is multiplied by a uniform random number between 0.5 and 1.5 to model differences in feature's noise
    scales, lscales: hyperparameters of the GP per factor
    G: Number of groups
    N: Number of time points/ samples per group
    """
    # # set scales to 0 for lengthscale of 0 to obtain identity
    # for k in range(K):
    #     if lscales[k] == 0:
    #         if scales[k] != 0:
    #             scales[k] = 0
    #             print("Setting scale of factor ", k, " to 0 as it has a lengthscale of 0 ")

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
        if not group_modelling:
            sample_cov = np.tile(sample_cov, (G,1))
        distC = SS.distance.pdist(sample_cov, 'euclidean')**2.
        distC = SS.distance.squareform(distC)

    else:
        lscales = [0]* K

    if group_modelling:
        rank_x = 2
        sigma = 0
        # x = np.random.binomial(1, 0.5, G * K * rank_x).reshape(K, rank_x ,G)
        # Gmat = np.zeros([K, G, G])
        # for k in range(K):
        #     Gmat[k,:,:] = np.dot(x[k,:,:].transpose(),x[k,:,:]) + sigma * np.eye(G)
        Gmat = np.random.binomial(1, 0.5, G * K * G).reshape(K, G ,G)
        for k in range(K):
            Gmat[k,:,:] = np.dot(Gmat[k,:,:], Gmat[k,:,:].transpose())
            np.fill_diagonal(Gmat[k,:,:], 1)

    # simulate Sigma
    Sigma =[]
    for k in range(K):
        if lscales[k] > 0:
            Kmat = scales[k] * np.exp(-distC / (2 * lscales[k] ** 2))
            if group_modelling:
                Kmat = np.kron(Gmat[k,:,:], Kmat)
            Sigma.append( Kmat + (1-scales[k]) * np.eye(N*G))
        elif lscales[k] == 0:
            Kmat = scales[k] * (distC == 0).astype(float)
            if group_modelling:
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

    # simulate same factor values for factors with non-zero lengthscales across groups, other can differ (avoids expanding the covaraicen matrix across groups)
    Zks = []
    for k in range(K):
        sig = Sigma[k]
        Zks.append(np.random.multivariate_normal(np.zeros(N * G), sig, 1))
    Zks = np.vstack(Zks).transpose()

    Z = []
    for g in range(G):
        Z.append(Zks[groupidx == g,])

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
        tau_m = stats.uniform.rvs(loc=0.5, scale=1, size=D[m]) * 1/noise_level # uniform between 0.5 and 1.5 scaled by noise level
        noise.append(np.random.multivariate_normal(np.zeros(D[m]), np.eye(D[m]) * 1 / np.lib.scimath.sqrt(tau_m), N))

    data = []
    for m in range(M):
        tmp = []
        for g in range(G):
            tmp.append(Z[g].dot(W[m].transpose()) + noise[m])
        data.append(tmp)

    # store as list of groups
    if not sample_cov is None:
        if not group_modelling:
            sample_cov = [sample_cov[groupidx == g] for g in range(G)]
        else:
            sample_cov = [sample_cov] * G

    if not group_modelling:
        return {'data': data, 'W': W, 'Z': Z, 'noise': noise, 'sample_cov': sample_cov, 'Sigma': Sigma,
                'views': views, 'lscales': lscales, 'N': N}
    else:
        return {'data': data, 'W': W, 'Z': Z, 'noise': noise, 'sample_cov': sample_cov, 'Sigma': Sigma,
                'views': views, 'lscales': lscales, 'N': N, 'Gmat' : Gmat }

# mask values at randomly sampled time points in each group and view
# perc of time points are drawn in each view and group independtly and all feature values at this time point set to NaN
# perc_all_views of time points are drawn in each group independently and all feature values in all views are set to NaN
def mask_samples(sim, perc = 0.2, perc_all_views = 0):
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
