
#' @title Simulate a data set using the generative model of MOFA
#' @name make_example_data
#' @description Function to simulate an example multi-view multi-group data set according to the generative model of MOFA2.
#' @param n_views number of views
#' @param n_features number of features in each view 
#' @param n_samples number of samples in each group
#' @param n_groups number of groups
#' @param n_factors number of factors
#' @param likelihood likelihood for each view, one of "gaussian" (default), "bernoulli", "poisson",
#'  or a character vector of length n_views
#' @param lscales vector of lengthscales, needs to be of length n_factors (default is 0 - no smooth factors)
#' @param sample_cov matrix of sample covariates for one group with covariates in rows and samples in columns 
#' or "equidistant" for sequential ordering, default is NULL (no smooth factors)
#' @param as.data.frame return data and covariates as long dataframe 
#' @return Returns a list containing the simulated data and simulation parameters.
#' @importFrom stats rnorm rbinom rpois
#' @importFrom dplyr left_join
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats dist
#' @export
#' @examples
#' # Generate a simulated data set
#' MOFAexample <- make_example_data()


make_example_data <- function(n_views=3, n_features=100, n_samples = 50, n_groups = 1,
                            n_factors = 5, likelihood = "gaussian",
                            lscales = 1, sample_cov = NULL, as.data.frame = FALSE) {
  
  # Sanity checks
  if (!all(likelihood %in% c("gaussian", "bernoulli", "poisson")))
    stop("Likelihood not implemented: Use either gaussian, bernoulli or poisson")
  
  if(length(lscales) == 1)
    lscales = rep(lscales, n_factors)
  if(!length(lscales) == n_factors)
    stop("Lengthscalces lscales need to be of length n_factors")
  if(all(lscales == 0)){
    sample_cov <- NULL
  }
  
  if (length(likelihood)==1) likelihood <- rep(likelihood, n_views) 
  if (!length(likelihood) == n_views) 
    stop("Likelihood needs to be a single string or matching the number of views!")
  
  if(!is.null(sample_cov)){
    if(sample_cov[1] == "equidistant") {
      sample_cov <- seq_len(n_samples)
    }
    if(is.null(dim(sample_cov))) sample_cov <- matrix(sample_cov, nrow = 1)
    if(ncol(sample_cov) != n_samples){
      stop("Number of columns in sample_cov must match number of samples n_samples.")
    }
  
    # Simulate covariance for factors
    Sigma = lapply(lscales, function(ls) {
      if(ls == 0) diag(1, n_samples)
      else (1) * exp(-as.matrix(stats::dist(t(sample_cov)))^2/(2*ls^2))
      # else (1-0.001) * exp(-as.matrix(stats::dist(t(sample_cov)))^2/(2*ls^2)) + diag(0.001, n_samples)
    })
  
    # simulate factors
    alpha_z <- NULL
    S_z <- lapply(seq_len(n_groups), function(vw) matrix(1, nrow=n_samples, ncol=n_factors))
    Z <-  vapply(seq_len(n_factors), function(fc) mvtnorm::rmvnorm(1, rep(0, n_samples), Sigma[[fc]]), numeric(n_samples))
    colnames(Z) <- paste0("simulated_factor_", 1:ncol(Z))
    Z <- lapply(seq_len(n_groups), function(gr) Z)
    sample_cov <- Reduce(cbind, lapply(seq_len(n_groups), function(gr) sample_cov))
  } else {
    # set sparsity for factors
    theta_z <- 0.5
    
    # set ARD prior for factors, each factor being active in at least one group
    alpha_z <- vapply(seq_len(n_factors), function(fc) {
      active_gw <- sample(seq_len(n_groups), 1)
      alpha_fc <- sample(c(1, 1000), n_groups, replace = TRUE)
      if(all(alpha_fc==1000)) alpha_fc[active_gw] <- 1
      alpha_fc
    }, numeric(n_groups))
    alpha_z <- matrix(alpha_z, nrow=n_factors, ncol=n_groups, byrow=TRUE)
    
    # simulate facors 
    S_z <- lapply(seq_len(n_groups), function(vw) matrix(rbinom(n_samples * n_factors, 1, theta_z),
                                                         nrow=n_samples, ncol=n_factors))
    Z <- lapply(seq_len(n_groups), function(vw) vapply(seq_len(n_factors), function(fc) rnorm(n_samples, 0, sqrt(1/alpha_z[fc,vw])), numeric(n_samples)))
  }
  
  # set sparsity for weights
  theta_w <- 0.5
  
  # set ARD prior, each factor being active in at least one view
  alpha_w <- vapply(seq_len(n_factors), function(fc) {
    active_vw <- sample(seq_len(n_views), 1)
    alpha_fc <- sample(c(1, 1000), n_views, replace = TRUE)
    if(all(alpha_fc==1000)) alpha_fc[active_vw] <- 1
    alpha_fc
  }, numeric(n_views))
  alpha_w <- matrix(alpha_w, nrow=n_factors, ncol=n_views, byrow=TRUE)
  
  # simulate weights 
  S_w <- lapply(seq_len(n_views), function(vw) matrix(rbinom(n_features*n_factors, 1, theta_w),
                                             nrow=n_features, ncol=n_factors))
  W <- lapply(seq_len(n_views), function(vw) vapply(seq_len(n_factors), function(fc) rnorm(n_features, 0, sqrt(1/alpha_w[fc,vw])), numeric(n_features)))
  
  # set noise level (for gaussian likelihood)
  tau <- 10
  
  # pre-compute linear term and rbind groups
  mu <- lapply(seq_len(n_views), function(vw) lapply(seq_len(n_groups), function(gw)  (S_z[[gw]]*Z[[gw]]) %*% t(S_w[[vw]]*W[[vw]])))
  mu <- lapply(mu, function(l) Reduce(rbind, l))
  groups <- rep(paste("group",seq_len(n_groups), sep = "_"), each = n_samples)
  
  # simulate data according to the likelihood
  data <- lapply(seq_len(n_views), function(vw){
    lk <- likelihood[vw]
    if (lk == "gaussian"){
      dd <- t(mu[[vw]] + rnorm(length(mu[[vw]]),0,sqrt(1/tau)))
    }
    else if (lk == "poisson"){
      term <- log(1+exp(mu[[vw]]))
      dd <- t(apply(term, 2, function(tt) rpois(length(tt),tt)))
    }
    else if (lk == "bernoulli") {
      term <- 1/(1+exp(-mu[[vw]]))
      dd <- t(apply(term, 2, function(tt) rbinom(length(tt),1,tt)))
    }
    colnames(dd) <- paste0("sample_", seq_len(ncol(dd)))
    rownames(dd) <- paste0("feature_", seq_len(nrow(dd)),"_view", vw)
    dd
  })

  if(!is.null(sample_cov)) {
    colnames(sample_cov) <- colnames(data[[1]])
    rownames(sample_cov) <- paste0("covariate_", seq_len(nrow(sample_cov)))
  }
 
  names(data) <- paste0("view_", seq_len(n_views))
  
  if(as.data.frame){
    gr_df <- data.frame(group = groups, sample = colnames(data[[1]]))
    dat <- lapply(names(data), function(vw){
            tmp <- data[[vw]]
            df <- melt(tmp, varnames = c("feature", "sample"))
            df$view <- vw
            df
    })
    data <- bind_rows(dat)
    data <- dplyr::left_join(data, gr_df, by = "sample")
    
    sample_cov <- melt(sample_cov, varnames = c("covariate", "sample"))
  }
  return(list(data = data, groups = groups, alpha_w=alpha_w, alpha_z =alpha_z,
              lscales = lscales, sample_cov = sample_cov, Z = Z))
}
