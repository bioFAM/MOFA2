
#######################################################
## Functions to prepare a MOFA object for training ##
#######################################################

#' @title prepare a MOFA for training
#' @name prepare_mofa
#' @description Function to prepare a \code{\link{MOFA}} object for training. 
#' It requires defining data, model and training options.
#' @param object an untrained \code{\link{MOFA}}
#' @param data_options list of data_options (see \code{\link{get_default_data_options}} details). 
#' If NULL, default options are used.
#' @param model_options list of model options (see \code{\link{get_default_model_options}} for details). 
#' If NULL, default options are used.
#' @param training_options list of training options (see \code{\link{get_default_training_options}} for details). 
#' If NULL, default options are used.
#' @param stochastic_options list of options for stochastic variational inference (see \code{\link{get_default_stochastic_options}} for details). 
#' If NULL, default options are used.
#' @return Returns an untrained \code{\link{MOFA}} with specified options filled in the corresponding slots
#' @export
prepare_mofa <- function(object, data_options = NULL, model_options = NULL, training_options = NULL, stochastic_options = NULL,
                           regress_covariates = NULL ) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (any(object@dimensions$N<25)) warning("It is not a hard requirement, but with less than ~25 samples (per group) it is hard to learn meaningful factors...")
  if (any(object@dimensions$D<25)) warning("Some views have less than 25 features, they might be underrepresented...")
  
  # Get data options
  message("Checking data options...")
  if (is.null(data_options)) {
    message("No data options specified, using default...")
    object@data_options <- get_default_data_options(object)
  } else {
    if (!is(data_options,"list") || !setequal(names(data_options), names(get_default_data_options(object)) ))
      stop("data_options are incorrectly specified, please read the documentation in get_default_data_options")
    object@data_options <- data_options
  }
  if (any(nchar(unlist(samples_names(object)))>50))
    warning("Due to string size limitations in the HDF5 format, sample names will be trimmed to less than 50 characters")
  
  # Get training options
  if (is.null(training_options)) {
    message("No training options specified, using default...")
    object@training_options <- get_default_training_options(object)
  } else {
    message("Checking training options...")
    if (!is(training_options,"list") || !setequal(names(training_options), names(get_default_training_options(object)) ))
      stop("training_options are incorrectly specified, please read the documentation in get_default_training_options")
    object@training_options <- training_options
    
    # if (is.na(object@training_options$drop_factor_threshold))
    #   object@training_options$drop_factor_threshold <- -1
    # if (object@training_options$drop_factor_threshold>0.1)
    #   warning("Fraction of variance explained to drop factors is very high...")
    if (object@training_options$maxiter<=100)
      warning("Maximum number of iterations is very small\n")
    if (object@training_options$startELBO<1) object@training_options$startELBO <- 1
    if (object@training_options$freqELBO<1) object@training_options$freqELBO <- 1
    if (!object@training_options$convergence_mode %in% c("fast","medium","slow")) 
      stop("Convergence mode has to be either 'fast', 'medium', or 'slow'")
  }
  
  # Get stochastic options
  if (object@training_options$stochastic) {
    if (sum(object@dimensions$N) < 1e4) warning("Stochastic inference is only recommended when you have a lot of samples (at least N>10,000)")
      
    if (is.null(stochastic_options)) {
      message("No stochastic options specified, using default...")
      object@stochastic_options <- get_default_stochastic_options(object)
    } else {
        message("Checking stochastic inference options...")
      if(!is(stochastic_options,"list") || !setequal(names(stochastic_options), names(get_default_stochastic_options(object)) ))
        stop("stochastic_options are incorrectly specified, please read the documentation in get_default_stochastic_options")
      
      if (!stochastic_options %in% c(0.05,0.10,0.15,0.20,0.25,0.50))
        stop("Batch size has to be one of the following values: 0.05, 0.10, 0.15, 0.20, 0.25, 0.50")
      if (stochastic_options$batch_size==1)
        warning("A batch size equal to 1 is equivalent to non-stochastic inference. Please set object@train_options$stochastic <- FALSE")
      if (stochastic_options$learning_rate<=0 || stochastic_options$learning_rate>1)
        stop("The learning rate has to be a value between 0 and 1")
      if (stochastic_options$forgetting_rate<=0 || stochastic_options$forgetting_rate>1)
        stop("The forgetting rate has to be a value between 0 and 1")
      
      object@stochastic_options <- stochastic_options
    }
  }
  
  # Get model options
  if (is.null(model_options)) {
    message("No model options specified, using default...")
    object@model_options <- get_default_model_options(object)
  } else {
    message("Checking model options...")
    if (!is(model_options,"list") || !setequal(names(model_options), names(get_default_model_options(object)) ))
      stop("model_options are incorrectly specified, please read the documentation in get_default_model_options")
    object@model_options <- model_options
  }
  if (object@model_options$num_factors > 50) warning("The number of factors is very large, training will be slow...")
  if (!object@model_options$ard_weights) warning("model_options$ard_weights should always be set to TRUE")
  
  # Center the data
  # message("Centering the features (per group, this is a mandatory requirement)...")
  # for (m in views_names(object)) {
  #   if (model_options$likelihoods[[m]] == "gaussian") {
  #     for (g in groups_names(object)) {
  #       object@data[[m]][[g]] <- scale(object@data[[m]][[g]], center=T, scale=F)
  #     }
  #   }
  # }
  
  # Regress out covariates
  # object <- .regress_covariates(object, covariates)

  # Transform sparse matrices into dense ones
  # See https://github.com/rstudio/reticulate/issues/72
  for (m in views_names(object)) {
    for (g in groups_names(object)) {
      if (is(object@data[[m]][[g]], "dgCMatrix") || is(object@data[[m]][[g]], "dgTMatrix"))
        object@data[[m]][[g]] <- as(object@data[[m]][[g]], "matrix")
    }
  }
  
  return(object)
}



#' @title Get default training options
#' @name get_default_training_options
#' @description Function to obtain the default training options.
#' @param object an untrained \code{\link{MOFA}}
#' @details The training options are the following: \cr
#' \itemize{
#'  \item{\strong{maxiter}:}{ numeric value indicating the maximum number of iterations. 
#'  Default is 5000 (a lot). Convergence is assessed using the ELBO statistic.}
#'  \item{\strong{drop_factor_threshold}:}{ (not functional yet) numeric indicating the threshold on fraction of variance explained to consider a factor inactive and drop it from the model.
#'  For example, a value of 0.01 implies that factors explaining less than 1\% of variance (in each view) will be dropped.}
#'  \item{\strong{convergence_mode}:}{ character indicating the convergence criteria, either "slow", "medium" or "fast".}
#'  \item{\strong{verbose}:}{ logical indicating whether to generate a verbose output.}
#'  \item{\strong{startELBO}:}{ integer indicating the first iteration to compute the ELBO}
#'  \item{\strong{freqELBO}:}{ integer indicating the first iteration to compute the ELBO}
#'  \item{\strong{stochastic}:}{ logical indicating whether to use stochastic variational inference (only required for very big data sets).}
#'  \item{\strong{gpu_mode}:}{ logical indicating whether to use GPUs (see details).}
#'  \item{\strong{seed}:}{ numeric indicating the seed for reproducibility (default is 0, random seed).}
#' }
#' @return Returns a list with default training options
#' @importFrom utils modifyList
#' @export
get_default_training_options <- function(object) {
  
  # Get default train options
  training_options <- list(
    maxiter = 5000,                # (numeric) Maximum number of iterations
    convergence_mode = 'medium',   # (string) Convergence mode based on change in the ELBO ("slow","medium","fast")
    # drop_factor_threshold = -1,    # (numeric) Threshold on fraction of variance explained to drop a factor
    verbose = FALSE,               # (logical) verbosity
    startELBO = 1,                 # First iteration to compute the ELBO
    freqELBO = 5,                  # Frequency of ELBO calculation
    stochastic = FALSE,            # (logical) Do stochastic variational inference?
    gpu_mode = FALSE,              # (logical) Use GPU?
    seed = 0                       # (numeric) random seed
  )
  
  # if training_options already exist, replace the default values but keep the additional ones
  if (length(object@training_options)>0)
    training_options <- modifyList(training_options, object@training_options)
  
  return(training_options)
}


#' @title Get default data options
#' @name get_default_data_options
#' @description Function to obtain the default data options.
#' @param object an untrained \code{\link{MOFA}} object
#' @details The data options are the following: \cr
#' \itemize{
#'  \item{\strong{scale_views}:}{ logical indicating whether to scale views to have the same unit variance. 
#'  As long as the scale differences between the data sets is not too high, this is not required. Default is FALSE.}
#' }
#' @return Returns a list with the default data options.
#' @importFrom utils modifyList
#' @export
get_default_data_options <- function(object) {
  
  # Define default data options
  data_options <- list(
    scale_views = FALSE                # (logical) Scale views to unit variance
  )
  
  # if data_options already exists, replace the default values but keep the additional ones
  if (length(object@data_options)>0)
    data_options <- modifyList(data_options, object@data_options)
  
  return(data_options)
}

#' @title Get default model options
#' @name get_default_model_options
#' @description Function to obtain the default model options.
#' @param object an untrained \code{\link{MOFA}} object
#' @details The model options are the following: \cr
#' \itemize{
#'  \item{\strong{likelihoods}:}{ character vector with data likelihoods per view: 
#'  'gaussian' for continuous data, 'bernoulli' for binary data and 'poisson' for count data.
#'  By default, they are guessed internally.}
#'  \item{\strong{num_factors}:}{ numeric value indicating the (initial) number of factors. Default is 15.}
#'  \item{\strong{spikeslab_factors}:}{ logical indicating whether to use spike and slab sparsity on the factors (Default is FALSE)}
#'  \item{\strong{spikeslab_weights}:}{ logical indicating whether to use spike and slab sparsity on the weights (Default is TRUE)}
#'  \item{\strong{ard_factors}:}{ logical indicating whether to use ARD sparsity on the factors (Default is TRUE only if using multiple groups)}
#'  \item{\strong{ard_weights}:}{ logical indicating whether to use ARD sparsity on the weights (Default is TRUE)}
#'  }
#' @return Returns a list with the default model options.
#' @importFrom utils modifyList
#' @export
get_default_model_options <- function(object) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (!.hasSlot(object,"dimensions") | length(object@dimensions) == 0) 
    stop("dimensions of object need to be defined before getting the model options")
  if (.hasSlot(object,"data")) {
    if (length(object@data)==0) stop("data slot is empty")
  } else {
    stop("data slot not found")
  }
  
  # Guess likelihoods from the data
  likelihoods <- .infer_likelihoods(object)
  # likelihoods <- rep(x="gaussian", times=object@dimensions$M)
  names(likelihoods) <- views_names(object)
  
  # Define default model options
  model_options <- list(
    likelihoods = likelihoods,    # (character vector) likelihood per view [gaussian/bernoulli/poisson]
    num_factors = 15,            # (numeric) initial number of latent factors
    spikeslab_factors = FALSE,         # Spike and Slab sparsity on the factors
    spikeslab_weights = TRUE,          # Spike and Slab sparsity on the loadins
    ard_factors = TRUE,              # Group-wise ARD sparsity on the factors
    ard_weights = TRUE                # Group-wise ARD sparsity on the loadings
  )
  
  # Group-wise ARD sparsity on the factors only if there are multiple groups
  if (object@dimensions$G==1)
    model_options$ard_factors <-FALSE

  # if model_options already exist, replace the default values but keep the additional ones
  if (length(object@model_options)>0)
    model_options <- modifyList(model_options, object@model_options)
  
  return(model_options)
}





.regress_covariates <- function(object, covariates, min_observations = 10) {

  # First round of sanity checks
  if (!is(object, "MOFA")) 
    stop("'object' has to be an instance of MOFA")
  if (length(object@data)==0)
    stop("Input data has not been provided")
  
  # Fetch data
  views <- names(covariates)
  groups <- names(covariates[[1]])
  Y <- sapply(object@data[views], function(x) x[groups], simplify=F, USE.NAMES=T)
  # Y <- get_data(object, views=views, groups=groups)
  
  # Second round of sanity checks
  if (any(object@model_options$likelihoods[views]!="gaussian")) 
    stop("Some of the specified views contains discrete data. \nRegressing out covariates only works in views with continuous (gaussian) data")
  
  # Prepare data.frame with covariates
  if (!is(covariates,"list"))
    stop("Covariates has to be a list of vectors (for one covariate) or a list of data.frames (for multiple covariates")
  for (m in names(covariates)) {
    for (g in names(covariates[[m]])) {
      if (!is(covariates[[m]][[g]],"data.frame"))
        covariates[[m]][[g]] <- data.frame(x=covariates[[m]][[g]])
      stopifnot(nrow(covariates[[m]][[g]])==object@dimensions$N[g])
    }
  }
  
  print("Regressing out the specified covariates...")
  
  Y_regressed <- list()
  for (m in views) {
    Y_regressed[[m]] <- list()
    for (g in groups) {
      if (any(rowSums(!is.na(Y[[m]][[g]]))<min_observations) ) stop(sprintf("Some features do not have enough observations (N=%s) to fit the linear model",min_observations))
      Y_regressed[[m]][[g]] <- apply(Y[[m]][[g]], 2, function(y) {
        
        # Fit linear model
        df <- cbind(y,covariates[[m]][[g]])
        lm.out <- lm(y~., data=df)
        residuals <- lm.out[["residuals"]]
        
        # Fill missing values
        all_samples <- rownames(Y[[m]][[g]])
        missing_samples <- all_samples[!all_samples %in% names(residuals)]
        residuals[missing_samples] <- NA
        residuals[all_samples]
      })
    }
  }
  object@data[views] <- Y_regressed
  
  return(object)
}

#' @title Get default stochastic options
#' @name get_default_stochastic_options
#' @description Function to obtain the default options for stochastic variational inference.
#' @param object an untrained \code{\link{MOFA}}
#' @details The training options are the following: \cr
#' \itemize{
#'  \item{\strong{batch_size}:}{ numeric value indicating the batch size (as a fraction)}. 
#'  Default is 0.5 (half of the data set).
#'  \item{\strong{learning_rate}:}{ numeric value indicating the learning rate. }
#'  Default is 0.75 
#'  \item{\strong{forgetting_rate}:}{ numeric indicating the forgetting rate.}
#'  Default is 1.0
#'  }
#' @return Returns a list with default options
#' @importFrom utils modifyList
#' @export
get_default_stochastic_options <- function(object) {
  
  # Get default stochastic options
  stochastic_options <- list(
    batch_size = 0.5,       # Batch size (as a fraction)
    learning_rate = 0.75,   # Starting learning rate
    forgetting_rate = 0.1   # Forgetting rate
  )
  
  # if stochastic_options already exist, replace the default values but keep the additional ones
  if (length(object@stochastic_options)>0)
    stochastic_options <- modifyList(stochastic_options, object@stochastic_options)
  
  return(stochastic_options)
}

