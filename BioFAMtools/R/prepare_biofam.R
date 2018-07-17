
####################################################
## Functions to prepare a BioFAM object for training ##
####################################################

#' @title prepare a BioFAM object for training
#' @name prepare_biofam
#' @description Function to prepare a \code{\link{BioFAModel}} object for training.
#' Here, data, input/output option are specified and data, model and training options can be set.
#' @param object an untrained \code{\link{BioFAModel}}
#' @param data_options list of data_options (see \code{\link{get_default_data_options}} details). 
#' If NULL, default data options are used.
#' @param model_options list of model_options (see \code{\link{get_default_model_options}} for details). 
#' If NULL, default model options are used.
#' @param training_options list of training_options (see \code{\link{get_default_training_options}} for details). 
#' If NULL, default training options are used.
#' @return Returns an untrained \code{\link{BioFAModel}} with specified data, model and training options.
#' Next step is to train the model with \code{\link{run_biofam}}
#' @export


prepare_biofam <- function(object, data_options = NULL, model_options = NULL, training_options = NULL) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) 
    stop("'object' has to be an instance of BioFAModel")
  
  # Get data options
  message("Checking data options...")
  if (is.null(data_options)) {
    message("No data options specified, using default...")
    object@data_options <- get_default_data_options()
  } else {
    if (!is(training_options, "list") & !all(names(training_options) == names(get_default_training_options())))
      stop("data_options are incorrectly specified, please read the documentation in get_default_data_options")
    object@data_options <- data_options
  }
  
  # Get training options
  message("Checking training options...")
  if (is.null(training_options)) {
    message("No training options specified, using default...")
    object@training_options <- get_default_training_options()
  } else {
    if(!is(training_options,"list") & !all(names(training_options) == names(get_default_training_options())))
      stop("training_options are incorrectly specified, please read the documentation in get_default_training_options")
    object@training_options <- training_options
  }
  
  # Get model options
  message("Checking model options...")
  if(is.null(model_options)) {
    message("No model options specified, using default...")
    object@model_options <- get_default_model_options(object)
  } else {
    # (To-do) Check that model_options is correct
    if(!is(model_options,"list") & !all(names(model_options) == names(get_default_model_options(object))))
      stop("model_options are incorrectly specified, please read the documentation in get_default_model_options")
    object@model_options <- model_options
  }
  
  return(object)
}



#' @title Get default training options
#' @name get_default_training_options
#' @description Function to obtain the default training options.
#' @details The training options are the following: \cr
#' \itemize{
#'  \item{\strong{maxiter}:}{ numeric value indicating the maximum number of iterations. 
#'  Default is 5000, but we recommend using the 'tolerance' as convergence criteria.}
#'  \item{\strong{tolerance}:}{ numeric value indicating the convergence threshold based
#'   on the change in Evidence Lower Bound (deltaELBO). 
#'  For quick exploration we recommend this to be around 1.0,
#'   and for a thorough training we recommend a value of 0.01. Default is 0.1}
#'  \item{\strong{DropFactorThreshold}:}{ numeric hyperparamter to automatically learn the number of factors.
#'  It indicates the threshold on fraction of variance explained to consider a factor inactive and 
#'  automatically drop it from the model during training. 
#'  For example, a value of 0.01 implies that factors explaining less
#'  than 1\% of variance (in each view) will be dropped.
#'  Default is 0, which implies that only factors that explain no variance at all will be removed
#'  }
#'  \item{\strong{verbose}:}{ logical indicating whether to generate a verbose output.}
#'  \item{\strong{seed}:}{ random seed for reproducibility (default is NULL, which samples a random seed).}
#' }
#' @return Returns a list with default training options, which have to be passed
#'as an argument to \code{\link{prepare_biofam}}
#' @export
#' @examples 
#' training_options <- get_default_training_options()
#' training_options

get_default_training_options <- function() {
  training_options <- list(
    maxiter = 5000,               # (numeric) Maximum number of iterations
    tolerance = 0.1,              # (numeric) Convergence threshold based on change in the evidence lower bound
    drop_factor_threshold = 0.00, # (numeric) Threshold on fraction of variance explained to drop a factor
    verbose = FALSE,              # (logical) verbosity?
    seed = NULL                   # (numeric or NULL) random seed
  )
  return(training_options)
}


#' @title Get default data options
#' @name get_default_data_options
#' @description Function to obtain the default data options.
#' @details The data options are the following: \cr
#' \itemize{
#'  \item{\strong{centerFeatures}:}{ logical indicating whether to center the features to zero mean. 
#'  This only works for gaussian data. Default is TRUE.}
#'  \item{\strong{scaleViews}:}{ logical indicating whether to scale views to have the same unit variance. 
#'  As long as the scale differences between the data sets is not too high, this is not required.
#'   Default is FALSE.}
#'  \item{\strong{removeIncompleteSamples}:}{ logical indicating whether to remove samples that
#'   are not profiled in all omics. We recommend this only for testing,
#'    as the model can cope with samples having missing assays. Default is FALSE.}
#' }
#' @return Returns a list with the default data options, which have to be passed as
#'  an argument to \code{\link{prepare_biofam}}
#' @export
#' @examples 
#' data_options <- get_default_data_options()
#' data_options

get_default_data_options <- function() {
  data_options <- list(
    center_features = TRUE,           # Center features to zero mean (does not apply to binary or count views)
    scale_views = FALSE,              # Scale views to unit variance (does not apply to binary or count views)
    remove_incomplete_samples = FALSE # Remove incomplete samples that are not profiled in all omics?
  )
  return(data_options)
}

#' @title Get default model options
#' @name get_default_model_optionss
#' @param object an untrained \code{\link{BioFAModel}} object
#' @description Function to obtain the default model options.
#' @details The model options are the following: \cr
#' \itemize{
#'  \item{\strong{likelihood}:}{ character vector with data likelihoods per view: 
#'  'gaussian' for continuous data, 'bernoulli' for binary data and 'poisson' for count data.
#'  By default, they are guessed internally.}
#'  \item{\strong{numFactors}:}{ numeric value indicating the initial number of factors. 
#'  If you want to learn the number of factors automatically we recommend
#'   setting this to a large value, default is 25.}
#'  \item{\strong{sparsity}:}{ logical indicating whether to use sparsity.
#'  This is always recommended, as it will make the loadings more interpretable. Default is TRUE.}
#'  \item{\strong{learnIntercept}:}{ logical indicating whether to learn an intercept
#'   term to capture differences in feature means. This is currently not functional. Therefore, we ask you to always center
#'    the data if it is gaussian (do not worry about the non-gaussian views) and keep learnIntercept to FALSE.}
#' }
#' @return Returns a list with the default model options, which have to be passed as
#'  an argument to \code{\link{prepare_biofam}}
#' @export

get_default_model_options <- function(object) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  if (!.hasSlot(object, "dimensions") | length(object@dimensions) == 0) 
    stop("dimensions of object need to be defined before getting model_optionss")
  if (!.hasSlot(object, "input_data")) 
    stop("input_data slot needs to be specified before getting model_optionss")
  if (!.hasSlot(object, "training_data")) 
    stop("training_data slot needs to be specified before getting model_optionss")
  
  # Guess likelihood type
  likelihood <- .infer_likelihoods(object)
  
  # Define default model options
  model_options <- list(
    likelihood = likelihood,    # (character vector) likelihood per view [gaussian/bernoulli/poisson]
    learn_intercept = FALSE,    # (bool) include a constant factor of 1s to learn the mean of features (intercept)?
    num_factors = 25,           # (numeric) initial number of latent factors
    sl_w = TRUE,                # (logical) use feature-wise sparsity?
    sl_z = FALSE,               # (logical) use sample-wise sparsity?
    ard_w = TRUE,               # (logical) use view-wise sparsity?
    ard_z = TRUE,               # (logical) use group-wise sparsity?
    noise_on = 'features'       # (character) noise on features or on samples
  )
  
  return(model_options)
}