
#######################################################
## Functions to prepare a bioFAM object for training ##
#######################################################

#' @title prepare a BioFAModel for training
#' @name prepare_biofam
#' @description Function to prepare a \code{\link{BioFAModel}} object for training. 
#' It requires defining data, model and training options.
#' @param object an untrained \code{\link{BioFAModel}}
#' @param dir_options list with Input/Output options, which must contain a 'data_dir' element where temporary text files will be stored 
#' and an 'outFile' with hdf5 extension where the final model will be stored.
#' @param data_options list of data_options (see \code{\link{get_default_data_options}} details). 
#' If NULL, default data options are used.
#' @param model_options list of model options (see \code{\link{get_default_model_options}} for details). 
#' If NULL, default model options are used.
#' @param training_options list of training options (see \code{\link{get_default_training_options}} for details). 
#' If NULL, default training options are used.
#' @return Returns an untrained \code{\link{BioFAModel}} with specified data, model and training options
#' @export
prepare_biofam <- function(object, data_options = NULL, model_options = NULL, training_options = NULL,
                           center=TRUE, regress_covariates=NULL ) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  
  # Create temporary folder to store data
  # dir.create(dir_options$data_dir, showWarnings = FALSE)
  
  # Get data options
  message("Checking data options...")
  if (is.null(data_options)) {
    message("No data options specified, using default...")
    object@data_options <- get_default_data_options(object)
  } else {
    if (!is(data_options,"list") | !setequal(names(data_options), names(get_default_data_options(object)) ))
      stop("data_options are incorrectly specified, please read the documentation in get_default_data_options")
    object@data_options <- data_options
  }
  if (any(nchar(unlist(samples_names(object)))>50))
    warning("Due to string size limitations in the HDF5 format, sample names will be trimmed to less than 50 characters")
  
  # Get training options
  message("Checking training options...")
  if (is.null(training_options)) {
    message("No training options specified, using default...")
    object@training_options <- get_default_training_options(object)
  } else {
    if(!is(training_options,"list") | !setequal(names(training_options), names(get_default_training_options(object)) ))
      stop("training_options are incorrectly specified, please read the documentation in get_default_training_options")
    object@training_options <- training_options
  }
  
  # Get model options
  message("Checking model options...")
  if (is.null(model_options)) {
    message("No model options specified, using default...")
    object@model_options <- get_default_model_options(object)
  } else {
    if (!is(model_options,"list") | !setequal(names(model_options), names(get_default_model_options(object)) ))
      stop("model_options are incorrectly specified, please read the documentation in get_default_model_options")
    object@model_options <- model_options
  }
  
  # Center the data
  if (center) {
    print("Centering the data (per group)...")
    for (m in views_names(object)) {
      for (g in groups_names(object)) {
        object@input_data[[m]][[g]] <- scale(object@input_data[[m]][[g]], center=T, scale=F)
      }
    }
  }
  
  # Regress out covariates
  object <- .regress_covariates(object, covariates)
  
  return(object)
}



#' @title Get default training options
#' @name get_default_training_options
#' @description Function to obtain the default training options.
#' @param object an untrained \code{\link{BioFAModel}}
#' @details The training options are the following: \cr
#' \itemize{
#'  \item{\strong{maxiter}:}{ numeric value indicating the maximum number of iterations. 
#'  Default is 5000, but we recommend using the 'tolerance' as convergence criteria.}
#'  \item{\strong{tolerance}:}{ numeric value indicating the convergence threshold based on the change in Evidence Lower Bound (deltaELBO). 
#'  For quick exploration we recommend this to be around 1.0, and for a thorough training we recommend a value of 0.01. Default is 0.1}
#'  \item{\strong{drop_factor_threshold}:}{ numeric indicating the threshold on fraction of variance explained to consider a factor inactive and drop it from the model.
#'  For example, a value of 0.01 implies that factors explaining less than 1\% of variance (in each view) will be dropped.}
#'  \item{\strong{verbose}:}{ logical indicating whether to generate a verbose output.}
#'  \item{\strong{seed}:}{ random seed for reproducibility (default is 0, random seed).}
#' }
#' @return Returns a list with default training options
#' @export
get_default_training_options <- function(object) {
  
  # Get default train options
  training_options <- list(
    maxiter = 5000,                # (numeric) Maximum number of iterations
    tolerance = 0.1,               # (numeric) Convergence threshold based on change in the evidence lower bound
    drop_factor_threshold = 0.02,  # (numeric) Threshold on fraction of variance explained to drop a factor
    verbose = FALSE,               # (logical) verbosity?
    seed = 0                       # (numeric or NULL) random seed
  )
  
  # if training_options already exist, replace the default values but keep the additional ones
  if (length(object@training_options)>0)
    training_options <- modifyList(training_options, object@training_options)
  
  return(training_options)
}




#' @title Get default data options
#' @name get_default_data_options
#' @description Function to obtain the default data options.
#' @param object an untrained \code{\link{BioFAModel}} object
#' @details The data options are the following: \cr
#' \itemize{
#'  \item{\strong{scale_views}:}{ logical indicating whether to scale views to have the same unit variance. 
#'  As long as the scale differences between the data sets is not too high, this is not required. Default is FALSE.}
#'  \item{\strong{center_features}:}{ logical indicating whether to center the features to zero mean. This only works for gaussian data. Default is TRUE.}
#' }
#' @return Returns a list with the default data options, which have to be passed as an argument to \code{\link{prepareMOFA}}
#' @export
get_default_data_options <- function(object) {
  
  # Define default data options
  data_options <- list(
    center_features_per_group = FALSE, # (logical) Center features to zero mean, for each group separately
    scale_views = FALSE                # (logical) Scale views to unit variance
  )
  
  # if data_options already exist, replace the default values but keep the additional ones
  if (length(object@data_options)>0)
    data_options <- modifyList(data_options, object@data_options)
  
  return(data_options)
}

#' @title Get default model options
#' @name get_default_model_options
#' @description Function to obtain the default model options.
#' @param object an untrained \code{\link{BioFAModel}} object
#' @details The model options are the following: \cr
#' \itemize{
#'  \item{\strong{likelihood}:}{ character vector with data likelihoods per view: 
#'  'gaussian' for continuous data, 'bernoulli' for binary data and 'poisson' for count data.
#'  By default, they are guessed internally.}
#'  \item{\strong{num_factors}:}{ numeric value indicating the initial number of factors. 
#'  If you want to learn the number of factors automatically we recommend setting this to a large value, around 50. Default is 25.}
#'  }
#' @return Returns a list with the default model options, which have to be passed as an argument to \code{\link{prepareMOFA}}
#' @export
get_default_model_options <- function(object) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  if (!.hasSlot(object,"dimensions") | length(object@dimensions) == 0) 
    stop("dimensions of object need to be defined before getting the model options")
  if (.hasSlot(object,"input_data")) {
    if (length(object@input_data)==0) stop("input_data slot is empty")
  } else {
    stop("input_data slot not found")
  }
  
  # Guess likelihoods from the data
  # likelihood <- .infer_likelihoods(object) THIS DOES NOT WORK 
  likelihood <- rep(x="gaussian", times=object@dimensions$M)
  names(likelihood) <- views_names(object)
  
  # Define default model options
  model_options <- list(
    likelihood = likelihood,    # (character vector) likelihood per view [gaussian/bernoulli/poisson]
    num_factors = 25,            # (numeric) initial number of latent factors
    sl_z = TRUE,
    sl_w = TRUE,
    ard_w = TRUE,
    ard_z = TRUE
  )

  # if model_options already exist, replace the default values but keep the additional ones
  if (length(object@model_options)>0)
    model_options <- modifyList(model_options, object@model_options)
  
  return(model_options)
}







.regress_covariates <- function(object, covariates, min_observations = 10) {

  # First round of sanity checks
  if (!is(object, "BioFAModel")) 
    stop("'object' has to be an instance of BioFAModel")
  if (length(object@input_data)==0)
    stop("Input data has not been provided")
  
  # Fetch data
  views <- names(covariates)
  groups <- names(covariates[[1]])
  Y <- sapply(object@input_data[views], function(x) x[groups], simplify=F, USE.NAMES=T)
  # Y <- get_input_data(object, views=views, groups=groups)
  
  # Second round of sanity checks
  if (any(object@model_options$likelihood[views]!="gaussian")) 
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
  
  Y_regressed <- list()
  for (m in views) {
    Y_regressed[[m]] <- list()
    for (g in groups) {
      if (any(rowSums(!is.na(Y[[m]][[g]]))<min_observations) ) stop(sprintf("Some features do not have enough observations (N=%s) to fit the linear model",min_observations))
      Y_regressed[[m]][[g]] <- t( apply(Y[[m]][[g]], 2, function(y) {
        
        # Fit linear model
        df <- cbind(y,covariates[[m]][[g]])
        lm.out <- lm(y~., data=df)
        residuals <- lm.out[["residuals"]]
        
        # Fill missing values
        all_samples <- rownames(Y[[m]][[g]])
        missing_samples <- all_samples[!all_samples %in% names(residuals)]
        residuals[missing_samples] <- NA
        residuals[all_samples]
      }))
    }
  }
  object@input_data[views] <- Y_regressed
  
  return(object)
}