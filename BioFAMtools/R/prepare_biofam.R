
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
prepare_biofam <- function(object, data_options = NULL, model_options = NULL, training_options = NULL) {
  
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
  
  # Store views as matrices in .txt files
  # message(sprintf("Storing input views in %s...", dir_options$data_dir))
  # for(view in viewNames(object)) {
  #   write.table(t(object@training_data[[view]]), file=file.path(dir_options$data_dir, paste0(view,".txt")),
  #               sep="\t", row.names=T, col.names=T, quote=F)
  # }
  
  # If output already exists, remove it
  # if (file.exists(dir_options$outfile))
  #   file.remove(dir_options$outfile)
  
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
    drop_factor_threshold = NA,    # (numeric) Threshold on fraction of variance explained to drop a factor
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
    center_features_per_group = TRUE,  # (logical) Center features to zero mean, for each group separately
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
  # TO-DO: RATHER THAN CHECKING IF SLOT EXISTS, WE NEED TO CHECK THAT THE IS NOT AN EMPTY LIST()
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  if (!.hasSlot(object,"dimensions") | length(object@dimensions) == 0) stop("dimensions of object need to be defined before getting the model options")
  if (!.hasSlot(object,"input_data")) stop("input_data slot needs to be specified before getting the model options")
  
  # Guess likelihoods from the data
  # likelihood <- .infer_likelihoods(object) THIS DOES NOT WORK 
  likelihood <- rep(x="gaussian", times=object@dimensions$M)
  
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
