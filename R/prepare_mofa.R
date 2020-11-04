
#######################################################
## Functions to prepare a MOFA object for training ##
#######################################################

#' @title Prepare a MOFA for training
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
#' @details This function is called after creating a \code{\link{MOFA}} object (using  \code{\link{create_mofa}}) 
#' and before starting the training (using \code{\link{run_mofa}}). Here, we can specify different options for
#' the data (data_options), the model (model_options) and the trainig (training_options, stochastic_options). Take a look at the
#' individual default options for an overview using the get_default_XXX_options functions above.
#' @export
#' @examples
#' # Using an existing simulated data with two groups and two views
#' file <- system.file("extdata", "test_data.RData", package = "MOFA2")
#' 
#' # Load data dt (in data.frame format)
#' load(file) 
#' 
#' # Create the MOFA object
#' MOFAmodel <- create_mofa(dt)
#' 
#' # Prepare MOFA object using default options
#' MOFAmodel <- prepare_mofa(MOFAmodel)
#' 
#' # Prepare MOFA object changing some of the default model options values
#' model_opts <- get_default_model_options(MOFAmodel)
#' model_opts$num_factors <- 10
#' MOFAmodel <- prepare_mofa(MOFAmodel, model_options = model_opts)
prepare_mofa <- function(object, data_options = NULL, model_options = NULL, training_options = NULL, stochastic_options = NULL) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (any(object@dimensions$N<15)) warning("Some group(s) have less than 15 samples, MOFA will have little power to learn meaningful factors for these group(s)...")
  if (any(object@dimensions$D<15)) warning("Some view(s) have less than 15 features, MOFA will have little power to to learn meaningful factors for these view(s)....")
  if (any(object@dimensions$D>1e4)) warning("Some view(s) have a lot of features, it is recommended to performa more stringent feature selection before creating the MOFA object....")
  
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
    
    if (object@training_options$maxiter<=100)
      warning("Maximum number of iterations is very small\n")
    if (object@training_options$startELBO<1) object@training_options$startELBO <- 1
    if (object@training_options$freqELBO<1) object@training_options$freqELBO <- 1
    if (!object@training_options$convergence_mode %in% c("fast","medium","slow")) 
      stop("Convergence mode has to be either 'fast', 'medium', or 'slow'")
  }
  
  # Get stochastic options
  if (is.null(stochastic_options)) {
    object@stochastic_options <- list()
  } else {
    if (isFALSE(object@training_options[["stochastic"]]))
      stop("stochastic_options have been provided but training_opts$stochastic is FALSE. If you want to use stochastic inference you have to set training_opts$stochastic = TRUE")
    # object@training_options$stochastic <- TRUE
  }
  
  if (object@training_options$stochastic) {
    message("Stochastic inference activated. Note that this is only recommended if you have a very large sample size (>1e4) and access to a GPU")
      
    if (is.null(stochastic_options)) {
      message("No stochastic options specified, using default...")
      object@stochastic_options <- get_default_stochastic_options(object)
    } else {
      object@training_options$stochastic <- TRUE
      message("Checking stochastic inference options...")
      if (!is(stochastic_options,"list") || !setequal(names(stochastic_options), names(get_default_stochastic_options(object)) ))
        stop("stochastic_options are incorrectly specified, please read the documentation in get_default_stochastic_options")
      
      if (!stochastic_options$batch_size %in% c(0.05,0.10,0.15,0.20,0.25,0.50))
        stop("Batch size has to be one of the following numeric values: 0.05, 0.10, 0.15, 0.20, 0.25, 0.50")
      if (stochastic_options$batch_size==1)
        warning("A batch size equal to 1 is equivalent to non-stochastic inference.")
      if (stochastic_options$learning_rate<=0 || stochastic_options$learning_rate>1)
        stop("The learning rate has to be a value between 0 and 1")
      if (stochastic_options$forgetting_rate<=0 || stochastic_options$forgetting_rate>1)
        stop("The forgetting rate has to be a value between 0 and 1")
  
      if (sum(object@dimensions$N)<1e4) warning("Stochastic inference is only recommended when you have a lot of samples (at least N>10,000))\n")
      
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
  # if (!object@model_options$ard_weights) warning("model_options$ard_weights should always be set to TRUE")
  if (min(object@dimensions$N) < 4 * object@model_options$num_factors) {
    warning(sprintf("The number of samples in at least one group is very small for learning %s factors.  
    Try to reduce the number of factors to obtain meaningful results. It should not exceed ~%s.",
    object@model_options$num_factors, floor(min(object@dimensions$N/4))))
  }

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
#' @details This function provides a default set of training options that can be modified and passed to the \code{\link{MOFA}} object
#' in the \code{\link{prepare_mofa}} step (see example), i.e. after creating a \code{\link{MOFA}} object
#'  (using \code{\link{create_mofa}}) and before starting the training (using \code{\link{run_mofa}})
#' The training options are the following: \cr
#' \itemize{
#'  \item{\strong{maxiter}:}{ numeric value indicating the maximum number of iterations. 
#'  Default is 1000. Convergence is assessed using the ELBO statistic.}
#'  \item{\strong{drop_factor_threshold}:}{ numeric indicating the threshold on fraction of variance explained to consider a factor inactive and drop it from the model.
#'  For example, a value of 0.01 implies that factors explaining less than 1\% of variance (in each view) will be dropped. Default is -1 (no dropping of factors)}
#'  \item{\strong{convergence_mode}:}{ character indicating the convergence criteria, either "slow", "medium" or "fast", corresponding to 5e-7\%, 5e-6\% or 5e-5\% deltaELBO change w.r.t. to the ELBO at the first iteration. }
#'  \item{\strong{verbose}:}{ logical indicating whether to generate a verbose output.}
#'  \item{\strong{startELBO}:}{ integer indicating the first iteration to compute the ELBO (default is 1). }
#'  \item{\strong{freqELBO}:}{ integer indicating the first iteration to compute the ELBO (default is 1). }
#'  \item{\strong{stochastic}:}{ logical indicating whether to use stochastic variational inference (only required for very big data sets, default is \code{FALSE}).}
#'  \item{\strong{gpu_mode}:}{ logical indicating whether to use GPUs (see details).}
#'  \item{\strong{seed}:}{ numeric indicating the seed for reproducibility (default is 42).}
#' }
#' @return Returns a list with default training options
#' @importFrom utils modifyList
#' @export
#' @examples
#' # Using an existing simulated data with two groups and two views
#' file <- system.file("extdata", "test_data.RData", package = "MOFA2")
#' 
#' # Load data dt (in data.frame format)
#' load(file) 
#' 
#' # Create the MOFA object
#' MOFAmodel <- create_mofa(dt)
#' 
#' # Load default training options
#' train_opts <- get_default_training_options(MOFAmodel)
#' 
#' # Edit some of the training options
#' train_opts$convergence_mode <- "medium"
#' train_opts$startELBO <- 100
#' train_opts$seed <- 42
#' 
#' # Prepare the MOFA object
#' MOFAmodel <- prepare_mofa(MOFAmodel, training_options = train_opts)
get_default_training_options <- function(object) {
  
  # Get default train options
  training_options <- list(
    maxiter = 1000,                # (numeric) Maximum number of iterations
    convergence_mode = 'fast',     # (string) Convergence mode based on change in the ELBO ("slow","medium","fast")
    drop_factor_threshold = -1,    # (numeric) Threshold on fraction of variance explained to drop a factor
    verbose = FALSE,               # (logical) Verbosity
    startELBO = 1,                 # (numeric) First iteration to compute the ELBO
    freqELBO = 1,                  # (numeric) Frequency of ELBO calculation
    stochastic = FALSE,            # (logical) Do stochastic variational inference?
    gpu_mode = FALSE,              # (logical) Use GPU?
    seed = 42,                     # (numeric) Random seed
    outfile = NULL,                # (string)  Output file name
    save_interrupted = FALSE       # (logical) Save partially trained model when training is interrupted?
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
#' @details This function provides a default set of data options that can be modified and passed to the \code{\link{MOFA}} object
#' in the \code{\link{prepare_mofa}} step (see example), i.e. after creating a \code{\link{MOFA}} object
#'  (using \code{\link{create_mofa}}) and before starting the training (using \code{\link{run_mofa}})
#' The data options are the following: \cr
#' \itemize{
#'  \item{\strong{scale_views}:}{ logical indicating whether to scale views to have the same unit variance. 
#'  As long as the scale differences between the views is not too high, this is not required. Default is FALSE.}
#'  \item{\strong{scale_groups}:}{ logical indicating whether to scale groups to have the same unit variance. 
#'  As long as the scale differences between the groups is not too high, this is not required. Default is FALSE.}
#' }
#' @return Returns a list with the default data options.
#' @importFrom utils modifyList
#' @export
#' @examples
#' # Using an existing simulated data with two groups and two views
#' file <- system.file("extdata", "test_data.RData", package = "MOFA2")
#' 
#' # Load data dt (in data.frame format)
#' load(file) 
#' 
#' # Create the MOFA object
#' MOFAmodel <- create_mofa(dt)
#' 
#' # Load default data options
#' data_opts <- get_default_data_options(MOFAmodel)
#' 
#' # Edit some of the data options
#' data_opts$scale_views <- TRUE
#' 
#' # Prepare the MOFA object
#' MOFAmodel <- prepare_mofa(MOFAmodel, data_options = data_opts)
get_default_data_options <- function(object) {
  
  # Define default data options
  data_options <- list(
    scale_views = FALSE,                 # (logical) Scale views to unit variance
    scale_groups = FALSE                 # (logical) Scale groups to unit variance
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
#' @details This function provides a default set of model options that can be modified and passed to the \code{\link{MOFA}} object
#' in the \code{\link{prepare_mofa}} step (see example), i.e. after creating a \code{\link{MOFA}} object
#'  (using \code{\link{create_mofa}}) and before starting the training (using \code{\link{run_mofa}})
#' The model options are the following: \cr
#' \itemize{
#'  \item{\strong{likelihoods}:}{ character vector with data likelihoods per view: 
#'  'gaussian' for continuous data (Default for all views), 'bernoulli' for binary data and 'poisson' for count data.}
#'  \item{\strong{num_factors}:}{ numeric value indicating the (initial) number of factors. Default is 15.}
#'  \item{\strong{spikeslab_factors}:}{ logical indicating whether to use spike and slab sparsity on the factors (Default is FALSE)}
#'  \item{\strong{spikeslab_weights}:}{ logical indicating whether to use spike and slab sparsity on the weights (Default is TRUE)}
#'  \item{\strong{ard_factors}:}{ logical indicating whether to use ARD sparsity on the factors (Default is TRUE only if using multiple groups)}
#'  \item{\strong{ard_weights}:}{ logical indicating whether to use ARD sparsity on the weights (Default is TRUE)}
#'  }
#' @return Returns a list with the default model options.
#' @importFrom utils modifyList
#' @export
#' @examples
#' # Using an existing simulated data with two groups and two views
#' file <- system.file("extdata", "test_data.RData", package = "MOFA2")
#' 
#' # Load data dt (in data.frame format)
#' load(file) 
#' 
#' # Create the MOFA object
#' MOFAmodel <- create_mofa(dt)
#' 
#' # Load default model options
#' model_opts <- get_default_model_options(MOFAmodel)
#' 
#' # Edit some of the model options
#' model_opts$num_factors <- 10
#' model_opts$spikeslab_weights <- FALSE
#' 
#' # Prepare the MOFA object
#' MOFAmodel <- prepare_mofa(MOFAmodel, model_options = model_opts)
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
  # likelihoods <- .infer_likelihoods(object)
  likelihoods <- rep(x="gaussian", times=object@dimensions$M)
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




#' @title Get default stochastic options
#' @name get_default_stochastic_options
#' @description Function to obtain the default options for stochastic variational inference.
#' @param object an untrained \code{\link{MOFA}}
#' @details This function provides a default set of stochastic inference options that can be modified and passed to the \code{\link{MOFA}} object
#' in the \code{\link{prepare_mofa}} step), i.e. after creating a \code{\link{MOFA}} object
#'  (using \code{\link{create_mofa}}) and before starting the training (using \code{\link{run_mofa}})
#' These options are only relevant when activating stochastic inference in training_options (see example).
#' The stochastic inference options are the following: \cr
#' \itemize{
#'  \item{\strong{batch_size}:}{ numeric value indicating the batch size (as a fraction)}. 
#'  Default is 0.5 (half of the data set).
#'  \item{\strong{learning_rate}:}{ numeric value indicating the learning rate. }
#'  Default is 1.0
#'  \item{\strong{forgetting_rate}:}{ numeric indicating the forgetting rate.}
#'  Default is 0.5
#'  \item{\strong{start_stochastic}:}{ integer indicating the first iteration to start stochastic inference}
#'  Default is 1
#'  }
#' @return Returns a list with default options
#' @importFrom utils modifyList
#' @export
#' @examples
#' # Using an existing simulated data with two groups and two views
#' file <- system.file("extdata", "test_data.RData", package = "MOFA2")
#' 
#' # Load data dt (in data.frame format)
#' load(file) 
#' 
#' # Create the MOFA object
#' MOFAmodel <- create_mofa(dt)
#' 
#' # activate stochastic inference in training options
#' train_opts <- get_default_training_options(MOFAmodel)
#' train_opts$stochastic <- TRUE
#' 
#' # Load default stochastic options
#' stochastic_opts <- get_default_stochastic_options(MOFAmodel)
#' 
#' # Edit some of the stochastic options
#' stochastic_opts$learning_rate <- 0.75
#' stochastic_opts$batch_size <- 0.25
#' 
#' # Prepare the MOFA object
#' MOFAmodel <- prepare_mofa(MOFAmodel, 
#'   training_options = train_opts,
#'   stochastic_options = stochastic_opts
#' )
#' 
get_default_stochastic_options <- function(object) {
  
  # Get default stochastic options
  stochastic_options <- list(
    batch_size = 0.5,        # Batch size (as a fraction)
    learning_rate = 1.0,     # Starting learning rate
    forgetting_rate = 0.5,   # Forgetting rate
    start_stochastic = 1     # First iteration to start stochastic inference
  )
  
  # if stochastic_options already exist, replace the default values but keep the additional ones
  if (length(object@stochastic_options)>0)
    stochastic_options <- modifyList(stochastic_options, object@stochastic_options)
  
  return(stochastic_options)
}

