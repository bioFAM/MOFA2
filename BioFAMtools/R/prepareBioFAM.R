
#######################################################
## Functions to prepare a bioFAM object for training ##
#######################################################

#' @title prepare a BioFAModel for training
#' @name prepareBioFAM
#' @description Function to prepare a \code{\link{BioFAModel}} object for training. 
#' It requires defining data, model and training options.
#' @param object an untrained \code{\link{BioFAModel}}
#' @param DirOptions list with Input/Output options, which must contain a 'data_dir' element where temporary text files will be stored 
#' and an 'outFile' with hdf5 extension where the final model will be stored.
#' @param DataOptions list of DataOptions (see \code{\link{getDefaultDataOptions}} details). 
#' If NULL, default data options are used.
#' @param ModelOptions list of ModelOptions (see \code{\link{getDefaultModelOptions}} for details). 
#' If NULL, default model options are used.
#' @param TrainOptions list of TrainOptions (see \code{\link{getDefaultTrainOptions}} for details). 
#' If NULL, default training options are used.
#' @return Returns an untrained \code{\link{BioFAModel}} with specified data, model and training options
#' @export
#' @examples
#' # Generate a random data set
#' data <- list("gaussian_view"=matrix(rnorm(100,mean=0,sd=1),10,10), "poisson_view"=matrix(rpois(100, lambda=1),10,10))
#' # Create a BioFAM object
#' BioFAMobject <- createBioFAM(data)
#' # Define I/O options
#' DirOptions <- list("data_dir" = tempdir(), "output_file" = tempfile())
#' # Define Data Options
#' DataOptions <- getDefaultDataOptions()
#' # Define Model Options
#' ModelOptions <- getDefaultModelOptions(BioFAMobject)
#' ModelOptions$likelihood <- ("gaussian", "poisson")
#' ModelOptions <- getDefaultModelOptions(BioFAMobject)
#' # Define Training Options
#' TrainOptions <- getDefaultTrainOptions()
#' TrainOptions$maxiter <- 100
#' # Prepare BioFAM object for training
#' BioFAMobject <- prepareBioFAM(BioFAMobject, DirOptions, DataOptions, ModelOptions, TrainOptions)
#' # train a BioFAModel
#' BioFAMobject <- runBioFAM(BioFAMobject, DirOptions)
prepareBioFAM <- function(object, DirOptions, DataOptions = NULL, ModelOptions = NULL, TrainOptions = NULL) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  
  # Create temporary folder to store data
  dir.create(DirOptions$data_dir, showWarnings = FALSE)
  
  # Get data options
  message("Checking data options...")
  if (is.null(DataOptions)) {
    message("No data options specified, using default...")
    object@DataOptions <- getDefaultDataOptions()
  } else {
    if (!is(TrainOptions,"list") & !all(names(TrainOptions) == names(getDefaultTrainOptions())))
      stop("DataOptions are incorrectly specified, please read the documentation in getDefaultDataOptions")
    object@DataOptions <- DataOptions
  }
  
  # Get training options
  message("Checking training options...")
  if (is.null(TrainOptions)) {
    message("No training options specified, using default...")
    object@TrainOptions <- getDefaultTrainOptions()
  } else {
    if(!is(TrainOptions,"list") & !all(names(TrainOptions) == names(getDefaultTrainOptions())))
      stop("TrainOptions are incorrectly specified, please read the documentation in getDefaultTrainOptions")
    object@TrainOptions <- TrainOptions
  }
  
  # Get model options
  message("Checking model options...")
  if(is.null(ModelOptions)) {
    message("No model options specified, using default...")
    object@ModelOptions <- getDefaultModelOptions(object)
  } else {
    # (To-do) Check that ModelOptions is correct
    if(!is(ModelOptions,"list") & !all(names(ModelOptions) == names(getDefaultModelOptions(object))))
      stop("ModelOptions are incorrectly specified, please read the documentation in getDefaultModelOptions")
    object@ModelOptions <- ModelOptions
  }
  
  # Store views as matrices in .txt files
  message(sprintf("Storing input views in %s...", DirOptions$dataDir))
  for(view in viewNames(object)) {
    write.table(t(object@TrainData[[view]]), file=file.path(DirOptions$dataDir, paste0(view,".txt")),
                sep="\t", row.names=T, col.names=T, quote=F)
  }
  
  # If output already exists, remove it
  if (file.exists(DirOptions$output_file))
    file.remove(DirOptions$output_file)
  
  return(object)
}



#' @title Get default training options
#' @name getDefaultTrainOptions
#' @description Function to obtain the default training options.
#' @details The training options are the following: \cr
#' \itemize{
#'  \item{\strong{maxiter}:}{ numeric value indicating the maximum number of iterations. 
#'  Default is 5000, but we recommend using the 'tolerance' as convergence criteria.}
#'  \item{\strong{tolerance}:}{ numeric value indicating the convergence threshold based on the change in Evidence Lower Bound (deltaELBO). 
#'  For quick exploration we recommend this to be around 1.0, and for a thorough training we recommend a value of 0.01. Default is 0.1}
#'  \item{\strong{learn_factors}:}{ logical indicating whether to learn the number of factors. 
#'  The criteria to shut down factors is based on a minimum fraction of variance explained, defined in the \code{DropFactorThreshold} option}
#'  \item{\strong{drop_factor_threshold}:}{ numeric indicating the threshold on fraction of variance explained to consider a factor inactive and drop it from the model.
#'  For example, a value of 0.01 implies that factors explaining less than 1\% of variance (in each view) will be dropped.}
#'  \item{\strong{verbose}:}{ logical indicating whether to generate a verbose output.}
#'  \item{\strong{seed}:}{ random seed for reproducibility (default is NULL, which samples a random seed).}
#' }
#' @return Returns a list with default training options
#' @export
getDefaultTrainOptions <- function() {
  TrainOptions <- list(
    maxiter = 5000,                # (numeric) Maximum number of iterations
    tolerance = 0.1,               # (numeric) Convergence threshold based on change in the evidence lower bound
    learn_factors = TRUE,          # (logical) learn the number of factors?
    drop_factor_threshold = 0.02,  # (numeric) Threshold on fraction of variance explained to drop a factor
    verbose = FALSE,               # (logical) verbosity?
    seed = NULL                    # (numeric or NULL) random seed
  )
  return(TrainOptions)
}




#' @title Get default data options
#' @name getDefaultDataOptions
#' @description Function to obtain the default data options.
#' @details The data options are the following: \cr
#' \itemize{
#'  \item{\strong{sclae_views}:}{ logical indicating whether to scale views to have the same unit variance. 
#'  As long as the scale differences between the data sets is not too high, this is not required. Default is FALSE.}
#'  \item{\strong{center_features}:}{ logical indicating whether to center the features to zero mean. This only works for gaussian data. Default is TRUE.}
#' }
#' @return Returns a list with the default data options, which have to be passed as an argument to \code{\link{prepareMOFA}}
#' @export
getDefaultDataOptions <- function() {
  DataOptions <- list(
    center_features = TRUE,   # Center features to zero mean (only applies to continuous data)
    scale_views = FALSE,      # Scale views to unit variance (only applies to continuous data)
  )
  return(DataOptions)
}

#' @title Get default model options
#' @name getDefaultModelOptions
#' @param object an untrained \code{\link{MOFAmodel}} object
#' @description Function to obtain the default model options.
#' @details The model options are the following: \cr
#' \itemize{
#'  \item{\strong{likelihood}:}{ character vector with data likelihoods per view: 
#'  'gaussian' for continuous data, 'bernoulli' for binary data and 'poisson' for count data.
#'  By default, they are guessed internally.}
#'  \item{\strong{num_factors}:}{ numeric value indicating the initial number of factors. 
#'  If you want to learn the number of factors automatically we recommend setting this to a large value, around 50. Default is 25.}
#'  \item{\strong{learn_intercept}:}{ logical indicating whether to learn an intercept term to capture differences in feature means.
#'  This prevents you from having to center the data, so this option is always recommended. Default is TRUE.}
#' }
#' @return Returns a list with the default model options, which have to be passed as an argument to \code{\link{prepareMOFA}}
#' @export
getDefaultModelOptions <- function(object) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  if (!.hasSlot(object,"Dimensions") | length(object@Dimensions) == 0) stop("Dimensions of object need to be defined before getting ModelOptions")
  if (!.hasSlot(object,"InputData")) stop("InputData slot needs to be specified before getting ModelOptions")
  if (!.hasSlot(object,"TrainData")) stop("TrainData slot needs to be specified before getting ModelOptions")
  
  # Guess likelihood type
  likelihood <- .inferLikelihoods(object)
  
  # Define default model options
  ModelOptions <- list(
    likelihood = likelihood,    # (character vector) likelihood per view [gaussian/bernoulli/poisson]
    learn_intercept = TRUE,     # (bool) include a constant factor of 1s to learn the mean of features (intercept)? If not, you need to center the data
    num_factors = 25            # (numeric) initial number of latent factors
  )
  
  return(ModelOptions)
}
