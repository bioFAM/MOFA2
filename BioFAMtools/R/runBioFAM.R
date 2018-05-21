
###########################
## Functions to run MOFA ##
###########################

#' @title train a BioFAM model
#' @name run_biofam
#' @description Function to train an untrained \code{\link{bioFAMmodel}} object.
#' @details In this step the R package is calling the \code{biofam} Python package, where the the training is performed. \cr
#' The interface with Python is done with the \code{\link{reticulate}} package. 
#' If you have several versions of Python installed and Rstudio is not detecting the correct one, you can change it using
#' \code{reticulate::use_python}. \cr
#' This module is in beta testing so please, read our FAQ for troubleshooting and report any problems.
#' @param object an untrained \code{\link{MOFAmodel}} object
#' @param DirOptions list with I/O options, it should contain at least two entries: \cr
#' \itemize{
#'  \item{\strong{data_dir}:}{ input directory where the matrices as stored as text files with row and column names and tab delimiter. }
#'  \item{\strong{outfile}:}{ output file where the model is going to be stored as an hdf5 file.}
#' }
#' @return a trained \code{\link{MOFAmodel}} object
#' @import reticulate
#' @export
run_biofam <- function(object, DirOptions) {
  
  # Sanity checks
  if (!is(object, "bioFAMmodel")) 
    stop("'object' has to be an instance of MOFAmodel")
  
  stopifnot(all(c("data_dir","outfile") %in% names(DirOptions)))
  
  if (object@Status=="trained") 
    stop("The model is already trained! If you want to retrain, create a new untrained MOFAmodel")
  
  if (object@TrainOptions$learn_factors == FALSE) {
    object@TrainOptions$drop_factor_threshold <- 0.00
  } else {
    if (object@TrainOptions$drop_factor_threshold==0) { 
      print("Warning: learn_factors is set to TRUE but drop_factor_threshold is 0, this is contradictory.")
      print("Please read the documentation in prepareMOFA about learning the number of factors.")
    }
  }
  
  if (file.exists(DirOptions$outfile))
    message("Warning: Output file already exists, it will be replaced")
  
  # Initiate reticulate
  mofa <- import("mofa")
  mofa_entrypoint <- mofa$core.init_asd2$entry_point()
  
  # Pass data options
  mofa_entrypoint$set_data_options(
    inFiles     = paste0(DirOptions$data_dir, "/", viewNames(object), ".txt"), 
    outFile     = DirOptions$outfile, 
    views       = viewNames(object), 
    delimiter   = "\t", 
    header_cols = TRUE, 
    header_rows = TRUE
  )
  
  # Pass training options  
  mofa_entrypoint$set_train_options(
    iter       = object@TrainOptions$maxiter,
    tolerance  = object@TrainOptions$tolerance,
    dropR2     = object@TrainOptions$drop_factor_threshold,
    seed       = object@TrainOptions$seed, 
    verbose    = object@TrainOptions$verbose
  )
  
  # Pass model options 
  mofa_entrypoint$set_model_options(
    factors        = object@ModelOptions$num_factors,
    likelihoods    = unname(object@ModelOptions$likelihood),
    learn_intercept = object@ModelOptions$learn_intercept
  )
  
  # Pass data processing options
  mofa_entrypoint$set_dataprocessing_options(
    center_features         = object@DataOptions$center_features,
    scale_views             = object@DataOptions$scale_views
    
  )
  
  # Load data
  mofa_entrypoint$load_data()
  
  # Define Priors
  mofa_entrypoint$define_priors()
  
  # Initialise variational distributions
  mofa_entrypoint$initialise_variational()
  
  # Parse the intercept factor
  mofa_entrypoint$parse_intercept()
  
  # Train the model
  mofa_entrypoint$train_model()
  
  # Load the trained model
  object <- loadModel(DirOptions$outfile, object)
  
  return(object)
}
