#######################################
## Functions to train a BioFAM model ##
#######################################

#' @title train a BioFAM model
#' @name run_biofam
#' @description Function to train an untrained \code{\link{bioFAMmodel}} object.
#' @details In this step the R package is calling the \code{biofam} Python package, where the the training is performed. \cr
#' The interface with Python is done with the \code{\link{reticulate}} package. 
#' If you have several versions of Python installed and Rstudio is not detecting the correct one, you can change it using
#' \code{reticulate::use_python}. \cr
#' This module is in beta testing so please, read our FAQ for troubleshooting and report any problems.
#' @param object an untrained \code{\link{bioFAMmodel}} object
#' @param dir_options list with I/O options, it should contain at least two entries: \cr
#' \itemize{
#'  \item{\strong{data_dir}:}{ input directory where the matrices as stored as text files with row and column names and tab delimiter. }
#'  \item{\strong{outfile}:}{ output file where the model is going to be stored as an hdf5 file.}
#' }
#' @return a trained \code{\link{bioFAMmodel}} object
#' @import reticulate
#' @export
run_biofam <- function(object, dir_options, samples_groups = NULL) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) 
    stop("'object' has to be an instance of BioFAModel")
  
  # stopifnot(all(c("data_dir", "outfile") %in% names(dir_options)))
  stopifnot(all(c("outfile") %in% names(dir_options)))
  
  if (object@status=="trained") 
    stop("The model is already trained! If you want to retrain, create a new untrained BioFAModel")
  
  if (object@training_options$learn_factors == FALSE) {
    object@training_options$drop_factor_threshold <- 0.00
  } else {
    if (object@training_options$drop_factor_threshold==0) { 
      print("Warning: learn_factors is set to TRUE but drop_factor_threshold is 0, this is contradictory.")
      print("Please read the documentation in prepare_biofam about learning the number of factors.")
    }
  }
  
  if (file.exists(dir_options$outfile))
    message("Warning: Output file already exists, it will be replaced")
  
  # Initiate reticulate
  biofam <- import("biofam")
  
  # Call entry point
  biofam_entrypoint <- biofam$run.entry_point$entry_point()
  
  # Pass data
  
  # Set data options
  biofam_entrypoint$set_data_options(
    lik = unname(object@model_options$likelihood),
    center_features = object@data_options$center_features,
    center_features_per_group = object@data_options$center_features_per_group,
    scale_views = object@data_options$scale_views
  )
  
  # Set the data
  # biofam_entrypoint$set_data_df(data=unname(lapply(object@input_data, function(x) r_to_py(t(x)))))
  if (.hasSlot(object, "training_data")) {
    biofam_entrypoint$set_data_matrix(r_to_py(object@training_data), 
                                      # r_to_py(names(object@training_data)), r_to_py(names(object@training_data[[1]])),
                                      r_to_py(lapply(object@training_data[[1]], rownames)), r_to_py(lapply(object@training_data, function(m) colnames(m[[1]]))))
  } else {
    biofam_entrypoint$set_data_df(r_to_py(object@input_data))
  }
  
  
  # Set training options  
  biofam_entrypoint$set_train_options(
    iter       = object@training_options$maxiter,
    tolerance  = object@training_options$tolerance,
    dropR2     = object@training_options$drop_factor_threshold,
    seed       = object@training_options$seed, 
    verbose    = object@training_options$verbose
  )
  
  # Set model options 
  biofam_entrypoint$set_model_options(
    factors         = object@model_options$num_factors,
    likelihoods     = unname(object@model_options$likelihood),
    learn_intercept = object@model_options$learn_intercept,
    sl_z=object@model_options$sl_z, 
    sl_w=object@model_options$sl_w, 
    ard_w=object@model_options$ard_w, 
    ard_z=object@model_options$ard_z
  )
  
  # Build the model
  biofam_entrypoint$build()
  
  # Run the model
  biofam_entrypoint$run()
  
  # Save the model as an hdf5 file
  biofam_entrypoint$save(
    outfile = dir_options$outfile
  )
  
  # Load the trained model
  object <- load_model(dir_options$outfile, object)
  
  return(object)
}
