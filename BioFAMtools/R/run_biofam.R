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
#' @param outfile output file (use .hdf5 format)
#' @return a trained \code{\link{bioFAMmodel}} object
#' @import reticulate
#' @export
run_biofam <- function(object, outfile = NA) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) 
    stop("'object' has to be an instance of BioFAModel")
  if (object@status=="trained") 
    stop("The model is already trained! If you want to retrain, create a new untrained BioFAModel")

  # If not outfile is provided, store a file in the /temp folder with the respective timestamp
  if (is.na(outfile)) {
    outfile <- file.path("/tmp", paste0("biofam_", format(Sys.time(), format = "%Y%m%d-%H%M%S"), ".hdf5"))
    warning(paste0("No output filename provided. Using ", outfile, " to store the trained model.\n\n"))
  }
  if (file.exists(outfile))
    message("Warning: Output file already exists, it will be replaced")
  
  # Initiate reticulate
  biofam <- import("biofam")
  
  # Call entry point
  biofam_entrypoint <- biofam$run.entry_point$entry_point()
  
  # Set data options
  biofam_entrypoint$set_data_options(
    likelihoods = unname(object@model_options$likelihood),
    center_features_per_group = object@data_options$center_features_per_group,
    scale_views = object@data_options$scale_views
  )
  
  # Set the data
  # np <- import("numpy", convert = FALSE)
  biofam_entrypoint$set_data_matrix(
    data = unname(lapply(object@input_data, function(x) r_to_py(t(x)))),
    # data = lapply(unname(object@input_data), function(x) r_to_py(unname(lapply(x, function(y) np_array(t(y), dtype = np$float64) )))),
    samples_names_dict = r_to_py(lapply(object@input_data[[1]], rownames)),
    features_names_dict = r_to_py(lapply(object@input_data, function(m) colnames(m[[1]])))
  )
  
  # Set model options 
  biofam_entrypoint$set_model_options(
    factors     = object@model_options$num_factors,
    likelihoods = unname(object@model_options$likelihood),
    sl_z        = object@model_options$sl_z, 
    sl_w        = object@model_options$sl_w, 
    ard_w       = object@model_options$ard_w, 
    ard_z       = object@model_options$ard_z
  )
  
  # Set training options  
  biofam_entrypoint$set_train_options(
    iter             = object@training_options$maxiter,
    convergence_mode = object@training_options$convergence_mode,
    # dropR2           = object@training_options$drop_factor_threshold,
    startELBO        = object@training_options$startELBO,
    elbofreq         = object@training_options$freqELBO,
    seed             = object@training_options$seed, 
    gpu_mode          = object@training_options$gpu_mode,
    verbose          = object@training_options$verbose
  )
  
  # Set stochastic options
  if (object@training_options$stochastic) {
    biofam_entrypoint$set_stochastic_options(
      learning_rate     = object@stochastic_options$learning_rate,
      forgetting_rate     = object@stochastic_options$forgetting_rate,
      batch_size     = object@stochastic_options$batch_size
    )
  }
  
  # Build the model
  biofam_entrypoint$build()
  
  # Run the model
  biofam_entrypoint$run()
  
  # Save the model as an hdf5 file
  biofam_entrypoint$save(outfile)
  
  # Load the trained model
  object <- load_model(outfile)
  
  return(object)
}
