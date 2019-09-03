#######################################
## Functions to train a MOFA model ##
#######################################

#' @title train a MOFA model
#' @name run_mofa
#' @description Function to train an untrained \code{\link{MOFA}} object.
#' @details In this step the R package is calling the \code{mofapy2} Python package, where the the training is performed. \cr
#' The interface with Python is done with the \code{\link{reticulate}} package. 
#' If you have several versions of Python installed and Rstudio is not detecting the correct one, you can change it using
#' \code{reticulate::use_python}. \cr
#' This module is in beta testing so please, read our FAQ for troubleshooting and report any problems.
#' @param object an untrained \code{\link{MOFA}} object
#' @param outfile output file for the model (.hdf5 format). If NULL, a temporary file is created.
#' @return a trained \code{\link{MOFA}} object
#' @import reticulate
#' @export
run_mofa <- function(object, outfile = NA) {
  
  # Sanity checks
  if (!is(object, "MOFA")) 
    stop("'object' has to be an instance of MOFA")
  if (object@status=="trained") 
    stop("The model is already trained! If you want to retrain, create a new untrained MOFA")

  # If not outfile is provided, store a file in the /temp folder with the respective timestamp
  if (is.na(outfile)) {
    outfile <- file.path("/tmp", paste0("mofa_", format(Sys.time(), format = "%Y%m%d-%H%M%S"), ".hdf5"))
    warning(paste0("No output filename provided. Using ", outfile, " to store the trained model.\n\n"))
  }
  if (file.exists(outfile))
    message("Warning: Output file already exists, it will be replaced")
  
  # Initiate reticulate
  mofa <- import("mofapy2")
  
  # Call entry point
  mofa_entrypoint <- mofa$run.entry_point$entry_point()
  
  # Set data options
  mofa_entrypoint$set_data_options(
    likelihoods = unname(object@model_options$likelihoods),
    scale_views = object@data_options$scale_views
  )
  
  # Set the data
  mofa_entrypoint$set_data_matrix(
    data = unname(lapply(object@data, function(x) lapply(x, function(y) r_to_py(t(y)) ))),
    views = r_to_py(as.list(object@data_options$views)),
    groups = r_to_py(as.list(object@data_options$groups)),
    samples = r_to_py(unname(lapply(object@data[[1]], colnames))),
    features = r_to_py(unname(lapply(object@data, function(x) rownames(x[[1]]))))
  )
  
  # Set model options 
  mofa_entrypoint$set_model_options(
    factors     = object@model_options$num_factors,
    likelihoods = unname(object@model_options$likelihood),
    spikeslab_z = object@model_options$spikeslab_factors, 
    spikeslab_w = object@model_options$spikeslab_weights, 
    ard_z       = object@model_options$ard_factors,
    ard_w       = object@model_options$ard_weights 
  )
  
  # Set training options  
  mofa_entrypoint$set_train_options(
    iter             = object@training_options$maxiter,
    convergence_mode = object@training_options$convergence_mode,
    # dropR2           = object@training_options$drop_factor_threshold,
    startELBO        = object@training_options$startELBO,
    elbofreq         = object@training_options$freqELBO,
    seed             = object@training_options$seed, 
    gpu_mode         = object@training_options$gpu_mode,
    verbose          = object@training_options$verbose
  )
  
  # Set stochastic options
  if (object@training_options$stochastic) {
    mofa_entrypoint$set_stochastic_options(
      learning_rate    = object@stochastic_options$learning_rate,
      forgetting_rate  = object@stochastic_options$forgetting_rate,
      batch_size       = object@stochastic_options$batch_size
    )
  }
  
  # Build the model
  mofa_entrypoint$build()
  
  # Run the model
  mofa_entrypoint$run()
  
  # Save the model output as an hdf5 file
  mofa_entrypoint$save(outfile)
  
  # Load the trained model
  object <- load_model(outfile)
  
  return(object)
}
