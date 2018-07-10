
###########################
## Functions to run BioFAM ##
###########################

#' @title train a BioFAModel
#' @name run_biofam
#' @description Function to train an untrained \code{\link{BioFAModel}} object.
#' @details In this step the R package is calling the \code{biofam} Python package,
#'  where the the training is performed. \cr
#' The interface with Python is done with the \code{\link{reticulate}} package. 
#' If you have several versions of Python installed and Rstudio is not detecting
#'  the correct one, you can change it using' \code{reticulate::use_python}.
#' @param object an untrained \code{\link{BioFAModel}} object
#' @param outfile output .hdf5 file
#' @return a trained \code{\link{BioFAModel}} object
#' @import reticulate
#' @export

run_biofam <- function(object, outfile=NULL) {
  
  # Sanity checks on the model
  if (!is(object, "BioFAModel")) 
    stop("'object' has to be an instance of BioFAModel")
  
  # Sanity checks on the output file
  if (is.null(outfile)) {
    print("No output file provided, using a temporary file...")
    outfile <- tempfile()
  } else {
    if (!dir.exists(dirname(outfile))) {
      print("Output directory not found, creating it...")
      dir.create(dirname(outfile), recursive = T, showWarnings = T)
    }
  }
    
  if (object@status=="trained") 
    stop("The model is already trained!
         If you want to retrain, create a new untrained BioFAModel")
  
  # Initiate reticulate
  biofam <- import("biofam")
  entrypoint <- biofam$run.entry_point$entry_point()
  
  # Pass data
  entrypoint$set_data(data = unname(lapply(object@training_data, function(view) 
    lapply(view, function(x)
      r_to_py(t(x)))
    )
  ))

  
  # Pass data processing options
  entrypoint$set_data_options(
    lik                       = unname(object@model_options$likelihood),
    center_features           = object@data_options$center_features,
    center_features_per_group = object@data_options$center_features_per_group,
    scale_features            = object@data_options$scale_features,
    scale_views               = object@data_options$scale_views
  )
  
  
  # Pass model options 
  entrypoint$set_model_options(
    factors         = object@model_options$num_factors,
    likelihoods     = unname(object@model_options$likelihood),
    learn_intercept = object@model_options$learn_intercept,
    noise_on        = object@model_options$noise_on,
    ard_z           = object@model_options$ard_z,
    ard_w           = object@model_options$ard_w,
    sl_w            = object@model_options$sl_w,
    sl_z            = object@model_options$sl_z
  )
  
  # Pass training options  
  entrypoint$set_train_options(
    iter       = object@training_options$maxiter,
    tolerance  = object@training_options$tolerance,
    dropR2     = object@training_options$drop_factor_threshold,
    seed       = object@training_options$seed, 
    verbose    = object@training_options$verbose
  )

  # Build the model
  entrypoint$build()
  
  # Train the model
  entrypoint$run()
  
  # Save the model
  samples_names <- colnames(object@training_data[[1]])
  features_names <- unname(lapply(object@training_data,rownames))
  entrypoint$save_model(outfile,
                        samples_names=samples_names,
                        features_names=features_names)
  
  # Load the model back into R
  object <- load_model(outfile, object)
  
  return(object)
}