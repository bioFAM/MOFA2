#######################################
## Functions to train a MOFA model ##
#######################################

#' @title Train a MOFA model
#' @name run_mofa
#' @description Function to train an untrained \code{\link{MOFA}} object.
#' @details This function is called once a MOFA object has been prepared (using \code{\link{prepare_mofa}})
#' In this step the R package calls the \code{mofapy2} Python package, where model training is performed. \cr
#' The interface with Python is done with the \code{\link{reticulate}} package. 
#' If you have several versions of Python installed and R is not detecting the correct one, 
#' you can change it using \code{reticulate::use_python} when loading the R session. 
#' Alternatively, you can let us install mofapy2 for you using \code{basilisk} if you set use_basilisk to \code{TRUE}
#' @param object an untrained \code{\link{MOFA}} object
#' @param save_data logical indicating whether to save the training data in the hdf5 file. 
#'  This is useful for some downstream analysis (mainly functions with the prefix \code{plot_data}), but it can take a lot of disk space.
#' @param outfile output file for the model (.hdf5 format). If \code{NULL}, a temporary file is created.
#' @param use_basilisk use \code{basilisk} to automatically install a conda environment with mofapy2 and all dependencies? 
#' If \code{FALSE} (default), you should specify the right python binary when loading R with \code{reticulate::use_python(..., force=TRUE)}
#' or the right conda environment with \code{reticulate::use_condaenv(..., force=TRUE)}.
#' @return a trained \code{\link{MOFA}} object
#' @import reticulate
#' @import basilisk
#' @export
#' @examples
#' # Load data (in data.frame format)
#' file <- system.file("extdata", "test_data.RData", package = "MOFA2")
#' load(file) 
#' 
#' # Create the MOFA object
#' MOFAmodel <- create_mofa(dt)
#' 
#' # Prepare the MOFA object with default options
#' MOFAmodel <- prepare_mofa(MOFAmodel)
#' 
#' # Run the MOFA model
#' \dontrun{ MOFAmodel <- run_mofa(MOFAmodel, use_basilisk = TRUE) }
run_mofa <- function(object, outfile = NULL, save_data = TRUE, use_basilisk = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) 
    stop("'object' has to be an instance of MOFA")
  if (object@status=="trained") 
    stop("The model is already trained! If you want to retrain, create a new untrained MOFA")
  
  # If no outfile is provided, store a file in the /tmp folder with the respective timestamp
  if (is.null(outfile) || is.na(outfile) || (outfile == "")) {
    outfile <- object@training_options$outfile
    if (is.null(outfile) || is.na(outfile) || (outfile == "")) {
      outfile <- file.path("/tmp", paste0("mofa_", format(Sys.time(), format = "%Y%m%d-%H%M%S"), ".hdf5"))
      warning(paste0("No output filename provided. Using ", outfile, " to store the trained model.\n\n"))
    }
  }
  if (file.exists(outfile))
    message(paste0("Warning: Output file ", outfile, " already exists, it will be replaced"))
  
  # Connect to mofapy2 using reticulate (default)
  if (!use_basilisk) {

    message("Connecting to the mofapy2 python package using reticulate (use_basilisk = FALSE)... 
    Please make sure to manually specify the right python binary when loading R with reticulate::use_python(..., force=TRUE) or the right conda environment with reticulate::use_condaenv(..., force=TRUE)
    If you prefer to let us automatically install a conda environment with 'mofapy2' installed using the 'basilisk' package, please use the argument 'use_basilisk = TRUE'\n")
    
    # Sanity checks
    have_mofa2 <- py_module_available("mofapy2")
    if (have_mofa2) {
      .run_mofa_reticulate(object, outfile, save_data)
    } else {
      warning("mofapy2 is not detected in the specified python binary, see reticulate::py_config(). Setting use_basilisk = TRUE...")
      use_basilisk <- TRUE
    }
  }
    
  # Connect to mofapy2 using basilisk (optional)
  if (use_basilisk) {
    
    message("Connecting to the mofapy2 package using basilisk. 
    Set 'use_basilik' to FALSE if you prefer to manually set the python binary using 'reticulate'.")
    
    proc <- basiliskStart(mofa_env)
    on.exit(basiliskStop(proc))
    tmp <- basiliskRun(proc, function(object, outfile, save_data) {
      .run_mofa_reticulate(object, outfile, save_data)
    }, object=object, outfile=outfile, save_data=save_data)
  }
  
  # Load the trained model
  object <- load_model(outfile)
  
  return(object)
}



.run_mofa_reticulate <- function(object, outfile, save_data) {
  
  # sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package \"reticulate\" is required but is not installed.", call. = FALSE)
  }
  
  # Initiate reticulate
  mofa <- import("mofapy2")
  
  # Call entry point
  mofa_entrypoint <- mofa$run.entry_point$entry_point()
  
  # Set data options
  mofa_entrypoint$set_data_options(
    scale_views = object@data_options$scale_views,
    scale_groups = object@data_options$scale_groups
  )

  # Set metadata
  if (.hasSlot(object, "samples_metadata")) {
    mofa_entrypoint$data_opts$samples_metadata <- r_to_py(lapply(object@data_options$groups, function(g) object@samples_metadata[object@samples_metadata$group == g,]))
  }

  if (.hasSlot(object, "features_metadata")) {
    mofa_entrypoint$data_opts$features_metadata <- r_to_py(unname(lapply(object@data_options$views, function(m) object@features_metadata[object@features_metadata$view == m,])))
  }
  
  # Set the data
  mofa_entrypoint$set_data_matrix(
    data = r_to_py( unname(lapply(object@data, function(x) unname( lapply(x, function(y) r_to_py(t(y)) ))) ) ),
    likelihoods = unname(object@model_options$likelihoods),
    views_names = r_to_py(as.list(object@data_options$views)),
    groups_names = r_to_py(as.list(object@data_options$groups)),
    samples_names = r_to_py(unname(lapply(object@data[[1]], colnames))),
    features_names = r_to_py(unname(lapply(object@data, function(x) rownames(x[[1]]))))
  )
  
  # Set model options 
  mofa_entrypoint$set_model_options(
    factors     = object@model_options$num_factors,
    spikeslab_factors = object@model_options$spikeslab_factors, 
    spikeslab_weights = object@model_options$spikeslab_weights, 
    ard_factors       = object@model_options$ard_factors,
    ard_weights       = object@model_options$ard_weights 
  )
  
  # Set training options  
  mofa_entrypoint$set_train_options(
    iter             = object@training_options$maxiter,
    convergence_mode = object@training_options$convergence_mode,
    dropR2           = object@training_options$drop_factor_threshold,
    startELBO        = object@training_options$startELBO,
    freqELBO         = object@training_options$freqELBO,
    seed             = object@training_options$seed, 
    gpu_mode         = object@training_options$gpu_mode,
    verbose          = object@training_options$verbose,
    outfile          = object@training_options$outfile,
    save_interrupted = object@training_options$save_interrupted
  )
  
  # Set stochastic options
  if (object@training_options$stochastic) {
    mofa_entrypoint$set_stochastic_options(
      learning_rate    = object@stochastic_options$learning_rate,
      forgetting_rate  = object@stochastic_options$forgetting_rate,
      batch_size       = object@stochastic_options$batch_size,
      start_stochastic = object@stochastic_options$start_stochastic
    )
  }
  
  # Build the model
  mofa_entrypoint$build()
  
  # Run the model
  mofa_entrypoint$run()

  # Save the model output as an hdf5 file
  mofa_entrypoint$save(outfile, save_data = save_data)
}
