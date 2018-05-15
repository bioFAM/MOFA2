
##############################################
## Functions to run BioFAM from the R package ##
##############################################

#' @title train an untrained BioFAModel
#' @name run_bio_fam
#' @description Function to train an untrained \code{\link{BioFAModel}} object.
#' @details In this step the R package is calling the \code{BioFAM} Python package, where the the training is performed. \cr
#' The interface with Python is not great for now. Sometimes the training might freeze, but it is running from the background, just be patient! \cr
#' It is important the \code{biofam} executable is on your $PATH so that R can detect it. This will generally be the case, 
#' but if it is not, then you need to specify it manually in the \code{biofam_path} argument. Please, read our FAQ for troubleshooting. \cr
#' @param object an untrained \code{\link{BioFAModel}} object
#' @param dir_options list with I/O options, it should contain at least two entries: \cr
#' 'dataDir' where the input matrices as stored as text files. \cr
#' 'outFile' where the model is going to be stored as an hdf5 file.
#' @param ... Extra options to add to the biofam command
#' @param biofam_path Path the the biofam executable. To be modified if the \code{BioFAM} Python package is not recognised by Rstudio.
#' @return a trained \code{\link{BioFAModel}} object
#' @export
run_bio_fam <- function(object, dir_options, ..., biofam_path="biofam") {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  stopifnot(all(c("data_dir", "output_file") %in% names(dir_options)))
  if (object@Status=="trained") { stop("The model is already trained! If you want to retrain, create a new untrained BioFAModel") }
  
  arglist <- list(
    input_files           = paste0(dir_options$data_dir, "/", viewNames(object), ".txt"),
    output_file           = dir_options$output_file,
    header_cols           = TRUE,
    header_rows           = TRUE,
    delimiter             = object@data_options$delimiter,
    views                 = view_names(object),
    groups                = group_names(object),

    likelihoods           = object@model_options$likelihood,
    factors               = object@model_options$num_factors,
    feature_wise_noise    = object@model_options$feature_wise_noise,
    feature_wise_sparsity = object@model_options$feature_wise_sparsity,
    sample_wise_noise     = object@model_options$sample_wise_noise,
    sample_wise_sparsity  = object@model_options$sample_wise_sparsity,

    iter                  = object@training_options$maxiter,
    drop_r2               = object@training_options$drop_factor_threshold,
    tolerance             = object@training_options$tolerance
  )
  
  # Decide whether to learn factors
  if (!object@training_options$learn_factors) {
    arglist$drop_r2 <- 0.00
  } else {
    if (arglist$drop_r2==0) { 
      print("Warning: learn_factors is set to TRUE but drop_factor_threshold is 0, this is contradictory.")
      print("Please read the documentation in runBioFAM about how to learn the number of factors.")
    }
  }

  # Setting the below arguments to NULL doesn't actually add them to
  # the argument list, but reserves that argument name to prevent
  # extra.arglist from using it.
  if (!is.null(object@model_options$covariates)) {
    arglist$covariates_file <- file.path(dir_options$data_dir, "covariates.txt")
    arglist$scale_covariates <- rep(1, ncol(object@model_options$covariates))
  } else {
    arglist$covariates_file <- NULL
    arglist$scale_covariates <- NULL
  }
  arglist$learn_intercept <- as.logical(object@model_options$learn_intercept)
  if (! object@model_options$sparsity) {
    arglist$learn_theta <- rep(0, object@dimensions$M) 
  } else {
    arglist$learn_theta <- NULL
  }

  arglist$center_features <- as.logical(object@data_options$center_features)
  arglist$scale_views <- as.logical(object@data_options$scale_views)
  if (object@data_options$remove_incomplete_samples == T) { command <- paste(command, "--remove-incomplete-samples", sep=" ") }
  arglist$verbose <- as.logical(object@training_options$verbose)

  extra.arglist <- list(...)
  
  # WARNING: is.na() applied to non-(list or vector) of type 'NULL
  # if (any(is.na(names(extra.arglist)) | names(extra.arglist) == "")) {
  #   stop("All extra options must be named")
  # }

  # Remove leading "--" from extra arg names if present (it will be
  # added back later)
  names(arglist) <- sub("^--", "", names(arglist))
  conflicting.argnames <- intersect(names(extra.arglist), names(arglist))
  if (length(conflicting.argnames) > 0)
    stop(paste0("You cannot pass the following arguments as extra options to runBioFAM: ",
      deparse(conflicting.argnames)))

  # No conflicting argument names,
  arglist <- c(arglist, extra.arglist)

  argv <- character(0)
  for (argname in names(arglist)) {
    argval  <- arglist[[argname]]
    argname <- paste0("--", gsub('_', '-', argname))

    if (is.null(argval)) {
      # Placeholder option; don't add it
    }
    if (is.logical(argval)) {
      # Flag option
      if (length(argval) != 1) {
        stop(paste("Invalid argument value:", deprase(argval)))
      } else if (argval == FALSE || is.na(argval)) {
        # Unset flag: don't add it
      } else if (argval == TRUE) {
        # Set flag: add it
        argv <- c(argv, argname)
      }
    } else {
      # Option with arguments: add the option followed by it args
      argv <- c(argv, argname, argval)
    }
  }
  argv <- unlist(argv)

  if (length(biofam_path) != 1) stop("Invalid biofam_path")

  # If output already exists, remove it
  if (file.exists(dir_options$output_file)) {
    if (arglist$verbose) {
      message("Deleting old output file")
      }
    file.remove(dir_options$output_file)
  }

  if (arglist$verbose) {
    message("Running biofam command: ", paste(collapse=" ", shQuote(c(biofam_path, argv))))
  }
  # Run!
  exitcode <- system2(command=biofam_path, args=shQuote(argv), wait=T)
  if (exitcode != 0) {
    stop(paste("biofam command failed with exit code", exitcode))
  }
  
  # Load trained model
  object <- loadModel(dir_options$output_file, object)
  
  return(object)
}
