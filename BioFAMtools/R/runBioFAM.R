
##############################################
## Functions to run BioFAM from the R package ##
##############################################

#' @title train an untrained BioFAModel
#' @name runBioFAM
#' @description Function to train an untrained \code{\link{BioFAModel}} object.
#' @details In this step the R package is calling the \code{BioFAM} Python package, where the the training is performed. \cr
#' The interface with Python is not great for now. Sometimes the training might freeze, but it is running from the background, just be patient! \cr
#' It is important the \code{biofam} executable is on your $PATH so that R can detect it. This will generally be the case, 
#' but if it is not, then you need to specify it manually in the \code{biofamPath} argument. Please, read our FAQ for troubleshooting. \cr
#' @param object an untrained \code{\link{BioFAModel}} object
#' @param DirOptions list with I/O options, it should contain at least two entries: \cr
#' 'dataDir' where the input matrices as stored as text files. \cr
#' 'outFile' where the model is going to be stored as an hdf5 file.
#' @param ... Extra options to add to the biofam command
#' @param biofamPath Path the the biofam executable. To be modified if the \code{BioFAM} Python package is not recognised by Rstudio.
#' @return a trained \code{\link{BioFAModel}} object
#' @export
runBioFAM <- function(object, DirOptions, ..., biofamPath="biofam") {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  stopifnot(all(c("dataDir","outFile") %in% names(DirOptions)))
  if (object@Status=="trained") { stop("The model is already trained! If you want to retrain, create a new untrained BioFAModel") }
  
  arglist <- list(
    inFiles = paste0(DirOptions$dataDir, "/", viewNames(object), ".txt"),
    header_cols = TRUE,
    header_rows = TRUE,
    delimiter = object@DataOptions$delimiter,
    outFile = DirOptions$outFile,
    views = viewNames(object),
    likelihoods = object@ModelOptions$likelihood,
    factors = object@ModelOptions$numFactors,
    iter = object@TrainOptions$maxiter,
    dropR2 =  object@TrainOptions$DropFactorThreshold,
    tolerance = object@TrainOptions$tolerance
  )
  
  # Decide whether to learn factors
  if (object@TrainOptions$learnFactors == F) {
    arglist$dropR2 <- 0.00
  } else {
    if (arglist$dropR2==0) { 
      print("Warning: learnFactors is set to TRUE but dropFactorThreshold is 0, this is contradictory.")
      print("Please read the documentation in prepareBioFAM about how to learn the number of factors.")
    }
  }

  # Setting the below arguments to NULL doesn't actually add them to
  # the argument list, but reserves that argument name to prevent
  # extra.arglist from using it.
  if (!is.null(object@ModelOptions$covariates)) {
    arglist$covariatesFile <- file.path(DirOptions$dataDir, "covariates.txt")
    arglist$scale_covariates <- rep(1,ncol(object@ModelOptions$covariates))
  } else {
    arglist$covariatesFile <- NULL
    arglist$scale_covariates <- NULL
  }
  arglist$learnIntercept <- as.logical(object@ModelOptions$learnIntercept)
  if (! object@ModelOptions$sparsity) {
    arglist$learnTheta <- rep(0, object@Dimensions$M) 
  } else {
    arglist$learnTheta <- NULL
  }

  arglist$center_features <- as.logical(object@DataOptions$centerFeatures)
  arglist$scale_views <- as.logical(object@DataOptions$scaleViews)
  if (object@DataOptions$removeIncompleteSamples == T) { command <- paste(command, "--RemoveIncompleteSamples", sep=" ") }
  arglist$verbose <- as.logical(object@TrainOptions$verbose)

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
    argval <- arglist[[argname]]
    argname <- paste0("--", argname)

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

  if (length(biofamPath) != 1) stop("Invalid biofamPath")

  # If output already exists, remove it
  if (file.exists(DirOptions$outFile)) {
    if (arglist$verbose) {
      message("Deleting old output file")
      }
    file.remove(DirOptions$outFile)
  }

  if (arglist$verbose) {
    message("Running BioFAM command: ", paste(collapse=" ", shQuote(c(biofamPath, argv))))
  }
  # Run!
  exitcode <- system2(command=biofamPath, args=shQuote(argv), wait=T)
  if (exitcode != 0) {
    stop(paste("biofam command failed with exit code", exitcode))
  }
  
  # Load trained model
  object <- loadModel(DirOptions$outFile, object)
  
  return(object)
}
