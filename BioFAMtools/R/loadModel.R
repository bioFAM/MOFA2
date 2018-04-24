
############################################
## Functions to load a trained BioFAModel ##
############################################


#' @title loading a trained BioFAModel
#' @name loadModel
#' @description Method to load a trained BioFAModel. \cr
#' The training of BioFAM is done using a Python framework, and the model output is saved as an .hdf5 file, which has to be loaded in the R package.
#' @param file an hdf5 file saved by the BioFAM python framework.
#' @param object either NULL (default) or an an existing untrained BioFAM object. If NULL, the \code{\link{BioFAModel}} object is created from the scratch.
#' @param sortFactors boolean indicating whether factors should be sorted by variance explained (default is TRUE)
#' @return a \code{\link{BioFAModel}} model.
#' @importFrom rhdf5 h5read
#' @export

loadModel <- function(file, object = NULL, sortFactors = TRUE, multiView = NULL, multiBatch = NULL) {
  
  # message(paste0("Loading the following BioFAModel: ", file))
  
  if (is.null(object)) object <- new("BioFAModel")
  
  # if(.hasSlot(object,"Status") & length(object@Status) !=0)
  #   if (object@Status == "trained") warning("The specified object is already trained, over-writing training output with new results!")
  
  # Load expectations
  object@Expectations <- h5read(file, "expectations")
  tryCatch(object@Parameters <- h5read(file, "parameters"), error = function(e) { print(paste("No parameters found in ", file)) })
  object@Status <- "trained"

  # Unify names of the nodes with sparsity
  if ("SZ" %in% names(object@Expectations)) {
    names(object@Expectations)[which(names(object@Expectations)=="SZ")] <- "Z"
  }
  if ("SW" %in% names(object@Expectations)) {
    names(object@Expectations)[which(names(object@Expectations)=="SW")] <- "W"
  }
  
  
  # Load training statistics
  tryCatch( {
    object@TrainStats <- h5read(file, 'training_stats', read.attributes=T);
    colnames(object@TrainStats$elbo_terms) <- attr(h5read(file,"training_stats/elbo_terms", read.attributes=T),"colnames")
  }, error = function(x) { print("Training stats not found, not loading it...") })

  
  # Load training options
  if (length(object@TrainOpts) == 0) {
    tryCatch(object@TrainOpts <- as.list(h5read(file, 'training_opts', read.attributes=T)), error = function(x) { print("Training opts not found, not loading it...") })
  }
    
  # Load model options
  # COMMENTED BECAUSE We always need to load the model options, as h5py sort the views alphabetically
  # if (length(object@ModelOpts) == 0) {
  #   tryCatch(object@ModelOpts <- as.list(h5read(file, 'model_opts',read.attributes=T)), error = function(x) { print("Model opts not found, not loading it...") })
  # }
  tryCatch(object@ModelOpts <- as.list(h5read(file, 'model_opts', read.attributes=T)), error = function(x) { print("Model opts not found, not loading it...") })
  for (opt in names(object@ModelOpts)) {
    if (opt == "False" | opt == "True") {
      object@ModelOpts[opt] <- as.logical(object@ModelOpts[opt])
    }
  }

  
  # Load training data

  TrainData   <- h5read(file, "data")
  featureData <- h5read(file, "features")
  sampleData  <- h5read(file, "samples")


  # Define data options

  object@DataOpts <- list()

  tryCatch( {
    # Specify if multiple views are present
    if (is.null(multiView)) {
      object@DataOpts$multiView <- is.list(featureData)
    } else {
      object@DataOpts$multiView <- multiView
    }

    # Specify if multiple batches are present
    if (is.null(multiBatch)) {
      object@DataOpts$multiBatch <- is.list(sampleData)
    } else {
      object@DataOpts$multiBatch <- multiBatch
    }
  }, error = function(x) { print("Error defining data options...") })

  # To keep reverse-compatibility with models without batches (e.g. MOFA models)
  if (!object@DataOpts$multiBatch) {
    sampleData <- list("B1" = sampleData)
    # Samples should be in columns
    if (length(unique(sapply(TrainData, nrow))) == 1 && length(unique(sapply(TrainData, ncol))) > 1) {
      TrainData <- lapply(TrainData, t)
    }
    if ((ncol(TrainData[[1]]) != length(sampleData[[1]])) & (nrow(TrainData[[1]]) == length(sampleData[[1]]))) {
      TrainData <- lapply(TrainData, t)
    }
    TrainData  <- lapply(TrainData, function(e) list("B1" = e))
    # Expectations for multi-batch and multi-view nodes should be nested lists
    if (("E" %in% names(object@Expectations$Y[[1]]))    & (class(object@Expectations$Y[[1]]$E) != "list") | 
        (!("E" %in% names(object@Expectations$Y[[1]]))) & (class(object@Expectations$Y[[1]])   != "list")) {
      tmp <- lapply(names(TrainData), function(m) {
        list("B1" = object@Expectations$Y[[m]])
      })
      names(tmp) <- names(TrainData)
      object@Expectations$Y <- tmp
    }
    if (("E" %in% names(object@Expectations$Z))    & (class(object@Expectations$Z$E) != "list") | 
        (!("E" %in% names(object@Expectations$Z))) & (class(object@Expectations$Z)   != "list")) {
      object@Expectations$Z <- list("B1" = object@Expectations$Z)
    }
    object@DataOpts$multiBatch <- TRUE
  }

  # To keep reverse-compatibility with models having views as batches (e.g. transposed MOFA)
  if (!object@DataOpts$multiView) {
    featureData <- list("V1" = featureData)
    # Features should be in rows
    if (length(unique(sapply(TrainData, ncol))) == 1 && length(unique(sapply(TrainData, nrow))) > 1) {
      TrainData <- lapply(TrainData, t)
    }
    TrainData <- list("V1" = TrainData)
    # Expectations for multi-batch and multi-view nodes should be nested lists
    if (("E" %in% names(object@Expectations$Y[[1]]))    & (class(object@Expectations$Y[[1]]$E) != "list") | 
        (!("E" %in% names(object@Expectations$Y[[1]]))) & (class(object@Expectations$Y[[1]])   != "list")) {
      object@Expectations$Y <- list("V1" = object@Expectations$Y)
    }
    if (("E" %in% names(object@Expectations$W[[1]]))    & (class(object@Expectations$W[[1]]$E) != "list") | 
        (!("E" %in% names(object@Expectations$W[[1]]))) & (class(object@Expectations$W[[1]])   != "list")) {
      object@Expectations$W <- list("V1" = object@Expectations$W)
    }
    object@DataOpts$multiView <- TRUE
  }


  # Load dimensions
  object@Dimensions[["M"]] <- length(TrainData)
  object@Dimensions[["H"]] <- length(TrainData[[1]])
  object@Dimensions[["N"]] <- sapply(TrainData[[1]], ncol)
  object@Dimensions[["D"]] <- sapply(TrainData, function(e) nrow(e[[1]]))
  # K=tail(training_stats$activeK[!is.nan(training_stats$activeK)],n=1)
  object@Dimensions[["K"]] <- ncol(object@Expectations$Z[[1]]$E)


  # Fix sample and feature names is they are null
  if (is.null(sampleData)) {
    sampleData <- paste0("S", lapply(object@Dimensions[["N"]], function(n) as.character(1:n)))
  }
  if (is.null(featureData)) {
    featureData <- paste0("G", lapply(object@Dimensions[["D"]], function(d) as.character(1:d)))
  }  
  
  # Give corresponding names for rows (features) and columns (samples)
  tryCatch( {
    for (m in names(TrainData)) {  # there is always at least one view
      for (h in names(TrainData[[m]])) {  # there is always at least one batch
        rownames(TrainData[[m]][[h]]) <- featureData[[m]]
        colnames(TrainData[[m]][[h]]) <- sampleData[[h]]
      }
    }
    object@TrainData <- TrainData
  }, error = function(x) { cat("Error defining feature and sample names!..\n") })
  
  # Replace NaN by NA
  for (m in names(TrainData)) {
    for (h in names(TrainData[[m]])) {
      # object@Expectations[[m]][[h]][is.nan(object@Expectations[[m]][[h]])] <- NA
      TrainData[[m]][[h]][is.nan(TrainData[[m]][[h]])] <- NA
    }
  }
  
  # Sanity check on the order of the likelihoods
  if (!is.null(attr(TrainData, "likelihood"))) {
    lik <- attr(TrainData, "likelihood")
    if (!all(object@ModelOpts$likelihood == lik)) {
      object@ModelOpts$likelihood <- lik
      names(object@ModelOpts$likelihood) <- names(TrainData)
    }
  }
  
  # Set view and batch names
  if (is.null(names(object@TrainData))) {
    viewNames(object) <- paste0("V", as.character(1:object@Dimensions[["M"]]))
  } else {
    viewNames(object) <- names(object@TrainData)
  }
  
  if (is.null(names(object@TrainData[[1]]))) {
    batchNames(object) <- paste0("B", as.character(1:object@Dimensions[["H"]]))
  } else {
    batchNames(object) <- names(object@TrainData[[1]])
  }
  
  # Update old models
  object <- .updateOldModel(object)

  # Set sample, feature, and factor names
  sampleNames(object)  <- lapply(object@TrainData[[1]], colnames)
  featureNames(object) <- lapply(object@TrainData, function(e) rownames(e[[1]]))
  factorNames(object)  <- paste0("F", as.character(1:object@Dimensions[["K"]]))
  
  # Add names to likelihood vector
  if (!is.null(names(object@ModelOpts$likelihood))) {
    names(object@ModelOpts$likelihood) <- viewNames(object)
  }
  
  # Rename covariates, including intercept
  # if (object@ModelOpts$learnIntercept == TRUE) factorNames(object) <- c("intercept",as.character(1:(object@Dimensions[["K"]]-1)))
  # if (!is.null(object@ModelOpts$covariates)) {
  #   if (object@ModelOpts$learnIntercept == TRUE) {
  #     factorNames(object) <- c("intercept", colnames(object@ModelOpts$covariates), as.character((ncol(object@ModelOpts$covariates)+1:(object@Dimensions[["K"]]-1-ncol(object@ModelOpts$covariates)))))
  #   } else {
  #     factorNames(object) <- c(colnames(object@ModelOpts$covariates), as.character((ncol(object@ModelOpts$covariates)+1:(object@Dimensions[["K"]]-1))))
  #   }
  # }
  
  # Rename factors if intercept is included
  if (object@ModelOpts$learnIntercept) {
    intercept_idx <- names(which(sapply(apply(object@Expectations$Z, 2, unique),length)==1))
    factornames <- as.character(1:(object@Dimensions[["K"]]))
    factornames[factornames==intercept_idx] <- "intercept"
    factorNames(object) <- factornames
    # object@Dimensions[["K"]] <- object@Dimensions[["K"]] - 1
  }
  # if (!is.null(object@ModelOpts$covariates)) {
  #   stop("Covariates not working")
  # }
  
  # Parse factors: Mask passenger samples
  object <- detectPassengers(object)

  # Parse factors: order factors in order of variance explained
  if (sortFactors) {
   r2 <- rowSums(calculateVarianceExplained(object)$R2PerFactor)
   order_factors <- c(names(r2)[order(r2, decreasing = T)])
   if (object@ModelOpts$learnIntercept==T) { order_factors <- c("intercept",order_factors) }
   object <- subsetFactors(object,order_factors)
   if (object@ModelOpts$learnIntercept==T) { 
     factorNames(object) <- c("intercept",1:(object@Dimensions$K-1))
   } else {
     factorNames(object) <- c(1:object@Dimensions$K) 
   }
  }
  
  
  # Check for intercept factors
  # findInterceptFactors(object)
  
  # Do quality control on the model
  # qualityControl(object)

  return(object)
}

