
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

loadModel <- function(file, object = NULL, sortFactors = TRUE, sharedFeatures = NULL) {
  
  # message(paste0("Loading the following BioFAModel: ", file))
  
  if (is.null(object)) object <- new("BioFAModel")
  
  # if(.hasSlot(object,"Status") & length(object@Status) !=0)
  #   if (object@Status == "trained") warning("The specified object is already trained, over-writing training output with new results!")
  
  # Load expectations
  object@Expectations <- h5read(file, "expectations")
  object@Parameters <- h5read(file,"parameters")
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

  tryCatch( {
    TrainData   <- h5read(file, "data")
    featureData <- h5read(file, "features")
    sampleData  <- h5read(file, "samples")

    # Define data options (if features are in rows or in columns)
    if (is.null(sharedFeatures)) {
      if (is.list(sampleData) & (length(featureData) == nrow(TrainData[[1]]))) {
        object@DataOpts <- list(shared_features = TRUE)
      } else {
        object@DataOpts <- list(shared_features = FALSE)
      }
    } else {
      object@DataOpts <- list(shared_features = sharedFeatures)
    }
  }, error = function(x) { print("Error loading the training data...") })


  tryCatch( {
    
    if (object@DataOpts$shared_features) {
      for (m in names(TrainData)) {
        rownames(TrainData[[m]]) <- featureData
        colnames(TrainData[[m]]) <- sampleData[[m]]
      }
    } else {
      for (m in names(TrainData)) {
        rownames(TrainData[[m]]) <- sampleData
        colnames(TrainData[[m]]) <- featureData[[m]]
      }
    }
    TrainData <- lapply(TrainData, t)
    object@TrainData <- TrainData

  }, error = function(x) { cat(paste0("Error loading the training data!.. Try to load the model with sharedFeatures set to ", as.character(!object@DataOpts$shared_features), ".\n")) })
  
  # Replace NaN by NA
  for (m in names(TrainData)) {
    # object@Expectations[[m]][is.nan(object@Expectations[[m]])] <- NA
    TrainData[[m]][is.nan(TrainData[[m]])] <- NA
  }
  
  # Sanity check on the order of the likelihoods
  if (!is.null(attr(TrainData,"likelihood"))) {
    lik <- attr(TrainData,"likelihood")
    if (!all(object@ModelOpts$likelihood == lik)) {
      object@ModelOpts$likelihood <- lik
      names(object@ModelOpts$likelihood) <- names(TrainData)
    }
  }
  
  # Update old models
  object <- .updateOldModel(object)
  
  # Load dimensions
  object@Dimensions[["M"]] <- length(object@TrainData)
  object@Dimensions[["N"]] <- ncol(object@TrainData[[1]])
  object@Dimensions[["D"]] <- sapply(object@TrainData, nrow)
  # K=tail(training_stats$activeK[!is.nan(training_stats$activeK)],n=1)
  object@Dimensions[["K"]] <- ncol(object@Expectations$Z)
  
  # Set view, sample, feature and factor names
  if (is.null(names(object@TrainData))) {
    viewNames(object) <- paste0("V", as.character(1:object@Dimensions[["M"]]))
  } #else {
    #viewNames(object) <- names(object@TrainData)
  #}
  if (is.null(colnames(object@TrainData[[1]]))) {
    sampleNames(object) <- paste0("S", as.character(1:object@Dimensions[["N"]]))
  } else {
    sampleNames(object) <- colnames(object@TrainData[[1]])
  }
  featureNames(object) <- lapply(object@TrainData, rownames)
  factorNames(object) <- as.character(1:object@Dimensions[["K"]])
  
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
  if (object@ModelOpts$learnIntercept == TRUE) {
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
  # if (sortFactors == T) {
  #   r2 <- rowSums(calculateVarianceExplained(object)$R2PerFactor)
  #   order_factors <- c(names(r2)[order(r2, decreasing = T)])
  #   if (object@ModelOpts$learnIntercept==T) { order_factors <- c("intercept",order_factors) }
  #   object <- subsetFactors(object,order_factors)
  #   if (object@ModelOpts$learnIntercept==T) { 
  #     factorNames(object) <- c("intercept",1:(object@Dimensions$K-1))
  #   } else {
  #     factorNames(object) <- c(1:object@Dimensions$K) 
  #   }
  # }
  
  
  # Check for intercept factors
  # findInterceptFactors(object)
  
  # Do quality control on the model
  # qualityControl(object)
  
  return(object)
}

