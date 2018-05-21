
############################################
## Functions to load a trained BioFAModel ##
############################################


#' @title loading a trained BioFAModel
#' @name load_model
#' @description Method to load a trained BioFAModel \cr
#' The training of biofam is done using a Python framework, and the model output is saved as an .hdf5 file, which has to be loaded in the R package.
#' @param file an hdf5 file saved by the biofam Python framework
#' @param object either NULL (default) or an an existing untrained biofam object. If NULL, the \code{\link{BioFAModel}} object is created from the scratch.
#' @param sortFactors boolean inOdicating whether factors should be sorted by variance explained (default is TRUE)
#' @return a \code{\link{BioFAModel}} model
#' @importFrom rhdf5 h5read
#' @export

load_model <- function(file, object = NULL, sortFactors = TRUE, multiView = NULL, multiGroup = NULL) {
  
  # message(paste0("Loading the following BioFAModel: ", file))
  
  if (is.null(object)) object <- new("BioFAModel")
  
  # if(.hasSlot(object,"Status") & length(object@Status) !=0)
  #   if (object@Status == "trained") warning("The specified object is already trained, over-writing training output with new results!")
  
  # Load expectations
  object@expectations <- h5read(file, "expectations")
  tryCatch(object@parameters <- h5read(file, "parameters"), error = function(e) { print(paste("No parameters found in ", file)) })
  object@status <- "trained"

  # Unify names of the nodes with sparsity
  if ("SZ" %in% names(object@expectations)) {
    names(object@Expectations)[which(names(object@expectations)=="SZ")] <- "Z"
  }
  if ("SW" %in% names(object@expectations)) {
    names(object@expectations)[which(names(object@expectations)=="SW")] <- "W"
  }
  
  
  # Load training statistics
  tryCatch( {
    object@training_stats <- h5read(file, 'training_stats', read.attributes=T);
    colnames(object@training_stats$elbo_terms) <- attr(h5read(file,"training_stats/elbo_terms", read.attributes=T),"colnames")
  }, error = function(x) { print("Training stats not found, not loading it...") })

  
  # Load training options
  if (length(object@training_options) == 0) {
    tryCatch(object@training_options <- as.list(h5read(file, 'training_options', read.attributes=T)), error = function(x) { print("Training opts not found, not loading it...") })
  }
    
  # Load model options
  # COMMENTED BECAUSE We always need to load the model options, as h5py sort the views alphabetically
  # if (length(object@ModelOptions) == 0) {
  #   tryCatch(object@ModelOptions <- as.list(h5read(file, 'model_opts',read.attributes=T)), error = function(x) { print("Model opts not found, not loading it...") })
  # }
  tryCatch(object@model_options <- as.list(h5read(file, 'model_options', read.attributes=T)), error = function(x) { print("Model opts not found, not loading it...") })
  for (opt in names(object@model_options)) {
    if (object@model_options[opt] == "False" | object@model_options[opt] == "True") {
      object@model_options[opt] <- as.logical(object@model_options[opt])
    } else {
      object@model_options[opt] <- object@model_options[opt]
    }
  }
  names(object@model_options) <- sapply(names(object@model_options), function(e)
    paste0(sapply(strsplit(e, split = "_")[[1]], function(s) 
      paste(toupper(substring(s, 1, 1)), substring(s, 2), sep="", collapse=" ")), collapse="")
  )

  
  # Load training data

  training_data <- h5read(file, "data")
  feature_data  <- h5read(file, "features")
  sample_data   <- h5read(file, "samples")


  # Define data options

  object@data_options <- list()

  tryCatch( {
    # Specify if multiple views are present
    if (is.null(multi_view)) {
      object@data_options$multi_view <- is.list(feature_data)
    } else {
      object@data_options$multi_view <- multi_view
    }

    # Specify if multiple groups are present
    if (is.null(multi_group)) {
      object@data_options$multi_group <- is.list(sample_data)
    } else {
      object@data_options$multi_group <- multi_group
    }
  }, error = function(x) { print("Error defining data options...") })

  # To keep reverse-compatibility with models without groups (e.g. MOFA models)
  if (!object@data_options$multi_group) {
    sample_data <- list("T1" = sample_data)
    # Samples should be in columns
    if (length(unique(sapply(training_data, nrow))) == 1 && length(unique(sapply(training_data, ncol))) > 1) {
      training_data <- lapply(training_data, t)
    }
    if ((ncol(training_data[[1]]) != length(sample_data[[1]])) & (nrow(training_data[[1]]) == length(sample_data[[1]]))) {
      training_data <- lapply(training_data, t)
    }
    training_data  <- lapply(training_data, function(e) list("B1" = e))
    # Expectations for multi-group and multi-view nodes should be nested lists
    if (("E" %in% names(object@expectations$Y[[1]]))    & (class(object@expectations$Y[[1]]$E) != "list") | 
        (!("E" %in% names(object@expectations$Y[[1]]))) & (class(object@expectations$Y[[1]])   != "list")) {
      tmp <- lapply(names(training_data), function(m) {
        list("B1" = object@expectations$Y[[m]])
      })
      names(tmp) <- names(training_data)
      object@expectations$Y <- tmp
    }
    if (("E" %in% names(object@expectations$Z))    & (class(object@expectations$Z$E) != "list") | 
        (!("E" %in% names(object@expectations$Z))) & (class(object@expectations$Z)   != "list")) {
      object@expectations$Z <- list("B1" = object@expectations$Z)
    }
    object@data_options$multi_group <- TRUE
  }

  # To keep reverse-compatibility with models having views as groups (e.g. transposed MOFA)
  if (!object@data_options$multi_view) {
    feature_data <- list("V1" = feature_data)
    # Features should be in rows
    if (length(unique(sapply(training_data, ncol))) == 1 && length(unique(sapply(training_data, nrow))) > 1) {
      training_data <- lapply(training_data, t)
    }
    training_data <- list("V1" = training_data)
    # Expectations for multi-group and multi-view nodes should be nested lists
    if (("E" %in% names(object@expectations$Y[[1]]))    & (class(object@expectations$Y[[1]]$E) != "list") | 
        (!("E" %in% names(object@expectations$Y[[1]]))) & (class(object@expectations$Y[[1]])   != "list")) {
      object@Expectations$Y <- list("V1" = object@Expectations$Y)
    }
    if (("E" %in% names(object@expectations$W[[1]]))    & (class(object@expectations$W[[1]]$E) != "list") | 
        (!("E" %in% names(object@expectations$W[[1]]))) & (class(object@expectations$W[[1]])   != "list")) {
      object@expectations$W <- list("V1" = object@expectations$W)
    }
    object@data_options$multi_view <- TRUE
  }


  # Load dimensions
  object@dimensions[["M"]] <- length(training_data)
  object@dimensions[["P"]] <- length(training_data[[1]])
  object@dimensions[["N"]] <- sapply(training_data[[1]], ncol)
  object@dimensions[["D"]] <- sapply(training_data, function(e) nrow(e[[1]]))
  # K=tail(training_stats$activeK[!is.nan(training_stats$activeK)],n=1)
  object@dimensions[["K"]] <- ncol(object@expectations$Z[[1]]$E)


  # Fix sample and feature names is they are null
  if (is.null(sample_data)) {
    sample_data <- paste0("S", lapply(object@dimensions[["N"]], function(n) as.character(1:n)))
  }
  if (is.null(feature_data)) {
    feature_data <- paste0("G", lapply(object@dimensions[["D"]], function(d) as.character(1:d)))
  }  
  
  # Give corresponding names for rows (features) and columns (samples)
  tryCatch( {
    for (m in names(training_data)) {  # there is always at least one view
      for (p in names(training_data[[m]])) {  # there is always at least one group
        rownames(training_data[[m]][[p]]) <- feature_data[[m]]
        colnames(training_data[[m]][[p]]) <- sample_data[[p]]
      }
    }
    object@training_data <- training_data
  }, error = function(x) { cat("Error defining feature and sample names!..\n") })
  
  # Replace NaN by NA
  for (m in names(training_data)) {
    for (h in names(training_data[[m]])) {
      # object@Expectations[[m]][[h]][is.nan(object@Expectations[[m]][[h]])] <- NA
      training_data[[m]][[h]][is.nan(training_data[[m]][[h]])] <- NA
    }
  }
  
  # Sanity check on the order of the likelihoods
  if (!is.null(attr(training_data, "Likelihood"))) {
    lik <- attr(training_data, "Likelihood")
    if (!all(object@model_options$likelihood == lik)) {
      object@model_options$likelihood <- lik
      names(object@model_options$likelihood) <- names(training_data)
    }
  }
  
  # Set view and group names
  if (is.null(names(object@training_data))) {
    view_names(object) <- paste0("V", as.character(1:object@Dimensions[["M"]]))
  } else {
    view_names(object) <- names(object@training_data)
  }
  
  if (is.null(names(object@training_data[[1]]))) {
    group_names(object) <- paste0("T", as.character(1:object@dimensions[["P"]]))
  } else {
    group_names(object) <- names(object@training_data[[1]])
  }
  
  # Update old models
  object <- .update_old_model(object)

  # Set sample, feature, and factor names
  sample_names(object)  <- lapply(object@training_data[[1]], colnames)
  feature_names(object) <- lapply(object@training_data, function(e) rownames(e[[1]]))
  factor_names(object)  <- paste0("F", as.character(1:object@dimensions[["K"]]))
  
  # Add names to likelihood vector
  if (!is.null(names(object@model_options$likelihood))) {
    names(object@model_options$likelihood) <- view_names(object)
  }
  
  # Rename covariates, including intercept
  # if (object@ModelOptions$LearnIntercept == TRUE) factorNames(object) <- c("intercept",as.character(1:(object@Dimensions[["K"]]-1)))
  # if (!is.null(object@ModelOptions$covariates)) {
  #   if (object@ModelOptions$LearnIntercept == TRUE) {
  #     factorNames(object) <- c("intercept", colnames(object@ModelOptions$covariates), as.character((ncol(object@ModelOptions$covariates)+1:(object@Dimensions[["K"]]-1-ncol(object@ModelOptions$covariates)))))
  #   } else {
  #     factorNames(object) <- c(colnames(object@ModelOptions$covariates), as.character((ncol(object@ModelOptions$covariates)+1:(object@Dimensions[["K"]]-1))))
  #   }
  # }
  
  # Rename factors if intercept is included
  if (object@model_options$learn_intercept) {
    intercept_idx <- names(which(sapply(apply(object@expectations$Z, 2, unique),length)==1))
    factornames <- as.character(1:(object@dimensions[["K"]]))
    factornames[factornames==intercept_idx] <- "intercept"
    factor_names(object) <- factornames
    # object@Dimensions[["K"]] <- object@Dimensions[["K"]] - 1
  }
  # if (!is.null(object@ModelOptions$covariates)) {
  #   stop("Covariates not working")
  # }
  
  # Parse factors: Mask passenger samples
  object <- detect_passengers(object)

  # Parse factors: order factors in order of variance explained
  if (sort_factors) {
    r2 <- rowSums(sapply(calculate_variance_explained(object)$r2_per_factor, function(e) rowSums(e)))
    order_factors <- c(names(r2)[order(r2, decreasing = T)])
    if (object@model_options$learn_intercept) { order_factors <- c("intercept", order_factors) }
    object <- subsetFactors(object, order_factors)
    if (object@model_options$learn_intercept) { 
     factorNames(object) <- c("intercept", 1:(object@dimensions$K-1))
    } else {
     factorNames(object) <- c(1:object@dimensions$K) 
    }
  }
  
  
  # Check for intercept factors
  # findInterceptFactors(object)
  
  # Do quality control on the model
  # qualityControl(object)

  return(object)
}

