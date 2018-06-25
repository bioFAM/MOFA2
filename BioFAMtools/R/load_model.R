
############################################
## Functions to load a trained BioFAModel ##
############################################


#' @title loading a trained BioFAModel
#' @name load_model
#' @description TO-FILL
#' @param file TO-FILL
#' @param object TO-FILL
#' @param sort_factors TO-FILL
#' @return a \code{\link{BioFAModel}} model
#' @importFrom rhdf5 h5read
#' @export

load_model <- function(file, object = NULL, sort_factors = TRUE) {
  
  if (is.null(object)) object <- new("BioFAModel")
  
  if(.hasSlot(object, "status") & length(object@status) != 0)
    if (object@status == "trained") warning("The specified object is already trained, over-writing training output with new results.")
  

  # Load data

  # Load training data and the identity of features and samples
  training_data <- h5read(file, "data")
  features_data <- h5read(file, "features")
  samples_data  <- h5read(file, "samples")

  # Replace NaN by NA
  for (m in 1:length(training_data)) {
    for (h in 1:length(training_data[[m]])) {
      training_data[[m]][[h]][is.nan(training_data[[m]][[h]])] <- NA
    }
  }
  object@training_data <- training_data

  # Fix sample and feature names is they are null
  if (is.null(samples_data)) {
    samples_data <- paste0("S", lapply(object@dimensions[["N"]], function(n) as.character(1:n)))
  }
  if (is.null(features_data)) {
    features_data <- paste0("G", lapply(object@dimensions[["D"]], function(d) as.character(1:d)))
  }

  # Give corresponding names for rows (features) and columns (samples)
  tryCatch( {
    for (m in 1:length(training_data)) {  # there is always at least one view
      for (p in 1:length(training_data[[m]])) {  # there is always at least one group
        rownames(training_data[[m]][[p]]) <- features_data[[m]]
        colnames(training_data[[m]][[p]]) <- samples_data[[p]]
      }
    }
    object@training_data <- training_data
  }, error = function(x) { cat("Error defining feature and sample names!..\n") })


  # Load expectations
  object@expectations <- h5read(file, "expectations")
  tryCatch(object@parameters <- h5read(file, "parameters"), error = function(e) { print(paste("No parameters found in ", file)) })

  # Unify names of the nodes with sparsity
  if ("SZ" %in% names(object@expectations)) {
    names(object@expectations)[which(names(object@expectations)=="SZ")] <- "Z"
  }
  if ("SW" %in% names(object@expectations)) {
    names(object@expectations)[which(names(object@expectations)=="SW")] <- "W"
  }


  # Specify dimensionality of the data
  object@dimensions[["M"]] <- length(training_data)                           # number of views (groups of features)
  object@dimensions[["P"]] <- length(training_data[[1]])                      # number of groups (groups of samples)
  object@dimensions[["N"]] <- sapply(training_data[[1]], ncol)                # number of samples in every group
  object@dimensions[["D"]] <- sapply(training_data, function(e) nrow(e[[1]])) # number of features in every view
  object@dimensions[["K"]] <- ncol(object@expectations$W[[1]]$E)              # number of factors



  # Set view and group names
  if (is.null(names(object@training_data))) {
    views_names(object) <- paste0("V", as.character(1:object@dimensions[["M"]]))
  } else {
    views_names(object) <- names(object@training_data)
  }

  if (is.null(names(object@training_data[[1]]))) {
    groups_names(object) <- paste0("T", as.character(1:object@dimensions[["P"]]))
  } else {
    groups_names(object) <- names(object@training_data[[1]])
  }


  # Load model options
  
  tryCatch( {
    object@model_options <- as.list(h5read(file, 'model_options', read.attributes=T))
  }, error = function(x) { print("Model opts not found, not loading it...") })

  for (opt in names(object@model_options)) {
    if (object@model_options[opt] == "False" | object@model_options[opt] == "True") {
      object@model_options[opt] <- as.logical(object@model_options[opt])
    } else {
      object@model_options[opt] <- object@model_options[opt]
    }
  }


  # Define node types of the model
  object@model_options$nodes <- list(multiview_nodes  = c("W", "AlphaW", "ThetaW", "SigmaAlphaW"),
                                     multigroup_nodes = c("Z", "AlphaZ", "ThetaZ", "SigmaZ"),
                                     twodim_nodes     = c("Y", "Tau"))


  # Set sample, feature, and factor names for the data and all the expectations
  # NOTE: The side effect is also removing the extra nestedness level (exp$E -> exp)
  samples_names(object)  <- samples_data
  features_names(object) <- features_data
  factors_names(object)  <- paste0("F", as.character(1:object@dimensions[["K"]]))


  # Define data options


  # Sanity check on the order of the likelihoods
  if (!is.null(attr(training_data, "likelihood"))) {
    lik <- attr(training_data, "likelihood")
    if (!all(object@model_options$likelihood == lik)) {
      object@model_options$likelihood <- lik
      names(object@model_options$likelihood) <- names(training_data)
    }
  }

  # Add names to likelihood vector
  if (!is.null(names(object@model_options$likelihood))) {
    names(object@model_options$likelihood) <- views_names(object)
  }

  # Rename factors if intercept is included
  if (object@model_options$learn_intercept) {
    intercept_idx <- names(which(sapply(apply(object@expectations$Z, 2, unique),length)==1))
    new_factornames <- as.character(1:(object@dimensions[["K"]]))
    new_factornames[new_factornames == intercept_idx] <- "intercept"
    factors_names(object) <- new_factornames
    # object@Dimensions[["K"]] <- object@Dimensions[["K"]] - 1
  }

  # Parse factors: order factors in order of variance explained
  if (sort_factors) {
    r2 <- rowSums(sapply(calculate_variance_explained(object)$r2_per_factor, function(e) rowSums(e)))
    order_factors <- c(names(r2)[order(r2, decreasing = T)])
    if (object@model_options$learn_intercept) { order_factors <- c("intercept", order_factors) }
    object <- subset_factors(object, order_factors)
    print(order_factors)
    if (object@model_options$learn_intercept) { 
     factors_names(object) <- c("intercept", 1:(object@dimensions$K-1))
    } else {
     factors_names(object) <- c(1:object@dimensions$K) 
    }
  }
  


  # Load training options
  if (length(object@training_options) == 0) {
    tryCatch( {
      object@training_options <- as.list(h5read(file, 'training_opts', read.attributes=T))
    }, error = function(x) { print("Training opts not found, not loading it...") })
  }    


  # Load training statistics
  tryCatch( {
    object@training_stats <- h5read(file, 'training_stats', read.attributes=T)
    colnames(object@training_stats$elbo_terms) <- attr(h5read(file,"training_stats/elbo_terms", read.attributes=T),"colnames")
  }, error = function(x) { print("Training stats not found, not loading it...") })

  # Update old models
  object <- .update_old_model(object)
  
  # Rename covariates, including intercept
  # if (object@ModelOptions$LearnIntercept == TRUE) factorNames(object) <- c("intercept",as.character(1:(object@Dimensions[["K"]]-1)))
  # if (!is.null(object@ModelOptions$covariates)) {
  #   if (object@ModelOptions$LearnIntercept == TRUE) {
  #     factorNames(object) <- c("intercept", colnames(object@ModelOptions$covariates), as.character((ncol(object@ModelOptions$covariates)+1:(object@Dimensions[["K"]]-1-ncol(object@ModelOptions$covariates)))))
  #   } else {
  #     factorNames(object) <- c(colnames(object@ModelOptions$covariates), as.character((ncol(object@ModelOptions$covariates)+1:(object@Dimensions[["K"]]-1))))
  #   }
  # }
  
  
  # Parse factors: Mask passenger samples
  # object <- detect_passengers(object)
  
  # Do quality control on the model
  # qualityControl(object)

  return(object)
}

