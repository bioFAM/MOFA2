
################################
## Functions to do subsetting ##
################################

#' @title Subset groups
#' @name subset_groups
#' @description Method to subset (or sort) groups
#' @param object a \code{\link{MOFA}} object.
#' @param groups character vector with the groups names, numeric vector with the groups indices
#' or logical vector with the groups to be kept as TRUE.
#' @return A \code{\link{MOFA}} object
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Subset the first group
#' model <- subset_groups(model, groups = 1)
subset_groups <- function(object, groups) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(length(groups) <= object@dimensions[["G"]])
  
  # Define groups    
  groups <- .check_and_get_groups(object, groups)
  
  # Subset expectations
  if (length(object@expectations)>0) {
    if ("Z" %in% names(object@expectations) & length(object@expectations$Z)>0)
      object@expectations$Z <- object@expectations$Z[groups]
    if ("Y" %in% names(object@expectations) & length(object@expectations$Y)>0)
      object@expectations$Y <- sapply(object@expectations$Y, function(x) x[groups], simplify = FALSE, USE.NAMES = TRUE) 
  }
  
  # Subset data
  if (length(object@data)>0) {
    object@data <- sapply(object@data, function(x) x[groups], simplify = FALSE, USE.NAMES = TRUE) 
  }

  # Subset imputed data
  if (length(object@imputed_data)>0) {
    object@imputed_data <- sapply(object@imputed_data, function(x) x[groups], simplify = FALSE, USE.NAMES = TRUE) 
  }
  
  # Subset intercepts
  if (length(object@intercepts[[1]])>0) {
    object@intercepts <- sapply(object@intercepts, function(x) x[groups], simplify = FALSE, USE.NAMES = TRUE) 
  }
  
  # Update dimensionality
  object@dimensions[["G"]] <- length(groups)
  object@dimensions[["N"]] <- object@dimensions[["N"]][groups]
  
  # Subset sample metadata
  stopifnot(groups%in%unique(object@samples_metadata$group))
  object@samples_metadata <- object@samples_metadata[object@samples_metadata$group %in% groups,]
  object@samples_metadata$group <- factor(object@samples_metadata$group, levels=groups)
  
  # Re-order samples
  samples <- unname(unlist(lapply(object@data[[1]],colnames)))
  object@samples_metadata <- object@samples_metadata[match(samples, object@samples_metadata$sample),]
  
  # Sanity checks
  stopifnot(object@samples_metadata$sample == unlist(lapply(object@data[[1]],colnames)))
  stopifnot(object@samples_metadata$sample == unlist(lapply(object@expectations$Z,rownames)))
  stopifnot(object@samples_metadata$sample == unlist(lapply(object@expectations$Y,colnames)))
  
  # Update groups names
  # groups_names(object) <- groups # don't need to run this
  object@data_options$groups <- groups
  
  # Subset variance explained
  object@cache[["variance_explained"]]$r2_per_factor <- object@cache[["variance_explained"]]$r2_per_factor[groups]
  object@cache[["variance_explained"]]$r2_total <- object@cache[["variance_explained"]]$r2_total[groups]
  
  return(object)
}


#' @title Subset views
#' @name subset_views
#' @description Method to subset (or sort) views
#' @param object a \code{\link{MOFA}} object.
#' @param views character vector with the views names, numeric vector with the views indices,
#' or logical vector with the views to be kept as TRUE.
#' @return A \code{\link{MOFA}} object
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Subset the first view
#' model <- subset_views(model, views = 1)
subset_views <- function(object, views) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(length(views) <= object@dimensions[["M"]])
  # warning("Removing views a posteriori is fine for an exploratory analysis, but you should removing them before training!")
  
  # Define views
  views <- .check_and_get_views(object, views)
  
  # Subset relevant slots
  if (length(object@expectations)>0) {
    object@expectations$W <- object@expectations$W[views]
    object@expectations$Y <- object@expectations$Y[views]
  }

  # Subset data
  if (length(object@data)>0) {
    object@data <- object@data[views]
  }
  
  # Subset imputed data
  if (length(object@imputed_data)>0) {
    object@imputed_data <- object@imputed_data[views]
  }
  
  # Subset intercepts
  if (length(object@intercepts[[1]])>0) {
    object@intercepts <- object@intercepts[views]
  }
  
  # Subset feature metadata
  if (length(object@features_metadata)>0) {
    object@features_metadata <- object@features_metadata[object@features_metadata$view %in% views,]
  }
  
  # Subset likelihoods
  object@model_options$likelihoods <- object@model_options$likelihoods[views]
  # Update dimensionality
  object@dimensions[["M"]] <- length(views)
  object@dimensions[["D"]] <- object@dimensions[["D"]][views]
  
  # Update view names
  object@data_options$views <- views
  
  # Subset variance explained
  if ((methods::.hasSlot(object, "cache")) && ("variance_explained" %in% names(object@cache))) {
    object@cache[["variance_explained"]]$r2_per_factor <- lapply(object@cache[["variance_explained"]]$r2_per_factor, function(x) x[,views,drop=FALSE])
    object@cache[["variance_explained"]]$r2_total <- lapply(object@cache[["variance_explained"]]$r2_total, function(x) x[views])
  }
  
  return(object)
}


#' @title Subset factors
#' @name subset_factors
#' @description Method to subset (or sort) factors
#' @param object a \code{\link{MOFA}} object.
#' @param factors character vector with the factor names, or numeric vector with the index of the factors.
#' @export
#' @return A \code{\link{MOFA}} object
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Subset factors 1 to 3
#' model <- subset_factors(model, factors = 1:3)
subset_factors <- function(object, factors) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(length(factors) <= object@dimensions[["K"]])
  
  # Define factors
  factors <- .check_and_get_factors(object, factors)
  
  # Subset expectations
  nodes_with_factors <- list(nodes = c("Z", "W", "AlphaZ", "AlphaW", "ThetaZ", "ThetaW"), axes = c(2, 2, 0, 0, 0, 0))
  stopifnot(all(nodes_with_factors$axes %in% c(0, 1, 2)))

  if (length(object@expectations)>0) {
    for (i in seq_len(length(nodes_with_factors$nodes))) {
      node <- nodes_with_factors$nodes[i]
      axis <- nodes_with_factors$axes[i]
      if (node %in% names(object@expectations)) {
        if (axis == 1) {
          object@expectations[[node]] <- sapply(object@expectations[[node]], function(x) x[factors,,drop=FALSE], simplify = FALSE, USE.NAMES = TRUE)
        } else if (axis == 2) {
          object@expectations[[node]] <- sapply(object@expectations[[node]], function(x) x[,factors,drop=FALSE], simplify = FALSE, USE.NAMES = TRUE)
        } else {
          object@expectations[[node]] <- sapply(object@expectations[[node]], function(x) x[factors], simplify = FALSE, USE.NAMES = TRUE)
        }
      }
    }
  }
  
  # Subset per-factor variance explained estimates
  if ((methods::.hasSlot(object, "cache")) && ("variance_explained" %in% names(object@cache))) {
    object@cache[["variance_explained"]]$r2_per_factor <- lapply(object@cache[["variance_explained"]]$r2_per_factor, function(x) x[factors,,drop=FALSE])
  }
  
  # Relalculate total variance explained estimates  
  if (length(factors) < object@dimensions[["K"]]) {
    object@cache[["variance_explained"]]$r2_total <- lapply(object@cache[["variance_explained"]]$r2_per_factor, colSums)
  }
  
  # Update dimensionality
  object@dimensions[["K"]] <- length(factors)
  
  # Update factor names
  factors_names(object) <- paste0("Factor", as.character(seq_len(object@dimensions[["K"]])))
  
  
  return(object)
}



#' @title Subset samples
#' @name subset_samples
#' @description Method to subset (or sort) samples
#' @param object a \code{\link{MOFA}} object.
#' @param samples character vector with the sample names or numeric vector with the sample indices.
#' @export
#' @return A \code{\link{MOFA}} object
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # (TO-DO) Remove a specific sample from the model (an outlier)
subset_samples <- function(object, samples) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Define samples
  samples <- .check_and_get_samples(object, samples)
  
  # Check if an entire group needs to be removed
  groups <- as.character(unique(object@samples_metadata[match(samples, object@samples_metadata$sample),]$group))
  if (length(groups)<length(groups_names(object))) object <- subset_groups(object, groups)
  
  # Subset data and expectations
  groups <- groups_names(object)
  tmp <- lapply(groups, function(g) samples_names(object)[[g]][samples_names(object)[[g]] %in% samples])
  names(tmp) <- groups
  
  for (g in groups) {
    samples_g <- tmp[[g]]
    
    # Subset expectations
    if (length(object@expectations)>0) {
      if ("Z" %in% names(object@expectations) & length(object@expectations$Z)>0) {
        object@expectations$Z[[g]] <- object@expectations$Z[[g]][samples_g,, drop=FALSE]
      }
      
      if ("Y" %in% names(object@expectations) & length(object@expectations$Y)>0) {
        for (m in views_names(object)) {
          object@expectations$Y[[m]][[g]] <- object@expectations$Y[[m]][[g]][,samples_g,drop=FALSE]
        }  
      }
    }

    # Subset data
    if (length(object@data)>0) { 
      for (m in views_names(object)) {
        object@data[[m]][[g]] <- object@data[[m]][[g]][,samples_g,drop=FALSE]
      }
    }
    
    # Subset imputed data
    if (length(object@imputed_data)>0) { 
      for (m in views_names(object)) {
        object@imputed_data[[m]][[g]] <- object@imputed_data[[m]][[g]][,samples_g,drop=FALSE]
      }
    }
    
  }

  # Subset sample metadata
  object@samples_metadata <- object@samples_metadata[object@samples_metadata$sample %in% samples,]
  
  # Update dimensionality
  object@dimensions[["N"]] <- sapply(tmp, length)
  
  # Update sample names
  samples_names(object) <- tmp
  
  # Sanity checks
  stopifnot(object@samples_metadata$sample == unlist(lapply(object@data[[1]],colnames)))
  stopifnot(object@samples_metadata$sample == unlist(lapply(object@expectations$Z,rownames)))
  stopifnot(object@samples_metadata$sample == unlist(lapply(object@expectations$Y,colnames)))
  
  # Remove variance explained estimates  
  warning("After subsetting the samples the variance explained estimates are not valid anymore, removing them...")
  object@cache[["variance_explained"]] <- NULL
  
  return(object)
}


#' @title Subset features
#' @name subset_features
#' @description Method to subset (or sort) features
#' @param object a \code{\link{MOFA}} object.
#' @param view character vector with the view name or integer with the view index
#' @param features character vector with the sample names, numeric vector with the feature indices 
#' or logical vector with the samples to be kept as TRUE.
#' @return A \code{\link{MOFA}} object
#' @export

subset_features <- function(object, view, features) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(length(features) <= sapply(object@dimensions[["D"]], sum))
  warning("Removing features a posteriori is fine for an exploratory analysis, but we recommend removing them before training!")

  if (is.numeric(view)) view <- views_names(object)[view]
  stopifnot(all(view %in% views_names(object)))  

  # Define features
  if (is.character(features)) {
    stopifnot(all(features %in% features_names(object)[[view]]))
  } else {
    features <- features_names(object)[[view]][features]
  }
  
  # Subset relevant slots
  if (length(object@expectations)>0) {
    if ("W" %in% names(object@expectations) & length(object@expectations$W)>0)
      object@expectations$W <- lapply(object@expectations$W, function(x) x[features,, drop=FALSE])
    if ("Y" %in% names(object@expectations) & length(object@expectations$Y)>0)
      object@expectations$Y[[view]] <- lapply(object@expectations$Y[[view]], function(x) x[features,])

  if (length(object@data)>0)
    object@data <- lapply(object@data, function(x) sapply(x, function(y) y[features,], simplify = FALSE, USE.NAMES = TRUE))

  if (length(object@expectations)>0)
  object@intercepts <- lapply(object@intercepts, function(x) sapply(x, function(y) y[features], simplify = FALSE, USE.NAMES = TRUE))

  if (length(object@imputed_data) != 0) {
    stop()
    # object@imputed_data <- lapply(object@imputed_data, function(x) sapply(x, function(y) y[,samples], simplify = FALSE, USE.NAMES = TRUE)) 
    }
  }
  # Update dimensionality
  object@dimensions[["D"]][[view]] <- length(features)
  
  # Update features names
  features_names(object)[[view]] <- features
  
  # Remove variance explained estimates  
  warning("After subsetting the features the variance explained estimates are not valid anymore, removing them...")
  object@cache[["variance_explained"]] <- NULL
  
  return(object)
}
