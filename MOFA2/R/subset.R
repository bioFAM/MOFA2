
################################
## Functions to do subsetting ##
################################

#' @title Subset groups
#' @name subset_groups
#' @description Method to subset (or sort) groups
#' @param object a \code{\link{MOFA}} object.
#' @param groups character vector with the groups names, numeric vector with the groups indices
#' or logical vector with the groups to be kept as TRUE.
#' @export
subset_groups <- function(object, groups) {
  
  # Sanity checks
  if (class(object) != "MOFA") stop("'object' has to be an instance of MOFA")
  stopifnot(length(groups) <= object@dimensions[["G"]])
  
  if (is.numeric(groups) | is.logical(groups))  {
    groups <- groups_names(object)[groups] 
  } else {
    stopifnot(all(groups %in% groups_names(object)))
  }
  
  # Subset relevant slots
  if ("Z" %in% names(object@expectations) & length(object@expectations$Z)>0)
    object@expectations$Z <- object@expectations$Z[groups]
  if ("Y" %in% names(object@expectations) & length(object@expectations$Y)>0)
    object@expectations$Y <- sapply(object@expectations$Y, function(x) x[groups], simplify = F, USE.NAMES = T) 
  
  object@data <- sapply(object@data, function(x) x[groups], simplify = F, USE.NAMES = T) 
  
  # Update dimensionality
  object@dimensions[["G"]] <- length(groups)
  object@dimensions[["N"]] <- object@dimensions[["N"]][groups]
  
  # Subset sample metadata
  stopifnot(groups%in%unique(object@samples_metadata$group_name))
  object@samples_metadata <- object@samples_metadata[object@samples_metadata$group_name %in% groups,]
  object@samples_metadata$group_name <- factor(object@samples_metadata$group_name, levels=groups)
  
  # Re-order samples
  samples <- unname(unlist(lapply(object@data[[1]],colnames)))
  object@samples_metadata <- object@samples_metadata[match(samples, object@samples_metadata$sample_name),]
  
  # Sanity checks
  stopifnot(object@samples_metadata$sample_name == unlist(lapply(object@data[[1]],colnames)))
  stopifnot(object@samples_metadata$sample_name == unlist(lapply(object@expectations$Z,rownames)))
  stopifnot(object@samples_metadata$sample_name == unlist(lapply(object@expectations$Y,colnames)))
  
  # Update groups names
  # groups_names(object) <- groups
  object@data_options$groups <- groups
  
  # Re-compute variance explained
  object@cache[["variance_explained"]] <- calculate_variance_explained(object)
  
  return(object)
}


#' @title Subset views
#' @name subset_views
#' @description Method to subset (or sort) views
#' @param object a \code{\link{MOFA}} object.
#' @param views character vector with the views names, numeric vector with the views indices,
#' or logical vector with the views to be kept as TRUE.
#' @export
subset_views <- function(object, views) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(length(views) <= object@dimensions[["M"]])
  # warning("Removing views a posteriori is fine for an exploratory analysis, but you should removing them before training!")
  
  if (is.numeric(views) | is.logical(views))  {
    views <- views_names(object)[views] 
  } else {
    stopifnot(all(views %in% views_names(object)))
  }
  
  # Subset relevant slots
  object@expectations$W <- object@expectations$W[views]
  object@expectations$Y <- object@expectations$Y[views]
  object@expectations$Tau <- object@expectations$Tau[views]
  object@data <- object@data[views]
  # object@data <- object@data[,,,]
  # if (length(object@imputed_data) != 0) {
  #   object@imputed_data <- lapply(object@imputed_data, function(x) sapply(x, function(y) y[,samples], simplify = F, USE.NAMES = T)) 
  # }
  
  # Update dimensionality
  object@dimensions[["M"]] <- length(views)
  object@dimensions[["D"]] <- object@dimensions[["D"]][views]
  
  # Update view names
  views_names(object) <- views
  
  # Re-compute variance explained
  object@cache[["variance_explained"]] <- calculate_variance_explained(object)
  
  return(object)
}


#' @title Subset factors
#' @name subset_factors
#' @description Method to subset (or sort) factors
#' @param object a \code{\link{MOFA}} object.
#' @param factors character vector with the factor names, or numeric vector with the index of the factors, 
#' or logical vector with the views to be kept as TRUE.
#' @export
subset_factors <- function(object, factors) {
  
  # Sanity checks
  if (class(object) != "MOFA") stop("'object' has to be an instance of MOFA")
  stopifnot(length(factors) <= object@dimensions[["K"]])
  
  # Get factors
  if (is.numeric(factors) | is.logical(factors)) {
    stopifnot(length(factors) == length(unique(factors)))
    factors <- factors_names(object)[factors]
  } else if (is.logical(factors)) {
    factors <- factors_names(object)[factors]
  } else if (is.character(factors)) { 
    stopifnot(length(factors) == length(unique(factors)))
    stopifnot(all(factors %in% factors_names(object)))
  }
  
  # Subset expectations
  nodes_with_factors <- list(nodes = c("Z", "W", "AlphaZ", "AlphaW", "ThetaZ", "ThetaW"), axes = c(2, 2, 0, 0, 0, 0))
  stopifnot(all(nodes_with_factors$axes %in% c(0, 1, 2)))

  for (i in seq_len(length(nodes_with_factors$nodes))) {
    node <- nodes_with_factors$nodes[i]
    axis <- nodes_with_factors$axes[i]
    if (node %in% names(object@expectations)) {
      if (axis == 1) {
        object@expectations[[node]] <- sapply(object@expectations[[node]], function(x) x[factors,], simplify = F, USE.NAMES = T)
      } else if (axis == 2) {
        object@expectations[[node]] <- sapply(object@expectations[[node]], function(x) x[,factors], simplify = F, USE.NAMES = T)
      } else {
        object@expectations[[node]] <- sapply(object@expectations[[node]], function(x) x[factors], simplify = F, USE.NAMES = T)
      }
    }
  }
  
  # Update dimensionality
  object@dimensions[["K"]] <- length(factors)
  
  # Re-compute variance explained
  object@cache[["variance_explained"]] <- calculate_variance_explained(object)
  
  # Update factor names
  factors_names(object) <- paste0("Factor", as.character(seq_len(object@dimensions[["K"]])))
  
  
  return(object)
}



#' @title Subset samples
#' @name subset_samples
#' @description Method to subset (or sort) samples
#' @param object a \code{\link{MOFA}} object.
#' @param samples character vector with the sample names, numeric vector with the sample indices or logical vector with the samples to be kept as TRUE.
#' @export
subset_samples <- function(object, samples) {
  
  # Sanity checks
  if (class(object) != "MOFA") stop("'object' has to be an instance of MOFA")
  stopifnot(length(samples) <= sum(object@dimensions[["N"]]))
  stopifnot(all(samples %in% unlist(samples_names(object))))
  stopifnot(!any(duplicated(samples)))
  
  # Subset sample metadata
  # object@samples_metadata <- object@samples_metadata[match(samples, object@samples_metadata$sample_name),]
  # groups <- as.character(unique(object@samples_metadata$group_name))
  # object@samples_metadata$group_name <- factor(object@samples_metadata$group_name, levels=groups)
  
  # Check if an entire group needs to be removed
  groups <- as.character(unique(object@samples_metadata[match(samples, object@samples_metadata$sample_name),]$group_name))
  if (length(groups)<length(groups_names(object))) object <- subset_groups(object, groups)
  
  
  # Subset data and expectations
  
  tmp <- lapply(groups, function(g) samples_names(object)[[g]][samples_names(object)[[g]] %in% samples])
  names(tmp) <- groups
  
  for (g in groups) {
    samples_g <- tmp[[g]]
    
    if ("Z" %in% names(object@expectations) & length(object@expectations$Z)>0) {
      object@expectations$Z[[g]] <- object@expectations$Z[[g]][samples_g,, drop=F]
    }
    
    if ("Y" %in% names(object@expectations) & length(object@expectations$Y)>0) {
      for (m in views_names(object)) {
        object@expectations$Y[[m]][[g]] <- object@expectations$Y[[m]][[g]][,samples_g,drop=F]
      }  
    }

    if (length(object@data)>0) { 
      for (m in views_names(object)) {
        object@data[[m]][[g]] <- object@data[[m]][[g]][,samples_g,drop=F]
      }
    }
    
    if (length(object@imputed_data)>0) { 
      for (m in views_names(object)) {
        object@imputed_data[[m]][[g]] <- object@imputed_data[[m]][[g]][,samples_g,drop=F]
      }
    }
    
  }

  # Update dimensionality
  object@dimensions[["N"]] <- sapply(tmp, length)
  
  # Update sample names
  samples_names(object) <- tmp
  
  # Sanity checks
  stopifnot(object@samples_metadata$sample_name == unlist(lapply(object@data[[1]],colnames)))
  stopifnot(object@samples_metadata$sample_name == unlist(lapply(object@expectations$Z,rownames)))
  stopifnot(object@samples_metadata$sample_name == unlist(lapply(object@expectations$Y,colnames)))
  
  # Re-compute variance explained
  object@cache[["variance_explained"]] <- calculate_variance_explained(object)
  
  return(object)
}


#' @title Subset features
#' @name subset_features
#' @description Method to subset (or sort) features
#' @param object a \code{\link{MOFA}} object.
#' @param view character vector with the view name or integer with the view index
#' @param features character vector with the sample names, numeric vector with the feature indices 
#' or logical vector with the samples to be kept as TRUE.
#' @export
subset_features <- function(object, view, features) {
  
  # Sanity checks
  if (class(object) != "MOFA") stop("'object' has to be an instance of MOFA")
  stopifnot(length(features) <= sapply(object@dimensions[["D"]], sum))
  warning("Removing features a posteriori is fine for an exploratory analysis, but we recommend removing them before training!")

  if (is.numeric(view)) view <- views_names(object)[view]
  stopifnot(all(view %in% views_names(object)))  

  
  # Get samples
  if (is.character(features)) {
    stopifnot(all(features %in% features_names(object)[[view]]))
  } else {
    features <- features_names(object)[[view]][features]
  }
  
  # Subset relevant slots
  object@expectations$W <- lapply(object@expectations$W, function(x) x[features,, drop=F])
  object@expectations$Y[[view]] <- lapply(object@expectations$Y[[view]], function(x) x[features,])
  object@expectations$Tau[[view]] <- lapply(object@expectations$Tau[[view]], function(x) x[features,])
  object@data <- lapply(object@data, function(x) sapply(x, function(y) y[features,], simplify = F, USE.NAMES = T))
  object@intercepts <- lapply(object@intercepts, function(x) sapply(x, function(y) y[features], simplify = F, USE.NAMES = T))
  if (length(object@imputed_data) != 0) {
    object@imputed_data <- lapply(object@imputed_data, function(x) sapply(x, function(y) y[,samples], simplify = F, USE.NAMES = T)) 
  }

  # Update dimensionality
  object@dimensions[["D"]][[view]] <- length(features)
  
  # Update features names
  features_names(object)[[view]] <- features
  
  # Re-compute variance explained
  object@cache[["variance_explained"]] <- calculate_variance_explained(object)
  
  return(object)
}
