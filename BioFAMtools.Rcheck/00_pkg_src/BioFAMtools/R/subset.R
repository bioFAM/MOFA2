
################################
## Functions to do subsetting ##
################################

#' @title Subset groups
#' @name subset_groups
#' @description Method to subset (or sort) groups
#' @param object a \code{\link{BioFAModel}} object.
#' @param groups character vector with the groups names, numeric vector with the groups indices or logical vector with the groups to be kept as TRUE.
#' @export
subset_groups <- function(object, groups) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  stopifnot(length(groups) <= object@dimensions[["P"]])
  
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
  if ("Tau" %in% names(object@expectations) & length(object@expectations$Tau)>0)
    object@expectations$Tau <- sapply(object@expectations$Tau, function(x) x[groups], simplify = F, USE.NAMES = T)
  
  object@training_data <- sapply(object@training_data, function(x) x[groups], simplify = F, USE.NAMES = T) 
  
  # Update dimensionality
  object@dimensions[["P"]] <- length(groups)
  object@dimensions[["N"]] <- object@dimensions[["N"]][groups]
  
  # Update view names
  groups_names(object) <- groups
  
  # Re-compute variance explained
  object@cache[["variance_explained"]] <- calculate_variance_explained(object)
  
  return(object)
}


#' @title Subset views
#' @name subset_views
#' @description Method to subset (or sort) views
#' @param object a \code{\link{BioFAModel}} object.
#' @param views character vector with the views names, numeric vector with the views indices or logical vector with the views to be kept as TRUE.
#' @export
subset_views <- function(object, views) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
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
  object@training_data <- object@training_data[views]
  # object@input_data <- object@input_data[,,,]
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
#' @param object a \code{\link{BioFAModel}} object.
#' @param factors character vector with the factor names, or numeric vector with the index of the factors.
#' @export
subset_factors <- function(object, factors) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
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

  for (i in 1:length(nodes_with_factors$nodes)) {
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
  
  # Update factor names
  factors_names(object) <- paste0("Factor", as.character(1:object@dimensions[["K"]]))
  
  # Re-compute variance explained
  object@cache[["variance_explained"]] <- calculate_variance_explained(object)
  
  return(object)
}



#' @title Subset samples
#' @name subset_samples
#' @description Method to subset (or sort) samples
#' @param object a \code{\link{BioFAModel}} object.
#' @param samples character vector with the sample names, numeric vector with the sample indices or logical vector with the samples to be kept as TRUE.
#' @export
subset_samples <- function(object, samples) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  stopifnot(length(samples) <= sum(object@dimensions[["N"]]))
  stopifnot(all(samples %in% unlist(samples_names(object))))
  
  groups <- object@data_options$samples_groups
  new_samples_names <-  sapply(groups, function(g) samples_names(object)[[g]][samples_names(object)[[g]] %in% samples])
  
  # Subset relevant slots
  for (g in names(new_samples_names)) {
    samples_g <- new_samples_names[[g]]
    
    object@expectations$Z[[g]] <- object@expectations$Z[[g]][samples_g,, drop=F]
    for (m in views_names(object)) {
      object@expectations$Y[[m]][[g]] <- object@expectations$Y[[m]][[g]][,samples_g,drop=F]
      object@training_data[[m]][[g]] <- object@training_data[[m]][[g]][,samples_g,drop=F]
      
    }
  }
  # if (length(object@ImputedData)==0) { object@ImputedData <- sapply(object@ImputedData, function(x) sapply(x, function(y) y[,samples], simplify = F, USE.NAMES = T)) }

  # Update dimensionality
  object@dimensions[["N"]] <- sapply(new_samples_names, length)
  
  # Update sample names
  samples_names(object) <- new_samples_names
  
  # Re-compute variance explained
  object@cache[["variance_explained"]] <- calculate_variance_explained(object)
  
  return(object)
}


#' @title Subset features
#' @name subset_features
#' @description Method to subset (or sort) features
#' @param object a \code{\link{BioFAModel}} object.
#' @param view view name or view index
#' @param features character vector with the sample names, numeric vector with the feature indices or logical vector with the samples to be kept as TRUE.
#' @export
subset_features <- function(object, view, features) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
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
  object@training_data <- lapply(object@training_data, function(x) sapply(x, function(y) y[features,], simplify = F, USE.NAMES = T))
  # object@input_data <- object@input_data[,,,] 
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
