
.infer_likelihoods <- function(object) {
  
  # Gaussian by default
  likelihood <- rep(x="gaussian", times=object@dimensions$M)
  names(likelihood) <- views(object)
  
  for (m in views(object)) {
    # data <- get_data(object, views=m)[[1]][[1]]  # take only first group
    data <- object@data[[m]][[1]]
    
    # bernoulli
    if (length(unique(data[!is.na(data)]))==2) {
      likelihood[m] <- "bernoulli"
    # poisson
    } else if (all(data[!is.na(data)]%%1==0)) {
      likelihood[m] <- "poisson"
    }
  }
  
  return(likelihood)
}

# Set view names and group names for nested list objects (e.g. Y)
.name_views_and_groups <- function(nested_list, view_names, group_names) {
  names(nested_list) <- view_names
  for (view in view_names) { names(nested_list[[view]]) <- group_names }
  nested_list
}

#' @importFrom stats sd
.detect_outliers <- function(object, groups = "all", factors = "all") {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Define groups
  groups <- .check_and_get_groups(object, groups)
  H <- length(groups)
  
  # Define factors
  factors <- .check_and_get_factors(object, factors)
  
  for (k in factors) {
    for (g in groups) {
    
      Z <- get_factors(object, groups=g, factors=k)[[1]][,1]
      Z <- Z[!is.na(Z)]
      
      cutoff <- 2 * 1.96
      tmp <- abs(Z - mean(Z)) / sd(Z)

      outliers <- names(which(tmp>cutoff & abs(Z)>0.5))
      
      if (length(outliers)>0) {
        object@expectations$Z[[g]][,k][outliers] <- NA
      }
      
    }
  }
  
  # re-compute variance explained
  object@cache[["variance_explained"]] <- calculate_variance_explained(object)
  
  return(object)
}


.flip_factor <- function(model, factor){
  for(g in names(model@expectations$Z)) {
    model@expectations$Z[[g]][,factor] <- - model@expectations$Z[[g]][,factor]
  }
  for(m in names(model@expectations$W)) {
    model@expectations$W[[m]][,factor] <- -model@expectations$W[[m]][,factor]
  }
return(model)
}


.check_and_get_factors <- function(object, factors) {
  stopifnot(!any(duplicated(factors)))
  if (is.numeric(factors)) {
    stopifnot(all(factors <= object@dimensions$K))
    factors(object)[factors] 
  } else {
    if (paste0(factors, collapse = "") == "all") { 
      factors(object)
    } else {
      stopifnot(all(factors %in% factors(object)))
      factors
    }
  }
  
}

.check_and_get_views <- function(object, views) {
  stopifnot(!any(duplicated(views)))
  if (is.numeric(views)) {
    stopifnot(all(views <= object@dimensions$M))
    views(object)[views] 
  } else {
    if (paste0(views, sep = "", collapse = "") == "all") { 
      views(object)
    } else {
      stopifnot(all(views %in% views(object)))
      views
    }
  }
}


.check_and_get_groups <- function(object, groups) {
  stopifnot(!any(duplicated(groups)))
  if (is.numeric(groups)) {
    stopifnot(all(groups <= object@dimensions$G))
    groups(object)[groups] 
  } else {
    if (paste0(groups, collapse = "") == "all") { 
      groups(object)
    } else {
      stopifnot(all(groups %in% groups(object)))
      groups
    }
  }
}


.check_and_get_samples <- function(object, samples) {
  stopifnot(!any(duplicated(samples)))
  if (is.numeric(samples)) {
    stopifnot(all(samples <= sum(object@dimensions$N)))
    unlist(samples(object))[samples] 
  } else {
    if (paste0(samples, collapse = "") == "all") { 
      unlist(samples(object))
    } else {
      stopifnot(all(samples %in% unlist(samples(object))))
      samples
    }
  }
}

.check_and_get_features <- function(object, features) {
  stopifnot(!any(duplicated(features)))
  if (is.numeric(features)) {
    stopifnot(all(features <= sum(object@dimensions$D)))
    unlist(features(object))[features] 
  } else {
    if (paste0(features, collapse = "") == "all") { 
      unlist(features(object))
    } else {
      stopifnot(all(features %in% unlist(features(object))))
      features
    }
  }
}

.get_nodes_types <- function() {
  nodes_types <- list(
    multiview_nodes  = c("W", "AlphaW", "ThetaW"),
    multigroup_nodes = c("Z", "AlphaZ", "ThetaZ"),
    twodim_nodes     = c("Y", "Tau")
  )
}

setClass("matrix_placeholder", 
         slots=c(rownames = "ANY",
                 colnames = "ANY",
                 nrow     = "integer",
                 ncol     = "integer")
)

setMethod("rownames", "matrix_placeholder", function(x) { x@rownames })
setMethod("colnames", "matrix_placeholder", function(x) { x@colnames })
setMethod("nrow", "matrix_placeholder", function(x) { x@nrow })
setMethod("ncol", "matrix_placeholder", function(x) { x@ncol })

setReplaceMethod("rownames", signature(x = "matrix_placeholder"),
  function(x, value) { x@rownames <- value; x@nrow <- length(value); x })
setReplaceMethod("colnames", signature(x = "matrix_placeholder"),
  function(x, value) { x@colnames <- value; x@ncol <- length(value); x })

.create_matrix_placeholder <- function(rownames, colnames) {
  mx <- new("matrix_placeholder")
  mx@rownames <- rownames
  mx@colnames <- colnames
  mx@nrow <- length(rownames)
  mx@ncol <- length(colnames)
  mx
}
