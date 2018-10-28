
.infer_likelihoods <- function(object) {
  
  # Gaussian by default
  likelihood <- rep(x="gaussian", times=object@dimensions$M)
  names(likelihood) <- views_names(object)
  
  for (view in views_names(object)) {
    data <- get_training_data(object, view)[[1]][[1]]  # take only first group
    
    # bernoulli
    if (length(unique(data[!is.na(data)]))==2) {
      likelihood[view] <- "bernoulli"
    # poisson
    } else if (all(data[!is.na(data)]%%1==0)) {
      likelihood[view] <- "poisson"
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

detect_passengers <- function(object, views = "all", groups = "all", factors = "all", r2_threshold = 0.03) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  
  # Define views
  if (paste0(views, sep="", collapse="") == "all") { 
    views <- views_names(object) 
  } else {
    stopifnot(all(views %in% views_names(object)))  
  }
  M <- length(views)

  # Define groups
  if (paste0(groups, sep="", collapse="") == "all") { 
    groups <- groups_names(object) 
  } else {
    stopifnot(all(groups %in% groups_names(object)))  
  }
  H <- length(groups)
  
  # Define factors
  factors <- as.character(factors)
  if (paste0(factors, collapse="") == "all") { 
    factors <- factors_names(object)
  } else {
    stopifnot(all(factors %in% factors_names(object)))  
  }
  
  # Collect relevant data
  Z <- get_factors(object)
  
  # Identify factors unique to a single view by calculating relative R2 per factor
  r2 <- calculate_variance_explained(object, views = views, groups = groups, factors = factors)$r2_per_factor
  unique_factors <- unique(unlist(lapply(groups, function(p) names(which(rowSums(r2[[p]]>=r2_threshold)==1)) )))
  
  # Mask samples that are unique in the unique factors
  missing <- lapply(get_training_data(object, views, groups), function(views) {
    lapply(views, function(group) {
      names(which(apply(group,2,function(x) all(is.na(x)))))
    })
  })

  missing <- .name_views_and_groups(missing, views_names(object), groups_names(object))
  for (fctr in unique_factors) {
    # view <- names(which(r2[fctr,]>=r2_threshold))
    for (p in groups) {
      view <- colnames(r2[[p]][,which(r2[[p]][fctr,]>=r2_threshold), drop=F])
      if (!is.null(view)) {
        missing_samples <- missing[[view]][[p]]
        if (length(missing_samples) > 0) {
          Z[[p]][missing_samples, fctr] <- rep(NA,length(missing_samples))
        }
      }
    }
  }
  
  # Replace the latent matrix
  object@expectations$Z <- Z
  
  return(object)
  
}


.flip_factor <- function(model, factor){
  for(groupnm in names(model@expectations$Z)) {
    model@expectations$Z[[groupnm]][,factor] <- - model@expectations$Z[[groupnm]][,factor]
  }
  for(viewnm in names(model@expectations$W)) {
    model@expectations$W[[viewnm]][,factor] <- -model@expectations$W[[viewnm]][,factor]
  }
return(model)
}




.check_and_get_views <- function(object, views) {
  if (is.numeric(views)) {
    stopifnot(all(views <= object@dimensions$M))
    views_names(object)[views] 
  } else {
    if (paste0(views, sep = "", collapse = "") == "all") { 
      views_names(object)
    } else {
      stopifnot(all(views %in% views_names(object)))
      views
    }
  }
}


.check_and_get_groups <- function(object, groups) {
  if (is.numeric(groups)) {
    stopifnot(all(groups <= object@dimensions$P))
    groups_names(object)[groups] 
  } else {
    if (paste0(groups, sep = "", collapse = "") == "all") { 
      groups_names(object)
    } else {
      stopifnot(all(groups %in% groups_names(object)))
      groups
    }
  }
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