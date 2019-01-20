
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

detect_passengers <- function(object, views = "all", groups = "all", factors = "all") {
  
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
  
  # Calculate variance explained
  r2 <- calculate_variance_explained(object, views = views, groups = groups, factors = factors)$r2_per_factor
  
  # (...)  
  data <- get_training_data(object, views, groups)
  missing <- lapply(groups, function(g) { 
    tmp <- lapply(views, function(m) {
      names(which(apply(data[[m]][[g]],2,function(x) all(is.na(x)))))
    }); names(tmp) <- views
    tmp
  }); names(missing) <- groups
  
  for (k in factors) {
    for (g in groups) {
      for (m in views) {
        samples <- samples_names(object)[[g]]
        if (any(samples%in%missing[[g]][[m]])) {
          foo <- predict(object, groups=g, factors=k, views=m)[[1]][[1]]
          
          # select top features 
          genes <- names(tail(sort(abs(get_weights(object, views=m, factors=k)[[1]][,1])), n=100))
          
          # Define foreground set
          bar <- apply(abs(foo[genes,]),2,mean,na.rm=T)
          bar <- bar[samples %in%missing[[g]][[m]]]
          
          # Define background set: group with the highest activity
          background_g <- names(which.max(sapply(r2, function(x) x[k,][[m]])))
          foo.background <- predict(object, groups=background_g, factors=k, views=m)[[1]][[1]]
          bar.background <- apply(abs(foo.background[genes,]),2,mean,na.rm=T)
          
          # detect outliers by comparing statistics of foreground versus background
          max.value <- max(bar.background,na.rm=T)
          outliers <- names(which(bar>max.value ))
          if (length(outliers)>0) {
            object@expectations$Z[[g]][,k][outliers] <- NA
          }
          
        }
      }
    }
  }
  
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
