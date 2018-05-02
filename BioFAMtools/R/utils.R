
.inferLikelihoods <- function(object) {
  likelihood <- rep(x="gaussian", times=object@Dimensions$M)
  names(likelihood) <- viewNames(object)
  
  for (view in viewNames(object)) {
    data <- getTrainData(object, view)[[1]]
    # if (all(data %in% c(0,1,NA))) {
    if (length(unique(data[!is.na(data)]))==2) {
      likelihood[view] <- "bernoulli"
    } else if (all(data[!is.na(data)]%%1==0)) {
      likelihood[view] <- "poisson"
    }
  }
  
  return(likelihood)
}

.updateOldModel <- function(object) {
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")  
  
  # Update node names
  if ("SW" %in% names(object@Expectations)) {
    # object@ModelOptions$schedule[object@ModelOptions$schedule == "SW"] <- "W" # schedule is depreciated from ModelOptions
    names(object@Expectations)[names(object@Expectations) == "SW"] <- "W"
    colnames(object@TrainStats$elbo_terms)[colnames(object@TrainStats$elbo_terms)=="SW"] <- "W"
  }
  if ("SZ" %in% names(object@Expectations)) {
    names(object@Expectations)[names(object@Expectations) == "SZ"] <- "Z"
    colnames(object@TrainStats$elbo_terms)[colnames(object@TrainStats$elbo_terms)=="SZ"] <- "Z"
  }

  # Update expectations
  if (is.list(object@Expectations$Z[[1]]) & ("E" %in% names(object@Expectations$Z[[1]]))) {
    for (m in viewNames(object)) {
      nodes <- c("W", "Tau", "AlphaW", "SigmaAlphaW", "ThetaW")
      for (node in nodes){
        if (node %in% names(object@Expectations)){
          object@Expectations[[node]][[m]] <- object@Expectations[[node]][[m]]$E
        }
      }
    }
    for (h in groupNames(object)) {
      nodes <- c("Z", "AlphaZ", "SigmaZ", "ThetaZ")
      for (node in nodes){
        if (node %in% names(object@Expectations)){
          object@Expectations[[node]][[h]] <- object@Expectations[[node]][[h]]$E
        }
      }
    }
    for (m in viewNames(object)) {
      for (h in groupNames(object)) {
        object@Expectations$Y[[m]][[h]] <- object@Expectations$Y[[m]][[h]]$E
      }
    }
  }
  
  
  # update learnMean to learnIntercept
  if ("learnMean" %in% names(object@ModelOptions)) {
    tmp <- names(object@ModelOptions)
    tmp[tmp=="learnMean"] <- "learnIntercept"
    names(object@ModelOptions) <- tmp
  }
  object@ModelOptions$learnIntercept <- as.logical(object@ModelOptions$learnIntercept)
  
  
  return(object)
}

# Set view names and group names for nested list objects (e.g. Y)
.name_views_and_groups <- function(nested_list, view_names, group_names) {
  names(nested_list) <- view_names
  for (view in view_names) { names(nested_list[[view]]) <- group_names }
  nested_list
}

# Function to find factors that act like an intercept term for the sample, 
# which means that they capture global mean effects
findInterceptFactors <- function(object, cor_threshold = 0.8) {
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")  
  
  data <- getTrainData(object)
  factors <- getFactors(object, include_intercept = F)
  
  r <- lapply(data, function(x) abs(cor(apply(x,2,mean),factors, use="complete.obs")))
  for (i in names(r)) {
    if (any(r[[i]]>cor_threshold))
      cat(paste0("Warning: factor ",which(r[[i]]>cor_threshold)," is capturing a size factor effect in ", i, " view, which indicates that input data might not be properly normalised...\n"))
  }
}


subset_augment <- function(mat, pats) {
  pats <- unique(pats)
  mat <- t(mat)
  aug_mat <- matrix(NA, ncol=ncol(mat), nrow=length(pats))
  aug_mat <- mat[match(pats,rownames(mat)),,drop=FALSE]
  rownames(aug_mat) <- pats
  colnames(aug_mat) <- colnames(mat)
  return(t(aug_mat))
}


detectPassengers <- function(object, views = "all", groups = "all", factors = "all", r2_threshold = 0.03) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  
  # Define views
  if (paste0(views, sep="", collapse="") == "all") { 
    views <- viewNames(object) 
  } else {
    stopifnot(all(views %in% viewNames(object)))  
  }
  M <- length(views)

  # Define groups
  if (paste0(groups, sep="", collapse="") == "all") { 
    groups <- groupNames(object) 
  } else {
    stopifnot(all(groups %in% groupNames(object)))  
  }
  H <- length(groups)
  
  # Define factors
  factors <- as.character(factors)
  if (paste0(factors, collapse="") == "all") { 
    factors <- factorNames(object)
  } else {
    stopifnot(all(factors %in% factorNames(object)))  
  }
  
  # Collect relevant data
  Z <- getFactors(object)
  
  # Identify factors unique to a single view by calculating relative R2 per factor
  r2 <- calculateVarianceExplained(object, views = views, groups = groups, factors = factors)$R2PerFactor
  unique_factors <- unique(unlist(lapply(groups, function(h) names(which(rowSums(r2[[h]]>=r2_threshold)==1)) )))
  
  # Mask samples that are unique in the unique factors
  missing <- lapply(getTrainData(object, views, groups), function(views) {
    lapply(views, function(group) {
      sampleNames(object)[apply(group, 2, function(x) all(is.na(x)))]
    })
  })
  missing <- .name_views_and_groups(missing, viewNames(object), groupNames(object))
  for (factor in unique_factors) {
    # view <- names(which(r2[factor,]>=r2_threshold))
    for (h in groups) {
      view <- colnames(r2[[h]][,which(r2[[h]][factor,]>=r2_threshold),drop=F])
      missing_samples <- missing[[view]][[h]]
      if (length(missing_samples) > 0) {
        Z[[h]][missing_samples,factor] <- NA
      }
    }
  }
  
  # Replace the latent matrix
  object@Expectations$Z <- Z
  
  return(object)
  
}


flip_factor <- function(model, factor){
  for(groupnm in names(model@Expectations$Z)) {
    model@Expectations$Z[[groupnm]][,factor] <- - model@Expectations$Z[[groupnm]][,factor]
  }
  for(viewnm in names(model@Expectations$W)) {
    model@Expectations$W[[viewnm]][,factor] <- -model@Expectations$W[[viewnm]][,factor]
  }
return(model)
}




.check_and_get_views <- function(object, views) {
  if (is.numeric(views)) {
    stopifnot(all(views <= object@Dimensions$M))
    viewNames(object)[views] 
  } else {
    if (paste0(views, sep = "", collapse = "") == "all") { 
      viewNames(object)
    } else {
      stopifnot(all(views %in% viewNames(object)))
      views
    }
  }
}


.check_and_get_groups <- function(object, groups) {
  if (is.numeric(groups)) {
    stopifnot(all(groups <= object@Dimensions$M))
    groupNames(object)[groups] 
  } else {
    if (paste0(groups, sep = "", collapse = "") == "all") { 
      groupNames(object)
    } else {
      stopifnot(all(groups %in% groupNames(object)))
      groups
    }
  }
}