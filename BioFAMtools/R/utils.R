
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
    # object@ModelOpts$schedule[object@ModelOpts$schedule == "SW"] <- "W" # schedule is depreciated from ModelOpts
    names(object@Expectations)[names(object@Expectations) == "SW"] <- "W"
    colnames(object@TrainStats$elbo_terms)[colnames(object@TrainStats$elbo_terms)=="SW"] <- "W"
  }
  if ("SZ" %in% names(object@Expectations)) {
    names(object@Expectations)[names(object@Expectations) == "SZ"] <- "Z"
    colnames(object@TrainStats$elbo_terms)[colnames(object@TrainStats$elbo_terms)=="SZ"] <- "Z"
  }
  
  # Update expectations
  if (is.list(object@Expectations$Z)) {
    object@Expectations$Z <- object@Expectations$Z$E
    
    nodes=c("AlphaZ","SigmaZ","ThetaZ")
    for (node in nodes){
      if (node %in% names(object@Expectations)){
        object@Expectations$node <- object@Expectations$node$E
      }
    }
    
    for (view in viewNames(object)) {
      
      object@Expectations$W[[view]] <- object@Expectations$W[[view]]$E
      object@Expectations$Y[[view]] <- object@Expectations$Y[[view]]$E
      object@Expectations$Tau[[view]] <- object@Expectations$Tau[[view]]$E
      
      nodes=c("AlphaW","SigmaAlphaW","ThetaW")
      for (node in nodes){
        if (node %in% names(object@Expectations)){
          object@Expectations$node[[view]] <- object@Expectations$node[[view]]$E
        }
      }
    }
    
  }
  
  # update learnMean to learnIntercept
  if ("learnMean" %in% names(object@ModelOpts)) {
    tmp <- names(object@ModelOpts)
    tmp[tmp=="learnMean"] <- "learnIntercept"
    names(object@ModelOpts) <- tmp
  }
  object@ModelOpts$learnIntercept <- as.logical(object@ModelOpts$learnIntercept)
  
  
  return(object)
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


detectPassengers <- function(object, views = "all", factors = "all", r2_threshold = 0.03, batches = "all") {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  
  # Define views
  if (paste0(views, sep="", collapse="") == "all") { 
    views <- viewNames(object) 
  } else {
    stopifnot(all(views %in% viewNames(object)))  
  }
  M <- length(views)

  # Define batches
  if (paste0(batches, sep="", collapse="") == "all") { 
    batches <- batchNames(object) 
  } else {
    stopifnot(all(batches %in% batchNames(object)))  
  }
  H <- length(batches)
  
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
  r2 <- calculateVarianceExplained(object, views = views, factors = factors, batches = batches)$R2PerFactor
  unique_factors <- names(which(rowSums(r2>=r2_threshold)==1))
  
  # Mask samples that are unique in the unique factors
  missing <- sapply(getTrainData(object,views), function(view) sampleNames(object)[apply(view, 2, function(x) all(is.na(x)))] )
  names(missing) <- viewNames(object)
  for (factor in unique_factors) {
    # view <- names(which(r2[factor,]>=r2_threshold))
    view <- colnames(r2[,which(r2[factor,]>=r2_threshold),drop=F])
    missing_samples <- missing[[view]]
    if (length(missing_samples)>0) {
      Z[missing_samples,factor] <- NA
    }
  }
  
  # Replace the latent matrix
  object@Expectations$Z <- Z
  
  return(object)
  
}


