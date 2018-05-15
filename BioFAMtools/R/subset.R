
################################
## Functions to do subsetting ##
################################

#' @title Subset factors
#' @name subsetFactors
#' @description Method to subset (or sort) factors
#' @param object a \code{\link{BioFAModel}} object.
#' @param factors character vector with the factor names, or numeric vector with the index of the factors.
#' @export

subsetFactors <- function(object, factors, keep_intercept=T) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  stopifnot(length(factors) <= object@Dimensions[["K"]])
  if (is.character(factors)) stopifnot(all(factors %in% factorNames(object)))
  


  # Get factors
  if(is.numeric(factors)) {
    if (object@ModelOptions$LearnIntercept == T) factors <- factorNames(object)[factors+1]
    else factors <- factorNames(object)[factors]
  } else { 
    stopifnot(all(factors %in% factorNames(object))) 
  }
  if (keep_intercept & object@ModelOptions$LearnIntercept & !"intercept" %in% factors) {
    factors <- c("intercept", factors)
  }
  
  # Subset expectations
  
  for (node in c("Z","W")){
    object@Expectations[[node]] <- sapply(object@Expectations[[node]], function(x) x[,factors], simplify = F, USE.NAMES = T)
  }
  
  for (node in c("ThetaZ","AlphaZ","ThetaW","AlphaW")){
    if (node %in% names(object@Expectations)){
      object@Expectations[[node]] <- sapply(object@Expectations[[node]], function(x) x[factors], simplify = F, USE.NAMES = T)
    }
  }
  
  # TODO for everything below : adapt for 2D
  
  # Reordering covariance hyperparameters and "spatial signifiances" (spatialFA)
  
  # Warning : below, does not work if factors are not convertible to integers
  if ("SigmaZ" %in% names(object@Expectations)){
    object@Expectations$SigmaZ <- object@Expectations$SigmaZ[,,as.integer(factors)]
  }
  
  # Warning : below, does not work if factors are not convertible to integers
  # Warning : below, does not work if SigmaAlphaW contains Alpha nodes in some views
  if ("SigmaAlphaW" %in% names(object@Expectations)){
    for (m in viewNames(object)){ 
      object@Expectations$SigmaAlphaW[[m]]<- object@Expectations$SigmaAlphaW[[m]][,,as.integer(factors)]
    }
  }
  
  if ("SigmaZ" %in% names(object@Parameters)){
    object@Parameters$SigmaZ$l <- object@Parameters$SigmaZ$l[factors]
    object@Parameters$SigmaZ$sig <- object@Parameters$SigmaZ$sig[factors]
  }
  
  # Warning : below, does not work if SigmaAlphaW contains Alpha nodes in some views
  if ("SigmaAlphaW" %in% names(object@Parameters)){
    for (m in viewNames(object)){ 
      object@Parameters$SigmaAlphaW[[m]]$l <- object@Parameters$SigmaAlphaW[[m]]$l[factors]
      object@Parameters$SigmaAlphaW[[m]]$sig <- object@Parameters$SigmaAlphaW[[m]]$sig[factors]
    }
  }
  
  # TODO : reorder parameters of other nodes
    
  # Modify dimensionality
  object@Dimensions[["K"]] <- length(factors)
  
  # Modify factor names
  factorNames(object) <- paste0("F", as.character(1:object@Dimensions[["K"]]))
  
  return(object)
}



#' @title Subset samples
#' @name subsetSamples
#' @description Method to subset (or sort) samples
#' @param object a \code{\link{BioFAModel}} object.
#' @param samples character vector with the sample names, numeric vector with the sample indices or logical vector with the samples to be kept as TRUE.
#' @export
subsetSamples <- function(object, samples) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  stopifnot(length(samples) <= object@Dimensions[["N"]])
  warning("Removing samples a posteriori is fine for an exploratory analysis, but we recommend removing them before training!")
  
  # Get samples
  if (is.character(samples)) {
    stopifnot(all(samples %in% sampleNames(object)))
  } else {
    samples <- sampleNames(object)[samples]
  }
  
  # Subset relevant slots
  object@Expectations$Z <- sapply(object@Expectations$Z, function(x) x[samples,, drop=F])
  object@Expectations$Y <- sapply(object@Expectations$Y, function(x) x[samples,], simplify = F, USE.NAMES = T)
  object@TrainData <- sapply(object@TrainData, function(x) sapply(x, function(y) y[,samples], simplify = F, USE.NAMES = T))
  object@InputData <- object@InputData[,samples,,] 
  if (length(object@ImputedData)==0) { object@ImputedData <- sapply(object@ImputedData, function(x) sapply(x, function(y) y[,samples], simplify = F, USE.NAMES = T)) }

  # Modify dimensionality
  object@Dimensions[["N"]] <- length(samples)
  
  # Modify sample names in the MOFAobject
  sampleNames(object) <- samples
  
  return(object)
}

