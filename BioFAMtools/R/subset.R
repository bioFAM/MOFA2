
################################
## Functions to do subsetting ##
################################

#' @title Subset factors
#' @name subsetFactors
#' @description Method to subset (or sort) factors
#' @param object a \code{\link{BioFAModel}} object.
#' @param factors character vector with the factor names, or numeric vector with the index of the factors.
#' @export

subsetFactors <- function(object, factors) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  stopifnot(length(factors) <= object@Dimensions[["K"]])
  if (is.character(factors)) stopifnot(all(factors %in% factorNames(object)))
  if (is.numeric(factors)) stopifnot(all(factors %in% 1:object@Dimensions[["K"]]))
  
  # Subset expectations
  object@Expectations$Z <- sapply(object@Expectations$Z, function(x) x[,factors], simplify = F, USE.NAMES = T)
  object@Expectations$W <- sapply(object@Expectations$W, function(x) x[,factors], simplify = F, USE.NAMES = T)
  
  nodes=c("AlphaZ","SigmaZ","ThetaZ","AlphaW","SigmaAlphaW","ThetaW")
  for (node in nodes){
    if (node %in% names(object@Expectations)){
      object@Expectations$node <- sapply(object@Expectations$node, function(x) x[factors], simplify = F, USE.NAMES = T)
    }
  }
  
  # TODO: adapt for 2D
  # Reordering covariance hyperparameters and "spatial signifiances" (spatialFA)

  if ("SigmaZ" %in% names(object@Parameters)){
    object@Parameters$SigmaZ$l <- object@Parameters$SigmaZ$l[factors]
    object@Parameters$SigmaZ$sig <- object@Parameters$SigmaZ$sig[factors]
  }

  if ("SigmaAlphaW" %in% names(object@Parameters)){
    for (m in viewNames(object)){ 
    object@Parameters$SigmaAlphaW[[m]]$l <- object@Parameters$SigmaAlphaW[[m]]$l[factors]
    object@Parameters$SigmaAlphaW[[m]]$sig <- object@Parameters$SigmaAlphaW[[m]]$sig[factors]
    }
  }
  
  # TODO : reorder parameters of other nodes ?
    
  # Modify dimensionality
  object@Dimensions[["K"]] <- length(factors)
  
  # Modify factor names
  factorNames(object) <- paste0("F", as.character(1:object@Dimensions[["K"]]))
  
  return(object)
}


