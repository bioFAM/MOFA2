
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
  object@Expectations$Z <- object@Expectations$Z[,factors]
  object@Expectations$AlphaW <- sapply(object@Expectations$AlphaW, function(x) x[factors], simplify = F, USE.NAMES = T)
  object@Expectations$W <- sapply(object@Expectations$W, function(x) x[,factors], simplify = F, USE.NAMES = T)
  object@Expectations$Theta <- sapply(object@Expectations$Theta, function(x) x[factors], simplify = F, USE.NAMES = T)

  # Modify dimensionality
  object@Dimensions[["K"]] <- length(factors)
  
  # Modify factor names
  factorNames(object) <- as.character(1:object@Dimensions[["K"]])
  
  return(object)
}

