#' @title qualityControl
#' @name qualityControl
#' @description Function to do quality control on a \code{\link{BioFAModel}} object. \cr
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param verbose logical indicating whether to generate a verbose output.
#' @export
#'
qualityControl <- function(object, verbose = F) {
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  if (object@Status != "trained") stop("This function only works in a trained BioFAModel")
  
  # Check that the model has view names
  if (verbose==T) message("Checking view names...")
  stopifnot(!is.null(viewNames(object)))
  
  # Check that the model has sample names
  if (verbose==T) message("Checking sample names...")
  stopifnot(!is.null(sampleNames(object)))
  
  # Check that the model has feature names
  if (verbose==T) message("Checking feature names...")
  stopifnot(!is.null(featureNames(object)))

  # DEPRECATED
  # Check that the model has the right node names
  # if (verbose==T) message("Checking nodes...")
  # stopifnot(identical(sort(c("W","Z","Theta","Tau","AlphaW","Y")), sort(names(object@Expectations))))
  
  # Check that all expectations are the correct object
  if (verbose==T) message("Checking expectations...")
  stopifnot(is.matrix(object@Expectations$Z))
  stopifnot(is.list(object@Expectations$W))
  stopifnot(all(sapply(object@Expectations$W, is.matrix)))
  stopifnot(is.list(object@Expectations$Y))
  stopifnot(all(sapply(object@Expectations$Y, is.matrix)))
  # stopifnot(is.list(object@Expectations$Theta))
  # stopifnot(all(sapply(object@Expectations$Theta, is.matrix)))
  stopifnot(is.list(object@Expectations$Tau))
  stopifnot(all(sapply(object@Expectations$Tau, is.numeric)))
  stopifnot(is.list(object@Expectations$AlphaW))
  stopifnot(all(sapply(object@Expectations$AlphaW, is.numeric)))
  
  # Check that the dimensionalities match
  # TO-DO...
  if (verbose==T) message("Checking dimensionalities...")
  
  # Check that there are no features with complete missing values
  if (verbose==T) message("Checking there are no features with complete missing values...")
  for (view in viewNames(object)) {
    # FIX THIS
    if (!all(apply(object@TrainData[[view]],1, function(x) mean(is.na(x))) < 1, na.rm=T)) {
      print("Warning: you have features which only contain missing values, consider removing them...")
    }
  }
  
  # Check that there are no features with zero variance
  if (verbose==T) message("Checking there are no features with zero variance...")
  for (view in viewNames(object)) {
    if (!all(apply(object@TrainData[[view]],1,var,na.rm=T) > 0, na.rm=T)) {
      print("Warning: you have features with zero variance, consider removing them...")
    }
  }
  
  # Check that the likelihoods match the data distribution
  if (verbose==T) message("Checking likelihooods...")
  predicted_lik <- .inferLikelihoods(object)
  for (view in viewNames(object)) {
    lk <- object@ModelOptions$likelihood[view]
    if (lk != predicted_lik[view])
      message(sprintf("Warning, view %s should follow a %s distribution rather than %s ", view, predicted_lik[view], lk))
  }
  
}
