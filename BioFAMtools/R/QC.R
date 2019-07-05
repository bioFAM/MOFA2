#' @title qualityControl
#' @name qualityControl
#' @description Function to do quality control on a \code{\link{BioFAModel}} object. \cr
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param verbose logical indicating whether to generate a verbose output.
#' @export
#'
qualityControl <- function(object, verbose = F) {
  
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  
  # Check views names
  if (verbose==T) message("Checking views names...")
  stopifnot(!is.null(views_names(object)))
  stopifnot(!duplicated(views_names(object)))
  
  # Check groups names
  if (verbose==T) message("Checking groups names...")
  stopifnot(!is.null(groups_names(object)))
  stopifnot(!duplicated(groups_names(object)))
  
  if (object@Status == "untrained") {
    
    # Check dimensionalities in the input data
    N <- object@dimensions$N
    D <- object@dimensions$D
    for (i in views_names(object)) {
      for (j in groups_names(object)) {
        stopifnot(nrow(object@data[[i]][[j]]) == N[[j]])
        stopifnot(ncol(object@data[[i]][[j]]) == D[[i]])
        stopifnot(length(rownames(object@data[[i]][[j]])) == N[[j]])
        stopifnot(length(colnames(object@data[[i]][[j]])) == D[[i]])
      }
    }

    # Check that there are no features with zero variance
    # if (verbose==T) message("Checking there are no features with zero variance...")
    # for (i in views_names(object)) {
    #   for (j in groups_names(object)) {
    #   if (!all(apply(object@data[[i]][[j]],1,var,na.rm=T) > 0, na.rm=T)) {
    #     print("Warning: you have features with zero variance, consider removing them...")
    #   }
    # }
    
    # Check that there are no features with complete missing values
    # if (verbose==T) message("Checking there are no features with complete missing values...")
    # for (i in views_names(object)) {
    #   for (j in groups_names(object)) {
    #   if (!all(apply(object@data[[i]][[j]],1, function(x) mean(is.na(x)))<1, na.rm=T)) {
    #     print("Warning: you have features which do not contain any observation, consider removing them...")
    #   }
    # }
    
  } else if (object@Status == "trained") {
    
    # Check dimensionalities in the input data
    N <- object@dimensions$N
    D <- object@dimensions$D
    for (i in views_names(object)) {
      for (j in groups_names(object)) {
        stopifnot(ncol(object@data[[i]][[j]]) == N[[j]])
        stopifnot(nrow(object@data[[i]][[j]]) == D[[i]])
        stopifnot(length(colnames(object@data[[i]][[j]])) == N[[j]])
        stopifnot(length(rownames(object@data[[i]][[j]])) == D[[i]])
      }
    }
    
    # Check samples names
    if (verbose==T) message("Checking samples names...")
    stopifnot(!is.null(samples_names(object)))
    stopifnot(!duplicated(unlist(samples_names(object))))
    
    # Check features names
    if (verbose==T) message("Checking features names...")
    stopifnot(!is.null(features_names(object)))
    stopifnot(!duplicated(unlist(features_names(object))))
    
    # Check expectations
    # stopifnot(identical(sort(c("W","Z","Theta","Tau","AlphaW","Y")), sort(names(object@Expectations))))
    stopifnot(all(c("W","Z") %in% names(object@expectations)))
    if (verbose==T) message("Checking expectations...")
    stopifnot(all(sapply(object@expectations$W, is.matrix)))
    stopifnot(all(sapply(object@expectations$Z, is.matrix)))

    # Check that the likelihoods match the data distribution
    # if (verbose==T) message("Checking likelihooods...")
    # predicted_lik <- .inferLikelihoods(object)
    # for (view in viewNames(object)) {
    #   lk <- object@ModelOptions$likelihood[view]
    #   if (lk != predicted_lik[view])
    #     message(sprintf("Warning, view %s should follow a %s distribution rather than %s ", view, predicted_lik[view], lk))
    # }
    
  }

}
