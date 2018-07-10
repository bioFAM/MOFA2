
##########################################################
## Functions to cluster samples based on latent factors ##
##########################################################

#' @title cluster_samples: K-means clustering on samples based on latent factors
#' @name cluster_samples
#' @description BioFAM factors are continuous in nature but they can be used to predict discrete clusters of samples, 
#' similar to the iCluster model (Shen, 2009). \cr
#' The clustering can be performed in a single factor, which is equivalent to setting a manual threshold; 
#' or using multiple factors, where multiple sources of variation are aggregated. \cr
#' Importantly, this type of clustering is not weighted and does not take into account the different importance of the latent factors. 
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param k number of clusters
#' @param factors character vector with the factor name(s), or numeric vector with the index of the factor(s) to use. 
#' Default is 'all'
#' @param ... extra arguments  passed to \code{\link{kmeans}}
#' @details In some cases, due to model technicalities, samples can have missing values in the latent factor space. 
#' In such a case, these samples are currently ignored in the clustering procedure.
#' @return output from \code{\link{kmeans}} function
#' @export
#' 
cluster_samples <- function(object, k, factors = "all", ...) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  

  # Define factors
  if (paste0(factors, collapse="") == "all") { factors <- factors_names(object) } 
  else if(is.numeric(factors)) {
    if (object@model_options$learn_intercept) factors <- factors_names(object)[factors+1]
    else factors <- factors_names(object)[factors]
  }
  else {
    stopifnot(all(factors %in% factors_names(object)))
  }
  
  # Collect relevant data
  Z <- get_factors(object, factors=factors, include_intercept=F)
  if (class(Z) == "list") Z <- do.call(rbind, Z)
  N <- nrow(Z)
  
  # For now remove sample with missing values on factors
  # (TO-DO) incorporate a clustering function that is able to cope with missing values
  haveAllZ <- apply(Z, 1, function(x) all(!is.na(x)))
  if(!all(haveAllZ)) warning(paste("Removing", sum(!haveAllZ), "samples with missing values on at least one factor"))
  Z <- Z[haveAllZ,]
  # Perform k-means clustering
  kmeans.out <- kmeans(Z, centers=k,  ...)

  return(kmeans.out)  

}
