
##########################################################
## Functions to cluster samples based on latent factors ##
##########################################################

#' @title K-means clustering on samples based on latent factors
#' @name cluster_samples
#' @description MOFA factors are continuous in nature but they can be used to predict discrete clusters of samples. \cr
#' The clustering can be performed in a single factor, which is equivalent to setting a manual threshold.
#' More interestingly, it can be done using multiple factors, where multiple sources of variation are aggregated. \cr
#' Importantly, this type of clustering is not weighted and does not take into account the different importance of the latent factors. 
#' @param object a trained \code{\link{MOFA}} object.
#' @param k number of clusters (integer).
#' @param factors character vector with the factor name(s), or numeric vector with the index of the factor(s) to use. 
#' Default is 'all'
#' @param ... extra arguments  passed to \code{\link{kmeans}}
#' @details In some cases, due to model technicalities, samples can have missing values in the latent factor space. 
#' In such a case, these samples are currently ignored in the clustering procedure.
#' @return output from \code{\link{kmeans}} function
#' @importFrom stats kmeans
#' @export 
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Cluster samples in the factor space using factors 1 to 3 and K=2 clusters 
#' clusters <- cluster_samples(model, k=2, factors=1:3)
cluster_samples <- function(object, k, factors = "all", ...) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Define factors
  factors <- .check_and_get_factors(object, factors)
  
  # Collect relevant data
  Z <- get_factors(object, factors=factors)
  if (is(Z, "list")) Z <- do.call(rbind, Z)
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
