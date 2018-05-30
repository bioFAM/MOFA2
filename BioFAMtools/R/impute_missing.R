
#######################################################
## Functions to perform imputation of missing values ##
#######################################################

#' @title Impute missing values from a fitted BioFAModel
#' @name impute_missing
#' @description This function uses the latent factors and the loadings infered in order to impute missing values.
#' @param object a \code{\link{BioFAM}} object.
#' @param views character vector with the view names, or numeric vector with view indexes.
#' @param factors character vector with the factor names, or numeric vector with the factor indexes.
#' @param type type of imputation. 
#' "response" gives mean for gaussian and poisson and probabilities for bernoulli,
#' "link" gives the linear predictions,
#' "inRange" (default) rounds the fitted values from "terms" for integer-valued distributions to the next integer.
#' @details matrix factorization models generate a denoised and condensed low-dimensional representation of the data which capture the main sources of heterogeneity of the data. 
#' These representation can be used to do predictions using the equation \code{Y = WX}. For more details read the supplementary methods of the manuscript. \cr
#' This method fills the \code{ImputedData} slot by replacing the missing values in the input data with the model predictions.
#' @export
impute_missing <- function(object, views = "all", groups = "all", factors = "all", type = c("inRange", "response", "link")) {
  
  # Get views and groups
  views  <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)
  
  # Select imputation type  
  type <- match.arg(type)
  
  # Do predictions
  predData <- predict(object, views = views, factors = factors, type = type)

  # Replace NAs with predicted values
  # RICARD: HOW TO REPLACE SUBSET VALUES WITH DELAYEDARRAYS??
  imputedData <- get_training_data(object, views = views, groups = groups)
  for (m in views) {
    for (g in groups) {
      non_observed <- which(is.na(imputedData[[m]][[g]]), arr.ind = T)
      imputedData[[m]][[g]][non_observed] <- predData[[m]][[g]][non_observed]
    }
  }

  # re- arrange list in accordance with other data slots in the model
  # names(imputedData) <- views
  # imputedData <- imputedData[views_names(object)]
  # names(imputedData) <- views_names(object)
  
  # Save imputed data in the corresponding slot  
  object@imputed_data <- imputedData
  
  return(object)
}
