
#######################################################
## Functions to perform imputation of missing values ##
#######################################################

#' @title Impute missing values from a fitted BioFAModel
#' @name impute
#' @description This function uses the latent factors and the loadings to impute missing values.
#' @param object a \code{\link{BioFAM}} object.
#' @param views character vector with the view name(s), or numeric vector with view index(es).
#' @param groups character vector with the group name(s), or numeric vector with group index(es).
#' @param factors character vector with the factor names, or numeric vector with the factor index(es).
#' @param type type of imputation. 
#' \itemize{
#' \item \strong{response}:{gives mean for gaussian and poisson and probabilities for bernoulli.}
#' \item \strong{link}: {gives the linear predictions.}
#' \item \strong{inRange}: {rounds the fitted values from "terms" for integer-valued distributions to the next integer (default).}
#' }
#' @details matrix factorization models generate a denoised and condensed low-dimensional representation of the data which capture the main sources of heterogeneity of the data. 
#' These representation can be used to do predictions using the equation \code{Y = WX}. For more details read the supplementary methods of the manuscript. \cr
#' This method fills the \code{ImputedData} slot by replacing the missing values in the input data with the model predictions.
#' @export
impute <- function(object, views = "all", groups = "all", factors = "all", type = c("inRange", "response", "link")) {
  
  # Get views and groups
  views  <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)
  
  # Select imputation type  
  type <- match.arg(type)
  
  # Do predictions
  pred <- predict(object, views=views, factors=factors, type=type)

  # Replace NAs with predicted values
  # TO-DO: HOW TO REPLACE  VALUES WITH DELAYEDARRAYS??
  imputed <- get_data(object, views=views, groups=groups, add_intercept = FALSE)
  for (m in views) {
    for (g in groups) {
      non_observed <- which(is.na(imputed[[m]][[g]]), arr.ind = T)
      imputed[[m]][[g]][non_observed] <- pred[[m]][[g]][non_observed]
    }
  }

  # Save imputed data in the corresponding slot  
  object@imputed_data <- imputed
  
  return(object)
}
