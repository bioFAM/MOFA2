
#######################################################
## Functions to perform imputation of missing values ##
#######################################################

#' @title Impute missing values from a fitted MOFA
#' @name impute
#' @description This function uses the latent factors and the loadings to impute missing values.
#' @param object a \code{\link{MOFA}} object.
#' @param views character vector with the view name(s), or numeric vector with view index(es).
#' @param groups character vector with the group name(s), or numeric vector with group index(es).
#' @param factors character vector with the factor names, or numeric vector with the factor index(es).
#' @param type type of imputation.
#' \itemize{
#' \item \strong{response}:{gives mean for gaussian and poisson and probabilities for bernoulli.}
#' \item \strong{link}: {gives the linear predictions.}
#' \item \strong{inRange}: {rounds the fitted values from "terms" for integer-valued distributions to the next integer (default).}
#' }
#' @details MOFA generates a denoised and condensed low-dimensional representation of the data that captures the main sources of heterogeneity of the data.
#' This representation can be used to reconstruct the data, simply using the equation \code{Y = WX}. 
#' For more details read the supplementary methods of the manuscript. \cr
#' Note that with \code{\link{impute}} you can only generate the point estimates (the means of the posterior distributions). 
#' If you want to add uncertainity estimates (the variance) you need to set \code{impute=TRUE} in the training options.
#' See \code{\link{get_default_training_options}}.
#' @return This method fills the \code{imputed_data} slot by replacing the missing values in the input data with the model predictions.
#' @export
impute <- function(object, views = "all", groups = "all", factors = "all", type = c("inRange", "response", "link")) {

  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (length(object@imputed_data)>0) warning("imputed_data slot is already filled. It will be replaced")
  
  # Get views and groups
  views  <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)

  # Select imputation type
  type <- match.arg(type)

  # Do predictions
  pred <- predict(object, views=views, factors=factors, type=type)

  # Replace NAs with predicted values
  imputed <- get_data(object, views=views, groups=groups, add_intercept = FALSE)
  for (m in views) {
    for (g in groups) {
      imputed[[m]][[g]] <- list("mean" = imputed[[m]][[g]])
      non_observed <- which(is.na(imputed[[m]][[g]]), arr.ind = T)
      imputed[[m]][[g]]$mean[non_observed] <- pred[[m]][[g]][non_observed]
      imputed[[m]][[g]]$variance <- list() # empty variance
    }
  }
  
  # Save imputed data in the corresponding slot
  object@imputed_data <- imputed

  return(object)
}

