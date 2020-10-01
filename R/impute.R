
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
#' \itemize{
#' \item \strong{response}:{gives mean for gaussian and poisson and probabilities for bernoulli.}
#' \item \strong{link}: {gives the linear predictions.}
#' \item \strong{inRange}: {rounds the fitted values from "terms" for integer-valued distributions to the next integer (default).}
#' }
#' @param add_intercept add feature intercepts to the imputation (default is TRUE).
#' @details MOFA generates a denoised and condensed low-dimensional representation of the data that captures the main sources of heterogeneity of the data.
#' This representation can be used to reconstruct the data, simply using the equation \code{Y = WX}. 
#' For more details read the supplementary methods of the manuscript. \cr
#' Note that with \code{\link{impute}} you can only generate the point estimates (the means of the posterior distributions). 
#' If you want to add uncertainity estimates (the variance) you need to set \code{impute=TRUE} in the training options.
#' See \code{\link{get_default_training_options}}.
#' @return This method fills the \code{imputed_data} slot by replacing the missing values in the input data with the model predictions.
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Impute missing values in all data modalities
#' imputed_data <- impute(model, views = "all")
#' 
#' # Impute missing values in all data modalities using factors 1:3
#' imputed_data <- impute(model, views = "all", factors = 1:3)
impute <- function(object, views = "all", groups = "all", factors = "all", 
                  add_intercept = TRUE) {

  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (length(object@imputed_data)>0) warning("imputed_data slot is already filled. It will be replaced and the variance estimates will be lost...")
  
  # Get views and groups
  views  <- .check_and_get_views(object, views, non_gaussian=FALSE)
  groups <- .check_and_get_groups(object, groups)


  # Do predictions
  pred <- predict(object, views=views, factors=factors, add_intercept=add_intercept)

  # Replace NAs with predicted values
  imputed <- get_data(object, views=views, groups=groups, add_intercept = add_intercept)
  for (m in views) {
    for (g in groups) {
      imputed[[m]][[g]] <- imputed[[m]][[g]]
      non_observed <- is.na(imputed[[m]][[g]])
      imputed[[m]][[g]][non_observed] <- pred[[m]][[g]][non_observed]
    }
  }
  
  # Save imputed data in the corresponding slot
  object@imputed_data <- imputed

  return(object)
}

