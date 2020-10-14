
######################################
## Functions to perform predictions ##
######################################

#' @title Do predictions using a fitted MOFA
#' @name predict
#' @description This function uses the latent factors and the weights to do data predictions.
#' @param object a \code{\link{MOFA}} object.
#' @param views character vector with the view name(s), or numeric vector with the view index(es).
#' Default is "all".
#' @param groups character vector with the group name(s), or numeric vector with the group index(es).
#' Default is "all".
#' @param factors character vector with the factor name(s) or numeric vector with the factor index(es).
#' Default is "all".
# #' @param type type of prediction returned, either:
# #' "response" gives the response vector, the mean for Gaussian and Poisson, and probabilities for Bernoulli,
# #' "link" gives the linear predictions,
# #' "inRange" rounds the fitted values integer-valued distributions to the next integer (default).
#' @param add_intercept add feature intercepts to the prediction (default is TRUE).
#' @details MOFA generates a denoised and condensed low-dimensional representation of the data that captures the main sources of heterogeneity of the data.
#' This representation can be used to reconstruct a denoised representation of the data, simply using the equation \code{Y = WX}. 
#' For more mathematical details read the supplementary methods of the manuscript.
#' @return Returns a list with the data reconstructed by the model predictions.
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Predict observations for all data modalities
#' predictions <- predict(model)
predict <- function(object, views = "all", groups = "all", factors = "all", add_intercept = TRUE) {

  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get views
  views <- .check_and_get_views(object, views, non_gaussian=FALSE)
  groups <- .check_and_get_groups(object, groups)

  # Sanity check
  if (any(views %in% names(which(object@model_options$likelihoods!="gaussian")))) stop("predict does not work for non-gaussian modalities")
  
  # Get factors
  if (paste0(factors, collapse="") == "all") {
    factors <- factors_names(object)
  } else if (is.numeric(factors)) {
    factors <- factors_names(object)[factors]
  } else {
    stopifnot(all(factors %in% factors_names(object)))
  }

  # Collect weights
  W <- get_weights(object, views = views, factors = factors)

  # Collect factors
  Z <- get_factors(object, groups = groups, factors = factors)
  Z[is.na(Z)] <- 0 # set missing values in Z to 0 to exclude from imputations

  # Do predictions
  predicted_data <- lapply(views, function(m) { lapply(groups, function(g) {

      # calculate terms based on linear model
      pred <- t(Z[[g]] %*% t(W[[m]]))

      # add feature-wise intercepts (i think this does not work for non-gaussian likelihhood, needs some verification)
      tryCatch( {
        if (add_intercept & length(object@intercepts[[1]])>0) {
          intercepts <- object@intercepts[[m]][[g]]
          intercepts[is.na(intercepts)] <- 0
          pred <- pred + object@intercepts[[m]][[g]]
        } }, error = function(e) { NULL })

      return(pred)
    })
  })

  predicted_data <- .name_views_and_groups(predicted_data, views, groups)

  return(predicted_data)
}
