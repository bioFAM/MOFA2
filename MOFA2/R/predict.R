
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
#' @param type type of prediction returned, either:
#' "response" gives the response vector, the mean for Gaussian and Poisson, and probabilities for Bernoulli,
#' "link" gives the linear predictions,
#' "inRange" rounds the fitted values integer-valued distributions to the next integer (default). \cr
#' @details the denoised and condensed low-dimensional representation of the data captures the main sources of heterogeneity of the data.
#' These representation can be used to do predictions using the equation Y = WX. This is the key step underlying imputation, see \code{\link{impute}} and Methods section of the article.
#' @return Returns a list with data predictions, each element corresponding to a view.
#' @export
predict <- function(object, views = "all", groups = "all", factors = "all",
                    type = c("inRange", "response", "link")) {

  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")

  # Get views
  views <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)

  # Get factors
  if (paste0(factors, collapse="") == "all") {
    factors <- factors_names(object)
  } else if (is.numeric(factors)) {
    factors <- factors_names(object)[factors]
  } else {
    stopifnot(all(factors %in% factors_names(object)))
  }

  # Get type of predictions wanted
  type = match.arg(type)

  # Collect weights
  W <- get_weights(object, views = views, factors = factors)

  # Collect factors
  Z <- get_factors(object, groups = groups, factors = factors)
  Z[is.na(Z)] <- 0 # set missing values in Z to 0 to exclude from imputations

  # Do predictions
  predicted_data <- lapply(views, function(m) { lapply(groups, function(g) {

      # calculate terms based on linear model
      pred <- t(Z[[g]] %*% t(W[[m]]))

      # make predicitons based on underlying likelihood
      lks <- object@model_options$likelihoods

      if (type != "link") {
        lk <- lks[m]

        if (lk == "gaussian") {
          pred <- pred
        } else if (lk == "bernoulli") {
          pred <- (exp(pred)/(1 + exp(pred)))
          if (type == "inRange") pred <- round(pred)
        } else if (lk == "poisson") {
          pred <- log(1 + exp(pred))
          if (type == "inRange") pred <- round(pred)
        }
      }
      pred
    })
  })

  predicted_data <- .name_views_and_groups(predicted_data, views, groups)

  return(predicted_data)
}
