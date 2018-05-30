
######################################
## Functions to perform predictions ##
######################################

#' @title Do predictions using a fitted BioFAModel
#' @name predict
#' @description This function uses the latent factors and the weights to do data predictions.
#' @param object a \code{\link{BioFAModel}} object.
#' @param views character vector with the view name(s), or numeric vector with the view index(es). 
#' Default is "all".
#' @param factors character vector with the factor name(s) or numeric vector with the factor index(es). 
#' Default is "all".
#' @param type type of prediction returned, either: 
#' "response" gives the response vector, the mean for Gaussian and Poisson, and probabilities for Bernoulli, 
#' "link" gives the linear predictions, 
#' "inRange" rounds the fitted values integer-valued distributions to the next integer. \cr
#' Default is "inRange".
#' @param include_intercept logical indicating whether to include the intercept factors for the prediction, if present.
#' Default is TRUE.
#' @details the denoised and condensed low-dimensional representation of the data captures the main sources of heterogeneity of the data. 
#' These representation can be used to do predictions using the equation Y = WX. This is the key step underlying imputation, see \code{\link{imputeMissing}} and Methods section of the article.
#' @return Returns a list with data predictions, each element corresponding to a view.
#' @export
predict <- function(object, views = "all", groups = "all", factors = "all", 
                    type = c("inRange", "response", "link"), 
                    include_intercept = TRUE) {

  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  
  # Get views  
  views <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)
  
  # Get factors
  if (paste0(factors, collapse="") == "all") { 
    factors <- factors_names(object) 
  } else if(is.numeric(factors)) {
      if (object@model_options$learn_intercept) factors <- factors_names(object)[factors+1]
      else factors <- factors_names(object)[factors]
  } else { 
    stopifnot(all(factors %in% factors_names(object))) 
  }

  # add/remove intercept factor for prediction
  if (include_intercept) {
    if (!"intercept" %in% factors & object@model_options$learn_intercept) {
      factors <- c("intercept", factors)  
    } else {
      factors <- factors[factors!="intercept"]
    }
    
  }
  
  # Get type of predictions wanted 
  type = match.arg(type)
  
  # Collect weights
  W <- get_weights(object, views = views, factors = factors)

  # Collect factors
  Z <- get_factors(object, groups = groups, factors = factors)
  Z[is.na(Z)] <- 0 # set missing values in Z to 0 to exclude from imputations
 
   # Coerce either W or Z to DelayedArrays and set HDF5array backend for on-disk operations
  if (object@on_disk) {
    Z <- lapply(Z, function(x) if (!is(x,"DelayedArray")) DelayedArray(x) else x )
    # W <- lapply(W, function(x) if (!is(x,"DelayedArray")) DelayedArray(x) else x )
    setRealizationBackend("HDF5Array")
  }
  
  # Predict data based on BioFAModel
  predictedData <- lapply(views, function(m) { lapply(groups, function(g) {
    
      # calculate terms based on linear model
      pred <- t(Z[[g]] %*% t(W[[m]])) 
      
      # make predicitons based on underlying likelihood
      lks <- object@model_options$likelihood
      names(lks) <- views_names(object)

      if (type != "link") {
        lk <- lks[m]
        if (lk == "gaussian") { 
          pred <- pred 
        }
        else if (lk == "bernoulli") { 
          pred <- (exp(pred)/(1 + exp(pred)))
          if (type == "inRange") pred <- round(pred)
        } else if (lk == "poisson") { 
          pred <- log(1 + exp(pred))
          if (type == "inRange") pred <- round(pred)
        }
        else { 
          stop(sprintf("Likelihood %s not implemented for imputation", lk)) 
        }
      }
      pred
      
    })
  })

  # RICARD: TO CHECK IF THIS IS OK IN DELAYEDARRAYS
  predictedData <- .name_views_and_groups(predictedData, views, groups)

  return(predictedData)
}