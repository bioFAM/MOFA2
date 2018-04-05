
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
predict <- function(object, views = "all", factors = "all", 
                    type = c("inRange","response", "link"), 
                    include_intercept = TRUE) {

  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  
  # Get views  
  if (is.numeric(views)) {
    stopifnot(all(views<=object@Dimensions$M))
    views <- viewNames(object)[views] 
  } else {
    if (paste0(views,sep="",collapse="") =="all") { 
      views = viewNames(object)
    } else {
      stopifnot(all(views%in%viewNames(object)))
    }
  }
  
  # Get factors
  if (paste0(factors,collapse="") == "all") { 
    factors <- factorNames(object) 
  } else if(is.numeric(factors)) {
      if (object@ModelOpts$learnIntercept == T) factors <- factorNames(object)[factors+1]
      else factors <- factorNames(object)[factors]
  } else { 
    stopifnot(all(factors %in% factorNames(object))) 
  }

  # add intercept factor for prediction
  if(!"intercept" %in% factors & object@ModelOpts$learnIntercept & include_intercept) factors <- c("intercept", factors)  
  if(!include_intercept & "intercept" %in% factors) factors <- factors[factors!="intercept"]
  
  # Get type of predictions wanted 
  type = match.arg(type)
  
  # Collect weights
  W <- getWeights(object, views=views, factors=factors)

  # Collect factors
  Z <- getFactors(object)[,factors]
  Z[is.na(Z)] <- 0 # set missing values in Z to 0 to exclude from imputations
 
  # Predict data based on BioFAModel
  # predictedData <- lapply(sapply(views, grep, viewNames(object)), function(viewidx){
  predictedData <- lapply(views, function(i){
    
    # calculate terms based on linear model
    predictedView <- t(Z%*% t(W[[i]])) 
    
    # make predicitons based on underlying likelihood
    if (type!="link") {
      lk <- object@ModelOpts$likelihood[i]
      if (lk == "gaussian") { 
        predictedView <- predictedView 
      }
      else if (lk == "bernoulli") { 
        predictedView <- (exp(predictedView)/(1+exp(predictedView)))
        if (type=="inRange") predictedView <- round(predictedView)
      } else if (lk == "poisson") { 
        # predictedView <- (exp(predictedView))
        predictedView <- log(1 + exp(predictedView))
        if(type=="inRange") predictedView <- round(predictedView)
      }
      else { 
        stop(sprintf("Likelihood %s not implemented for imputation",lk)) 
      }
    }
    predictedView
  })

  names(predictedData) <- views

  return(predictedData)
}