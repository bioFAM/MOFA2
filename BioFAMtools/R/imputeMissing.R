
#######################################################
## Functions to perform imputation of missing values ##
#######################################################

#' @title Impute missing values from a fitted BioFAModel
#' @name imputeMissing
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
imputeMissing <- function(object, views = "all", factors = "all", type = c("inRange","response", "link")) {
  
  # Get views  
  if (paste0(views,sep="",collapse="") =="all") { 
    views = viewNames(object)
  } else {
    stopifnot(all(views%in%viewNames(object)))
  }
  
  # Select imputation type  
  type = match.arg(type)
  
  # Do predictions
  predData <- predict(object, views=views, factors = factors, type = type)

  # WHY IS THIS COMMENTED? CHECK
  #replace NAs with predicted values
  # imputedData <- getTrainData(object, views = views)
  # imputedData <- lapply(names(imputedData), function(viewnm) {
  #     view <- imputedData[[viewnm]]
  #     non_observed <- which(is.na(view), arr.ind = T)
  #     if(viewnm %in% names(predData)) view[non_observed] <- predData[[viewnm]][non_observed]
  #     view
  # })

  # re- arrange list in accordance with other data slots in the model
  # names(imputedData) <- views
  # imputedData <- imputedData[viewNames(object)]
  # names(imputedData) <- viewNames(object)

  imputedData <- predData

  # Save imputed data in the corresponding slot  
  object@ImputedData <- imputedData
  
  return(object)
}
