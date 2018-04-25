
################################################
## Get functions to fetch data from the model ##
################################################

#' @title getDimensions 
#' @name getDimensions
#' @description Extract dimensionalities from the model. 
#' @details K indicates the number of factors, D indicates the number of features, 
#' N indicates the (total) number of samples and M indicates the number of views.
#' @param object a \code{\link{BioFAModel}} object.
#' @export
getDimensions <- function(object) {
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")  
  return(object@Dimensions)
}


#' @title getFactors
#' @name getFactors
#' @description Extract the latent factors from the model.
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param factors character vector with the factor name(s), or numeric vector with the factor index(es).
#' Default is "all".
#' @param include_intercept logical indicating where to include the intercept term of the model, if present. 
#' Default is \code{TRUE}.
#' @param as.data.frame logical indicating whether to return a long data frame instead of a matrix.
#' Default is \code{FALSE}.
#' @return By default it returns the latent factor matrix of dimensionality (N,K), where N is number of samples and K is number of factors. \cr
#' Alternatively, if \code{as.data.frame} is \code{TRUE}, returns a long-formatted data frame with columns (sample,factor,value).
#' @export
#' 
getFactors <- function(object, factors = "all", as.data.frame = FALSE, include_intercept = TRUE) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  
  # Get factors
  if (paste0(factors, collapse="") == "all") { factors <- factorNames(object) } 
  else if(is.numeric(factors)) {
    if (object@ModelOpts$learnIntercept == T) factors <- factorNames(object)[factors+1]
    else factors <- factorNames(object)[factors]
  }
  else { stopifnot(all(factors %in% factorNames(object))) }

  # Collect factors
  Z <- getExpectations(object, "Z", as.data.frame)
  if (as.data.frame == FALSE) {
    Z <- Z[,factors, drop=FALSE]
  } else {
    Z <- Z[Z$factor %in% factors,]
  }

  # Remove intercept
  if (include_intercept == FALSE) {
    if (as.data.frame==FALSE) {
      if ("intercept" %in% colnames(Z)) Z <- Z[,colnames(Z)!="intercept"]
    } else {
      if ("intercept" %in% unique(Z$factor)) Z <- Z[Z$factor!="intercept",]
    }
  }
  return(Z)
}


#' @title getWeights
#' @name getWeights
#' @description Extract the weights from the model.
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param views character vector with the view name(s), or numeric vector with the view index(es). 
#' Default is "all".
#' @param factors character vector with the factor name(s) or numeric vector with the factor index(es). \cr
#' Default is "all".
#' @param as.data.frame logical indicating whether to return a long data frame instead of a list of matrices. 
#' Default is \code{FALSE}.
#' @return By default it returns a list where each element is a loading matrix with dimensionality (D,K), 
#' where D is the number of features and K is the number of factors. \cr
#' Alternatively, if \code{as.data.frame} is \code{TRUE}, returns a long-formatted data frame with columns (view,feature,factor,value).
#' @export
#' 
getWeights <- function(object, views = "all", factors = "all", as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  
  # Get views and factors
  if (paste0(views,collapse="") == "all") { views <- viewNames(object) } else { stopifnot(all(views %in% viewNames(object))) }

  # Get factors
  if (paste0(factors,collapse="") == "all") { 
    factors <- factorNames(object) 
  } else if (is.numeric(factors)) {
      if (object@ModelOpts$learnIntercept == T) {
        factors <- factorNames(object)[factors+1]
      } else {
        factors <- factorNames(object)[factors]
      }
    } else { stopifnot(all(factors %in% factorNames(object))) }
        
  # Fetch weights
  weights <- getExpectations(object,"W",as.data.frame)
  if (as.data.frame==T) {
    weights <- weights[weights$view%in%views & weights$factor%in%factors, ]
  } else {
    weights <- lapply(views, function(m) weights[[m]][,factors,drop=F])
    names(weights) <-  views
    # if (length(views)==1) { weights <- weights[[1]] }
  }
  return(weights)
}


#' @title getTrainData
#' @name getTrainData
#' @description Fetch the training data
#' @param object a \code{\link{BioFAModel}} object.
#' @param views character vector with the view name(s), or numeric vector with the view index(es). 
#' Default is "all".
#' @param features list of character vectors with the feature names or list of numeric vectors with the feature indices. 
#' Default is "all"
#' @param as.data.frame logical indicating whether to return a long data frame instead of a list of matrices.
#' Default is \code{FALSE}.
#' @details By default this function returns a list where each element is a data matrix with dimensionality (D,N) 
#' where D is the number of features and N is the number of samples. \cr
#' Alternatively, if \code{as.data.frame} is \code{TRUE}, the function returns a long-formatted data frame with columns (view,feature,sample,value).
#' @export
getTrainData <- function(object, views = "all", features = "all", as.data.frame = F) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  
  # Get views
  if (paste0(views,collapse="") == "all") { views <- viewNames(object) } else { stopifnot(all(views %in% viewNames(object))) }
  
  # Get features
  if (class(features)=="list") {
    stopifnot(all(sapply(1:length(features), function(i) all(features[[i]] %in% featureNames(object)[[views[i]]]))))
  } else {
    if (paste0(features,collapse="") == "all") { 
      features <- featureNames(object)[views]
    } else {
      stop("features not recognised, please read the documentation")
    }
  }
  
  # Fetch data
  trainData <- object@TrainData[views]
  trainData <- lapply(1:length(trainData), function(m) trainData[[m]][features[[m]],,drop=F]); names(trainData) <- views
  
  # Convert to long data frame
  if (as.data.frame==T) {
    tmp <- lapply(views, function(m) { tmp <- reshape2::melt(trainData[[m]]); colnames(tmp) <- c("feature","sample","value"); tmp <- cbind(view=m,tmp); return(tmp) })
    trainData <- do.call(rbind,tmp)
    trainData[,c("view","feature","sample")] <- sapply(trainData[,c("view","feature","sample")], as.character)
  }# else if ((length(views)==1) && (as.data.frame==F)) {
  #  trainData <- trainData[[views]]
  #}
  
  return(trainData)
}


#' @title getImputedData
#' @name getImputedData
#' @description Function to get the imputed data. It requires the previous use of the \code{\link{imputeMissing}} method.
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param views character vector with the view name(s), or numeric vector with the view index(es). 
#' Default is "all".
#' @param features list of character vectors with the feature names or list of numeric vectors with the feature indices. 
#' Default is "all"
#' @param as.data.frame logical indicating whether to return a long-formatted data frame instead of a list of matrices. 
#' Default is \code{FALSE}.
#' @return By default returns a list where each element is a matrix with dimensionality (D,N), where D is the number of features in this view and N is the number of samples. \cr
#' Alternatively, if \code{as.data.frame} is \code{TRUE}, returns a long-formatted data frame with columns (view,feature,sample,value).
#' @export
getImputedData <- function(object, views = "all", features = "all", as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  
  # Get views
  if (paste0(views,collapse="") == "all") { views <- viewNames(object) } else { stopifnot(all(views %in% viewNames(object))) }
  
  # Get features
  if (class(features)=="list") {
    stopifnot(all(sapply(1:length(features), function(i) all(features[[i]] %in% featureNames(object)[[views[i]]]))))
  } else {
    if (paste0(features,collapse="") == "all") { 
      features <- featureNames(object)[views]
    } else {
      stop("features not recognised, please read the documentation")
    }
  }
  
  # Fetch data
  ImputedData <- object@ImputedData[views]
  ImputedData <- lapply(1:length(ImputedData), function(m) ImputedData[[m]][features[[m]],,drop=F]); names(ImputedData) <- views
  
  # Convert to long data frame
  if (as.data.frame==T) {
    tmp <- lapply(views, function(m) { tmp <- reshape2::melt(ImputedData[[m]]); colnames(tmp) <- c("feature","sample","value"); tmp <- cbind(view=m,tmp); return(tmp) })
    ImputedData <- do.call(rbind,tmp)
    ImputedData[,c("view","feature","sample")] <- sapply(ImputedData[,c("view","feature","sample")], as.character)
  } else if ((length(views)==1) && (as.data.frame==F)) {
    ImputedData <- ImputedData[[views]]
  }
  
  return(ImputedData)
}

#' @name getCovariates
#' @title getCovariates
#' @description This function extracts covariates from the \code{colData} in the input \code{MultiAssayExperiment} object. \cr
#' Note that if you did not use \code{MultiAssayExperiment} to create your \code{\link{createBioFAMobject}}, this function will not work.
#' @param object a \code{\link{BioFAModel}} object.
#' @param covariates names of the covariates
#' @export
#' 
getCovariates <- function(object, covariates) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  if(class(object@InputData) != "MultiAssayExperiment") stop("To work with covariates, InputData has to be specified in form of a MultiAssayExperiment")  
  stopifnot(all(covariates %in% colnames(colData(object@InputData))))
  
  # Get covariates
  covariates <- colData(object@InputData)[,covariates]
  
  return(covariates)
}

#' @title getExpectations
#' @name getExpectations
#' @description Function to extract the expectations from the (variational) posterior distributions of a trained \code{\link{BioFAModel}} object.
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param variable variable name: 'Z' for factors, 'W' for weights, 'Tau' for noise,
#' 'Y' for pseudodata, 'Theta' for feature-wise spike-and-slab sparsity, 'AlphaW' for view and factor-wise ARD sparsity
#' @param as.data.frame logical indicating whether to output the result as a long data frame, default is \code{FALSE}.
#' @details Technical note: BioFAM is a Bayesian model where each variable has a prior distribution and a posterior distribution. 
#' In particular, to achieve scalability we used the variational inference framework, thus true posterior distributions are replaced by approximated variational distributions.
#' This function extracts the expectations of the variational distributions, which can be used as final point estimates to analyse the results of the model. \cr 
#' The priors and variational distributions of each variable are extensively described in the supplementary methods of the original paper.
#' @return the output varies depending on the variable of interest: \cr
#' \itemize{
#'  \item{"Z"}{a matrix with dimensions (samples,factors). If \code{as.data.frame} is \code{TRUE}, a long-formatted data frame with columns (sample,factor,value)}
#'  \item{"W"}{a list of length (views) where each element is a matrix with dimensions (features,factors). If \code{as.data.frame} is \code{TRUE}, a long-formatted data frame with columns (view,feature,factor,value)}
#'  \item{"Y"}{a list of length (views) where each element is a matrix with dimensions (features,samples). If \code{as.data.frame} is \code{TRUE}, a long-formatted data frame with columns (view,feature,sample,value)}
#'  \item{"Theta"}{}
#'  \item{"Tau"}{}
#' }
#' @export
getExpectations <- function(object, variable, as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  stopifnot(variable %in% names(object@Expectations))
  
  # Get expectations in single matrix or list of matrices (for multi-view nodes)
  exp <- object@Expectations[[variable]]
  # if (variable=="Z") {
  #   exp <- object@Expectations$Z
  # } else {
  #   exp <- lapply(object@Expectations[[variable]], function(x) x$E)
  # }
  
  # Convert to long data frame
  if (as.data.frame==T) {
    if (variable=="Z") {
      tmp <- reshape2::melt(exp)
      colnames(tmp) <- c("sample", "factor", "value")
      tmp[c("sample", "factor")] <- sapply(tmp[c("sample", "factor")], as.character)
    }
    else if (variable=="W") {
      tmp <- lapply(names(exp), function(m) { 
        tmp <- reshape2::melt(exp[[m]])
        colnames(tmp) <- c("feature", "factor", "value")
        tmp$view <- m
        tmp[c("view","feature","factor")] <- sapply(tmp[c("view", "feature", "factor")], as.character)
        return(tmp)
      })
      tmp <- do.call(rbind.data.frame,tmp)
    }
    else if (variable=="Y") {
      tmp <- lapply(names(exp), function(m) { 
        tmp <- reshape2::melt(exp[[m]])
        colnames(tmp) <- c("sample","feature","value")
        tmp$view <- m
        tmp[c("view","feature","factor")] <- sapply(tmp[c("view", "feature", "factor")], as.character)
        return(tmp) 
      })
      tmp <- do.call(rbind,tmp)
    }
    else if (variable=="Tau") {
      stop("Not implemented")
      # tmp <- lapply(names(exp), function(m) { 
      #   data.frame(view=m, feature=names(exp[[m]]), value=unname(exp[[m]]))
      #   tmp[c("view","feature","factor")] <- sapply(tmp[c("view","feature","factor")], as.character)
      #   return(tmp) 
      # })
      # tmp <- do.call(rbind,tmp)
    }
    else if (variable=="AlphaW" | variable=="AlphaZ") {
      tmp <- lapply(names(exp), function(m) { 
        tmp <- data.frame(view=m, factor=names(exp[[m]]), value=unname(exp[[m]]))
        tmp[c("view","feature","factor")] <- sapply(tmp[c("view","feature","factor")], as.character)
        return(tmp) 
      })
      tmp <- do.call(rbind,tmp)
    }
    else if (variable=="ThetaW" | variable=="ThetaZ") {
      stop("Not implemented")
      # tmp <- lapply(names(exp), function(m) { tmp <- reshape2::melt(exp[[m]]); colnames(tmp) <- c("sample","feature","value"); tmp$view <- m; tmp[c("view","feature","factor")] <- sapply(tmp[c("view","feature","factor")], as.character); return(tmp) })
      # tmp <- do.call(rbind,tmp)
    }
    else if (variable=="SigmaAlphaW" | variable=="SigmaZ") {
      stop("Not implemented")
    }
    exp <- tmp
  }
  return(exp)
}


