
################################################
## Get functions to fetch data from the model ##
################################################

#' @title get_dimensions
#' @name get_dimensions
#' @description Extract dimensionalities from the model. 
#' @details K indicates the number of factors, D indicates the number of features, 
#' N indicates the (total) number of samples and M indicates the number of views.
#' @param object a \code{\link{BioFAModel}} object.
#' @export
get_dimensions <- function(object) {
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")  
  return(object@dimensions)
}


#' @title get_factors
#' @name get_factors
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
get_factors <- function(object, groups = "all", factors = "all", as.data.frame = FALSE, include_intercept = TRUE) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  
  # Get groups
  groups <- .check_and_get_groups(object, groups)

  # Get factors
  if (paste0(factors, collapse="") == "all") { factors <- factors_names(object) } 
  else if(is.numeric(factors)) {
    if (object@model_options$learn_intercept) factors <- factors_names(object)[factors+1]
    else factors <- factors_names(object)[factors]
  }
  else {
    stopifnot(all(factors %in% factors_names(object)))
  }
  
  # Collect factors
  Z <- get_expectations(object, "Z", as.data.frame)
  if (as.data.frame) {
    Z <- Z[Z$factor %in% factors,]
  } else {
    Z <- lapply(Z, function(z) z[,factors, drop=FALSE])
    names(Z) <- groups
  }

  # Remove intercept
  if (include_intercept == FALSE) {
    if (as.data.frame == FALSE) {
      Z <- lapply(names(Z), function(p) {
        Z[[p]][,colnames(Z[[p]])!="intercept"]
      })
      names(Z) <- groups
    } else {
      if ("intercept" %in% unique(Z$factor)) {
        Z <- Z[Z$factor!="intercept",]
      }
    }
  }
  return(Z)
}


#' @title get_weights
#' @name get_weights
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
get_weights <- function(object, views = "all", factors = "all", as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  
  # Get views
  views <- .check_and_get_views(object, views)
  
  # Get factors
  if (paste0(factors, collapse = "") == "all") { 
    factors <- factors_names(object)
  } else if (is.numeric(factors)) {
    if (object@model_options$learn_intercept) {
      factors <- factors_names(object)[factors+1]
    } else {
      factors <- factors_names(object)[factors]
    }
  } else { 
    stopifnot(all(factors %in% factors_names(object)))
  }
  
  # Fetch weights
  weights <- get_expectations(object, "W", as.data.frame)
  if (as.data.frame) {
    weights <- weights[weights$view %in% views & weights$factor %in% factors, ]
  } else {
    weights <- lapply(views, function(m) weights[[m]][,factors,drop=F])
    names(weights) <- views
    # if (length(views)==1) { weights <- weights[[1]] }
  }
  return(weights)
}


#' @title get_training_data
#' @name get_training_data
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
get_training_data <- function(object, views = "all", groups = "all", features = "all", as.data.frame = F) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  
  # Get views and groups
  views  <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)
  
  # Get features
  if (class(features) == "list") {
    stopifnot(all(sapply(1:length(features), function(i) all(features[[i]] %in% features_names(object)[[views[i]]]))))
  } else {
    if (paste0(features, collapse="") == "all") { 
      features <- features_names(object)[views]
    } else {
      stop("features not recognised, please read the documentation")
    }
  }

  # Fetch data
  train_data <- lapply(object@training_data[views], function(e) e[groups])
  train_data <- lapply(1:length(train_data), function(m) lapply(1:length(train_data[[1]]), function(p) train_data[[m]][[p]][features[[m]],,drop=F]))
  train_data <- .name_views_and_groups(train_data, views, groups)
  
  # Convert to long data frame
  if (as.data.frame==T) {
    tmp <- lapply(views, function(m) { 
      lapply(groups, function(p) { 
        tmp <- reshape2::melt(train_data[[m]][[p]])
        colnames(tmp) <- c("feature", "sample", "value")
        tmp <- cbind(view = m, group = p, tmp)
        return(tmp) 
      })
    })
    train_data <- do.call(rbind, do.call(rbind, tmp))
    train_data[,c("view","group","feature","sample")] <- sapply(train_data[,c("view","group","feature","sample")], as.character)
  }# else if ((length(views)==1) && (as.data.frame==F)) {
  #  train_data <- train_data[[views]]
  #}
  
  return(train_data)
}


#' @title getimputed_data
#' @name getimputed_data
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
getimputed_data <- function(object, views = "all", groups = "all", features = "all", as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  
  # Get views and groups
  if (paste0(views, collapse="") == "all") { views <- views_names(object) } else { stopifnot(all(views %in% views_names(object))) }
  if (paste0(groups, collapse="") == "all") { groups <- groups_names(object) } else { stopifnot(all(groups %in% groups_names(object))) }
  
  # Get features
  if (class(features)=="list") {
    stopifnot(all(sapply(1:length(features), function(i) all(features[[i]] %in% features_names(object)[[views[i]]]))))
  } else {
    if (paste0(features,collapse="") == "all") { 
      features <- features_names(object)[views]
    } else {
      stop("features not recognised, please read the documentation")
    }
  }
  

  # Fetch data
  imputed_data <- lapply(object@imputed_data[views], function(e) e[groups])
  imputed_data <- lapply(1:length(imputed_data), function(m) lapply(1:length(imputed_data[[1]]), function(h) imputed_data[[m]][[h]][features[[m]],,drop=F]))
  imputed_data <- .setViewAndgroups_names(imputed_data, views, groups)
  
  # Convert to long data frame
  if (as.data.frame==T) {
    tmp <- lapply(views, function(m) { 
      lapply(groups, function(h) { 
        tmp <- reshape2::melt(imputed_data[[m]][[h]])
        colnames(tmp) <- c("feature", "sample", "value")
        tmp <- cbind(view = m, group = h, tmp)
        return(tmp) 
      })
    })
    imputed_data <- do.call(rbind, tmp)
    imputed_data[,c("view","group","feature","sample")] <- sapply(imputed_data[,c("view","group","feature","sample")], as.character)
  } else if ((length(views)==1) && (as.data.frame==F)) {
    imputed_data <- imputed_data[[views]]
  }
  
  return(imputed_data)
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
  if(class(object@input_data) != "MultiAssayExperiment") stop("To work with covariates, InputData has to be specified in form of a MultiAssayExperiment")  
  stopifnot(all(covariates %in% colnames(colData(object@input_data))))
  
  # Get covariates
  covariates <- colData(object@input_data)[,covariates]
  
  return(covariates)
}

#' @title get_expectations
#' @name get_expectations
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
get_expectations <- function(object, variable, as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  stopifnot(variable %in% names(object@expectations))
  
  # Get expectations in single matrix or list of matrices (for multi-view nodes)
  exp <- object@expectations[[variable]]
  # if (variable=="Z") {
  #   exp <- object@expectations$Z
  # } else {
  #   exp <- lapply(object@expectations[[variable]], function(x) x$E)
  # }
  
  # Convert to long data frame
  if (as.data.frame) {
    if (variable=="Z") {
      tmp <- reshape2::melt(exp)
      colnames(tmp) <- c("sample", "factor", "value", "group")
      tmp[c("sample", "factor", "group")] <- sapply(tmp[c("sample", "factor", "group")], as.character)
    }
    else if (variable=="W") {
      tmp <- lapply(names(exp), function(m) { 
        tmp <- reshape2::melt(exp[[m]])
        colnames(tmp) <- c("feature", "factor", "value")
        tmp$view <- m
        tmp[c("view", "feature", "factor")] <- sapply(tmp[c("view", "feature", "factor")], as.character)
        return(tmp)
      })
      tmp <- do.call(rbind.data.frame,tmp)
    }
    else if (variable=="Y") {
      tmp <- lapply(names(exp), function(m) {
        tmp <- lapply(names(exp[[m]]), function(h) {
          tmp <- reshape2::melt(exp[[m]][[h]])
          colnames(tmp) <- c("sample", "feature", "value")
          tmp$view <- m
          tmp$group <- h
          tmp[c("view", "group", "feature", "factor")] <- sapply(tmp[c("view", "group", "feature", "factor")], as.character)
          return(tmp) 
        })
      })
      tmp <- do.call(rbind, tmp)
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


#' @title get_elbo
#' @name get_elbo
#' @description Extract the value of the ELBO statistics after model training. This can be useful for model selection.
#' @param object a \code{\link{BioFAModel}} object.
#' @export
get_elbo <- function(object) {
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")  
  return(tail(object@TrainStats$elbo, 1))
}

get_groups_annotation <- function(object){
  samples_list <- samples_names(object)
  if(class(samples_list) == "list") samples <- Reduce(c, samples_list) else samples <- samples_list
  data.frame(sample = samples, group = rep(names(samples_list), times = sapply(samples_list, length)))
}
