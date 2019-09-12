
################################################
## Get functions to fetch data from the model ##
################################################

#' @title Get dimensions
#' @name get_dimensions
#' @description Extract dimensionalities from the model. 
#' @details K indicates the number of factors, D indicates the number of features, 
#' N indicates the (total) number of samples and M indicates the number of views.
#' @param object a \code{\link{MOFA}} object.
#' @export
get_dimensions <- function(object) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  return(object@dimensions)
}

#' @title Get ELBO
#' @name get_elbo
#' @description Extract the value of the ELBO statistics after model training. This can be useful for model selection.
#' @details This can be useful for model selection.
#' @param object a \code{\link{MOFA}} object.
#' @export
get_elbo <- function(object) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  return(max(object@training_stats$elbo, na.rm=TRUE))
}

#' @title Get fractors
#' @name get_factors
#' @description Extract the latent factors from the model.
#' @param object a trained \code{\link{MOFA}} object.
#' @param factors character vector with the factor name(s), or numeric vector with the factor index(es).
#' Default is "all".
#' @param groups character vector with the group name(s), or numeric vector with the group index(es).
#' Default is "all".
#' @param as.data.frame logical indicating whether to return a long data frame instead of a matrix.
#' Default is \code{FALSE}.
#' @return By default it returns the latent factor matrix of dimensionality (N,K), where N is number of samples and K is number of factors. \cr
#' Alternatively, if \code{as.data.frame} is \code{TRUE}, returns a long-formatted data frame with columns (sample,factor,value).
#' @export
#' 
get_factors <- function(object, groups = "all", factors = "all", as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get factors and groups
  groups <- .check_and_get_groups(object, groups)
  factors <- .check_and_get_factors(object, factors)

  # Collect factors
  Z <- get_expectations(object, "Z", as.data.frame)
  if (as.data.frame) {
    Z <- Z[Z$factor%in%factors & Z$group%in%groups,]
  } else {
    Z <- lapply(Z[groups], function(z) z[,factors, drop=FALSE])
    names(Z) <- groups
  }

  return(Z)
}


#' @title Get weights
#' @name get_weights
#' @description Extract the weights from the model.
#' @param object a trained \code{\link{MOFA}} object.
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
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get views
  views <- .check_and_get_views(object, views)
  factors <- .check_and_get_factors(object, factors)
  
  # Fetch weights
  weights <- get_expectations(object, "W", as.data.frame)
  if (as.data.frame) {
    weights <- weights[weights$view %in% views & weights$factor %in% factors, ]
  } else {
    weights <- lapply(views, function(m) weights[[m]][,factors,drop=FALSE])
    names(weights) <- views
  }
  return(weights)
}


#' @title Get data
#' @name get_data
#' @description Fetch the input data
#' @param object a \code{\link{MOFA}} object.
#' @param views character vector with the view name(s), or numeric vector with the view index(es). 
#' Default is "all".
#' @param features a *named* list of character vectors. Example: list("view1"=c("feature_1","feature_2"), "view2"=c("feature_3","feature_4"))
#' Default is "all" If this is used, the argument views is ignored.
#' @param as.data.frame logical indicating whether to return a long data frame instead of a list of matrices.
#' @param add_intercept logical indicating whether to add feature intercepts to the data. Default is \code{TRUE}.
#' @details By default this function returns a list where each element is a data matrix with dimensionality (D,N) 
#' where D is the number of features and N is the number of samples. \cr
#' Alternatively, if \code{as.data.frame} is \code{TRUE}, the function returns a long-formatted data frame with columns (view,feature,sample,value).
#' @export
get_data <- function(object, views = "all", groups = "all", features = "all", as.data.frame = FALSE, add_intercept = TRUE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get views and groups
  views  <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)
  
  # Get features
  if (is(features, "list")) {
    if (is.null(names(features))) stop("features has to be a *named* list of character vectors. Please see the documentation")
    if (!all(sapply(names(features), function(i) all(features[[i]] %in% features(object)[[i]]) ))) stop("features not recognised")
    if (!all(sapply(features,length)>=1)) stop("features not recognised, please read the documentation")
    views <- names(features)
  } else {
    if (paste0(features, collapse="") == "all") { 
      features <- features(object)[views]
    } else {
      stop("features not recognised, please read the documentation")
    }
  }

  # Fetch data
  data <- lapply(object@data[views], function(x) x[groups])
  data <- lapply(seq_len(length(data)), function(m) lapply(seq_len(length(data[[1]])), function(p) data[[m]][[p]][as.character(features[[m]]),,drop=FALSE]))
  data <- .name_views_and_groups(data, views, groups)
  
  # Add feature intercepts (only for gaussian likelihoods)
  tryCatch( {
    
    if (add_intercept & length(object@intercepts[[1]])>0) {
      intercepts <- lapply(object@intercepts[views], function(x) x[groups]) 
      intercepts <- lapply(seq_len(length(intercepts)), function(m) lapply(seq_len(length(intercepts[[1]])), function(p) intercepts[[m]][[p]][as.character(features[[m]])]))
      intercepts <- .name_views_and_groups(intercepts, views, groups)
      
      for (m in names(data)) {
        if (object@model_options$likelihoods[[m]]=="gaussian") {
          for (g in names(data[[m]])) {
              data[[m]][[g]] <- data[[m]][[g]] + intercepts[[m]][[g]][as.character(features[[m]])]
          }
        }
      }
    } }, error = function(e) { NULL })

  # Convert to long data frame
  if (as.data.frame) {
    tmp <- lapply(views, function(m) { 
      lapply(groups, function(p) { 
        tmp <- reshape2::melt(data[[m]][[p]], na.rm=T)
        colnames(tmp) <- c("feature", "sample", "value")
        tmp <- cbind(view = m, group = p, tmp)
        return(tmp) 
      })
    })
    data <- do.call(rbind, do.call(rbind, tmp))
    factor.cols <- c("view","group","feature","sample")
    data[factor.cols] <- lapply(data[factor.cols], factor)
    
  }
  
  return(data)
}


#' @title Get imputed data
#' @name get_imputed_data
#' @description Function to get the imputed data. It requires the previous use of the \code{\link{impute}} method.
#' @param object a trained \code{\link{MOFA}} object.
#' @param views character vector with the view name(s), or numeric vector with the view index(es). 
#' Default is "all".
#' @param groups character vector with the group name(s), or numeric vector with the group index(es).
#' Default is "all".
#' @param features list of character vectors with the feature names or list of numeric vectors with the feature indices. 
#' Default is "all".
#' @param as.data.frame logical indicating whether to return a long-formatted data frame instead of a list of matrices. 
#' Default is \code{FALSE}.
#' @return By default returns a list where each element is a matrix with dimensionality (D,N), where D is the number of features in this view and N is the number of samples. \cr
#' Alternatively, if \code{as.data.frame} is \code{TRUE}, returns a long-formatted data frame with columns (view,feature,sample,value).
#' @export
get_imputed_data <- function(object, views = "all", groups = "all", features = "all", as.data.frame = FALSE, add_intercept = TRUE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (length(object@imputed_data)==0) stop("imputed data not found, did you run: 'object <- impute(object)'?")
  
  # Get views and groups
  views <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)
  
  # Get features
  if (is(features, "list")) {
    stopifnot(all(sapply(seq_len(length(features)), function(i) all(features[[i]] %in% features(object)[[views[i]]]))))
    stopifnot(length(features)==length(views))
    if (is.null(names(features))) names(features) <- views
  } else {
    if (paste0(features, collapse="") == "all") { 
      features <- features(object)[views]
    } else {
      stop("features not recognised, please read the documentation")
    }
  }
  
  # Fetch data
  # imputed_data <- lapply(object@imputed_data[views], function(x) x[groups])
  imputed_data <- lapply(object@imputed_data[views], function(x) lapply(x[groups],"[[","mean"))
  # imputed_data <- lapply(1:length(imputed_data), function(m) lapply(1:length(imputed_data[[1]]), function(p) imputed_data[[m]][[p]][as.character(features[[m]]),,drop=F]))
  imputed_data <- lapply(seq_len(length(imputed_data)), function(m) lapply(seq_len(length(imputed_data[[1]])), function(p) imputed_data[[m]][[p]][as.character(features[[m]]),,drop=FALSE]))
  imputed_data <- .name_views_and_groups(imputed_data, views, groups)
  
  # Add feature intercepts
  tryCatch( {
    
    if (add_intercept & length(object@intercepts[[1]])>0) {
      intercepts <- lapply(object@intercepts[views], function(x) x[groups]) 
      intercepts <- .name_views_and_groups(intercepts, views, groups)
      
      for (m in names(imputed_data)) {
        for (g in names(imputed_data[[m]])) {
          imputed_data[[m]][[g]] <- imputed_data[[m]][[g]] + intercepts[[m]][[g]][as.character(features[[m]])]
        }
      }
    } }, error = function(e) { NULL })
  
  # Convert to long data frame
  if (as.data.frame) {
    tmp <- lapply(views, function(m) { 
      lapply(groups, function(g) { 
        tmp <- reshape2::melt(imputed_data[[m]][[g]], na.rm=T)
        colnames(tmp) <- c("feature", "sample", "value")
        tmp <- cbind(view = m, group = g, tmp)
        return(tmp) 
      })
    })
    imputed_data <- do.call(rbind, do.call(rbind, tmp))
    factor.cols <- c("view","group","feature","sample")
    imputed_data[factor.cols] <- lapply(imputed_data[factor.cols], factor)
  }
  
  return(imputed_data)
}





#' @title Get expectations
#' @name get_expectations
#' @description Function to extract the expectations from the (variational) posterior distributions of a trained \code{\link{MOFA}} object.
#' @param object a trained \code{\link{MOFA}} object.
#' @param variable variable name: 'Z' for factors and 'W' for weights.
#' @param as.data.frame logical indicating whether to output the result as a long data frame, default is \code{FALSE}.
#' @details Technical note: MOFA is a Bayesian model where each variable has a prior distribution and a posterior distribution. 
#' In particular, to achieve scalability we used the variational inference framework, thus true posterior distributions are replaced by approximated variational distributions.
#' This function extracts the expectations of the variational distributions, which can be used as final point estimates to analyse the results of the model. \cr 
#' The priors and variational distributions of each variable are extensively described in the supplementary methods of the original paper.
#' @return the output varies depending on the variable of interest: \cr
#' \itemize{
#'  \item{"Z"}{a matrix with dimensions (samples,factors). If \code{as.data.frame} is \code{TRUE}, a long-formatted data frame with columns (sample,factor,value)}
#'  \item{"W"}{a list of length (views) where each element is a matrix with dimensions (features,factors). If \code{as.data.frame} is \code{TRUE}, a long-formatted data frame with columns (view,feature,factor,value)}
#' }
#' @export
get_expectations <- function(object, variable, as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(variable %in% names(object@expectations))
  
  # Get expectations in single matrix or list of matrices (for multi-view nodes)
  exp <- object@expectations[[variable]]

  # For memory and space efficiency, Y expectations are not saved to the model file when using only gaussian likelihoods.
  if (variable == "Y") {
    if ((length(object@expectations$Y) == 0) && all(object@model_options$likelihood == "gaussian")) {
      # message("Using training data slot as Y expectations since all the likelihoods are gaussian.")
      exp <- object@data
    }
  }
  
  # Convert to long data frame
  if (as.data.frame) {
    
    # Z node
    if (variable=="Z") {
      tmp <- reshape2::melt(exp, na.rm=T)
      colnames(tmp) <- c("sample", "factor", "value", "group")
      
      factor.cols <- c("sample", "factor", "group")
      factor.cols[factor.cols] <- lapply(factor.cols[factor.cols], factor)
    }
    
    # W node
    else if (variable=="W") {
      tmp <- lapply(names(exp), function(m) { 
        tmp <- reshape2::melt(exp[[m]], na.rm=T)
        colnames(tmp) <- c("feature","factor","value")
        tmp$view <- m
        factor.cols <- c("view", "feature", "factor")
        tmp[factor.cols] <- lapply(tmp[factor.cols], factor)
        return(tmp)
      })
      tmp <- do.call(rbind.data.frame,tmp)
    }
    
    # Y node
    else if (variable=="Y") {
      tmp <- lapply(names(exp), function(m) {
        tmp <- lapply(names(exp[[m]]), function(g) {
          tmp <- reshape2::melt(exp[[m]][[g]], na.rm=T)
          colnames(tmp) <- c("sample", "feature", "value")
          tmp$view <- m
          tmp$group <- g
          factor.cols <- c("view", "group", "feature", "factor")
          tmp[factor.cols] <- lapply(tmp[factor.cols], factor)
          return(tmp) 
        })
      })
      tmp <- do.call(rbind, tmp)
    }
    
    exp <- tmp
  }
  return(exp)
}

