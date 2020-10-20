
################################################
## Get functions to fetch data from the model ##
################################################

#' @title Get dimensions
#' @name get_dimensions
#' @description Extract dimensionalities from the model. 
#' @details K indicates the number of factors, D indicates the number of features, 
#' N indicates the (total) number of samples, M indicates the number of views and C indicates the number of covariates.
#' @param object a \code{\link{MOFA}} object.
#' @return list containing the dimensionalities of the model
#' @export
#' @examples
#' # Using an existing trained model
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' dims <- get_dimensions(model)

get_dimensions <- function(object) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  return(object@dimensions)
}

#' @title Get ELBO
#' @name get_elbo
#' @description Extract the value of the ELBO statistics after model training. This can be useful for model selection.
#' @details This can be useful for model selection.
#' @param object a \code{\link{MOFA}} object.
#' @return Value of the ELBO
#' @export
#' @examples
#' # Using an existing trained model
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' elbo <- get_elbo(model)

get_elbo <- function(object) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  return(max(object@training_stats$elbo, na.rm=TRUE))
}

#' @title Get lengthscales
#' @name get_lengthscales
#' @description Extract the inferred lengthscale for each factor after model training. 
#' @details This can be used only if covariates are passed to the MOFAobject upon creation and GP_factors is set to True.
#' @param object a \code{\link{MOFA}} object.
#' @export
get_lengthscales <- function(object) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if(is.null(object@covariates)) stop("No covariates specified in 'object'")
  if(is.null(object@training_stats$length_scales)) stop("No lenghtscales saved in 'object' \n Make sure you specify the covariates and train setting the option 'GP_factors' to TRUE.")
  tmp <- object@training_stats$length_scales
  return(tmp)
}


#' @title Get scales
#' @name get_scales
#' @description Extract the inferred scale for each factor after model training. 
#' @details This can be used only if covariates are passed to the MOFAobject upon creation and GP_factors is set to True.
#' @param object a \code{\link{MOFA}} object.
#' @export
get_scales <- function(object) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if(is.null(object@covariates)) stop("No covariates specified in 'object'")
  if(is.null(object@training_stats$scales)) stop("No scales saved in 'object' \n Make sure you specify the covariates and train setting the option 'GP_factors' to TRUE.")
  tmp <- object@training_stats$scales
  return(tmp)
}

#' @title Get group covariance matrix
#' @name get_group_kernel
#' @description Extract the inferred group-group covariance matrix per factor
#' @details This can be used only if covariates are passed to the MOFAobject upon creation and GP_factors is set to True.
#' @param object a \code{\link{MOFA}} object.
#' @export
get_group_kernel <- function(object) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if(is.null(object@covariates)) stop("No covariates specified in 'object'")
  if(object@dimensions$G < 2) stop("No groups specified in 'object'")
  if(is.null(object@training_stats$Kg)) stop("No group kernel saved in 'object' \n Make sure you specify the covariates and train setting the option 'GP_factors' as well as 'model_groups' to TRUE.")
  tmp <- lapply(seq_len(dim(object@training_stats$Kg)[3]), function(x) {
    mat <- object@training_stats$Kg[ , , x]
    rownames(mat) <- colnames(mat) <- groups_names(object)
    mat
    })
  names(tmp) <- factors_names(object)
  return(tmp)
}

#' @title Get interpolated factor values
#' @name get_interpolated_factors
#' @description Extract the interpolated factor values
#' @details This can be used only if covariates are passed to the object upon creation, GP_factors is set to True and new covariates were passed for interpolation.
#' @param object a \code{\link{MOFA}} object
#' @param as.data.frame logical indicating whether to return data as a data.frame
#' @export
get_interpolated_factors <- function(object, as.data.frame = FALSE) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if(is.null(object@interpolated_Z)) stop("No interpolated factors present in 'object'")
  if(!as.data.frame){
    return(object@interpolated_Z)
  } else {
    preds <- lapply(object@interpolated_Z, function(l) l[names(l)[names(l) != "new_values"]])
    df_interpol <- melt(preds, varnames = c("factor", "sample_id"))
    df_interpol <- rename(df_interpol, group = L1, type = L2)
    
    if("new_values" %in% names(object@interpolated_Z[[1]])) {
      new_vals <- lapply(object@interpolated_Z, function(l) l[names(l)[names(l) == "new_values"]])
      new_vals <- melt(new_vals, varnames = c("covariate","sample_id"))
      new_vals <- mutate(new_vals, covariate = covariates_names(object)[covariate])
      new_vals <- rename(new_vals, group = L1, covariate_value = value)
      new_vals <- spread(new_vals, key = covariate, value = covariate_value)
      new_vals <- select(new_vals, -L2)
      df_interpol <- left_join(df_interpol, new_vals, by = c("group", "sample_id"))
      df_interpol <- select(df_interpol, -sample_id)
    } else { # compatibility to older objects
      df_interpol <- rename(df_interpol, covariate_value = sample_id)
      df_interpol <- mutate(df_interpol, covariate = covariates_names(object))
    }
    df_interpol <- mutate(df_interpol, factor = factors_names(object)[factor])
    df_interpol <- spread(df_interpol, key = type, value = value)
    return(df_interpol)
  }
}


#' @title Get factors
#' @name get_factors
#' @description Extract the latent factors from the model.
#' @param object a trained \code{\link{MOFA}} object.
#' @param factors character vector with the factor name(s), or numeric vector with the factor index(es).
#' Default is "all".
#' @param groups character vector with the group name(s), or numeric vector with the group index(es).
#' Default is "all".
#' @param scale logical indicating whether to scale factor values.
#' @param as.data.frame logical indicating whether to return a long data frame instead of a matrix.
#' Default is \code{FALSE}.
#' @return By default it returns the latent factor matrix of dimensionality (N,K), where N is number of samples and K is number of factors. \cr
#' Alternatively, if \code{as.data.frame} is \code{TRUE}, returns a long-formatted data frame with columns (sample,factor,value).
#' @export
#' 
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#'
#' # Fetch factors in matrix format (a list, one matrix per group)
#' factors <- get_factors(model)
#'
#' # Concatenate groups
#' factors <- do.call("rbind",factors)
#'
#' # Fetch factors in data.frame format instead of matrix format
#' factors <- get_factors(model, as.data.frame = TRUE)
get_factors <- function(object, groups = "all", factors = "all", scale = FALSE, as.data.frame = FALSE) {

  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get groups
  groups <- .check_and_get_groups(object, groups)
  # Get factors
  factors <- .check_and_get_factors(object, factors)
  
  # Collect factors
  Z <- get_expectations(object, "Z", as.data.frame)
  if (as.data.frame) {
    Z <- Z[Z$factor%in%factors & Z$group%in%groups,]
    if (scale) Z$value <- Z$value/max(abs(Z$value),na.rm=TRUE)
  } else {
    Z <- lapply(Z[groups], function(z) z[,factors, drop=FALSE])
    if (scale) Z <- lapply(Z, function(x) x/max(abs(x)) )
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
#' @param abs logical indicating whether to take the absolute value of the weights.
#' @param scale logical indicating whether to scale all weights from -1 to 1 (or from 0 to 1 if \code{abs=TRUE}).
#' @param as.data.frame logical indicating whether to return a long data frame instead of a list of matrices. 
#' Default is \code{FALSE}.
#' @return By default it returns a list where each element is a loading matrix with dimensionality (D,K), 
#' where D is the number of features and K is the number of factors. \cr
#' Alternatively, if \code{as.data.frame} is \code{TRUE}, returns a long-formatted data frame with columns (view,feature,factor,value).
#' @export
#' 
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#'
#' # Fetch weights in matrix format (a list, one matrix per view)
#' weights <- get_weights(model)
#'
#' # Fetch weights for factor 1 and 2 and view 1
#' weights <- get_weights(model, views = 1, factors = c(1,2))
#'
#' # Fetch weights in data.frame format
#' weights <- get_weights(model, as.data.frame = TRUE)

get_weights <- function(object, views = "all", factors = "all", abs = FALSE, scale = FALSE, as.data.frame = FALSE) {

  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get views
  views <- .check_and_get_views(object, views)
  factors <- .check_and_get_factors(object, factors)
  
  # Fetch weights
  weights <- get_expectations(object, "W", as.data.frame)
  
  if (as.data.frame) {
    weights <- weights[weights$view %in% views & weights$factor %in% factors, ]
    if (abs) weights$value <- abs(weights$value)
    if (scale) weights$value <- weights$value/max(abs(weights$value))
  } else {
    weights <- lapply(weights[views], function(x) x[,factors,drop=FALSE])
    if (abs) weights <- lapply(weights, abs)
    if (scale) weights <- lapply(weights, function(x) x/max(abs(x)) )
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
#' @param groups character vector with the group name(s), or numeric vector with the group index(es). 
#' Default is "all".
#' @param features a *named* list of character vectors. Example: list("view1"=c("feature_1","feature_2"), "view2"=c("feature_3","feature_4"))
#' Default is "all".
#' @param as.data.frame logical indicating whether to return a long data frame instead of a list of matrices. Default is \code{FALSE}.
#' @param add_intercept logical indicating whether to add feature intercepts to the data. Default is \code{TRUE}.
#' @param denoise logical indicating whether to return the denoised data (i.e. the model predictions). Default is \code{FALSE}.
#' @param na.rm remove NAs from the data.frame (only if as.data.frame is \code{TRUE}).
#' @details By default this function returns a list where each element is a data matrix with dimensionality (D,N) 
#' where D is the number of features and N is the number of samples. \cr
#' Alternatively, if \code{as.data.frame} is \code{TRUE}, the function returns a long-formatted data frame with columns (view,feature,sample,value).
#' Missing values are not included in the the long data.frame format by default. To include them use the argument \code{na.rm=FALSE}.
#' @return A  list of data matrices with dimensionality (D,N) or a \code{data.frame} (if \code{as.data.frame} is TRUE)
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#'
#' # Fetch data
#' data <- get_data(model)
#'
#' # Fetch a specific view
#' data <- get_data(model, views = "view_0")
#'
#' # Fetch data in data.frame format instead of matrix format
#' data <- get_data(model, as.data.frame = TRUE)
#'
#' # Fetch centered data (do not add the feature intercepts)
#' data <- get_data(model, as.data.frame = FALSE)
#' 
#' # Fetch denoised data (do not add the feature intercepts)
#' data <- get_data(model, denoise = TRUE)
get_data <- function(object, views = "all", groups = "all", features = "all", as.data.frame = FALSE, add_intercept = TRUE, denoise = FALSE, na.rm = TRUE) {

  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get views and groups
  views  <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)
  
  # Get features
  if (is(features, "list")) {
    if (is.null(names(features))) stop("features has to be a *named* list of character vectors. Please see the documentation")
    if (!(names(features)%in%views_names(object))) stop("Views not recognised")
    if (!all(sapply(names(features), function(i) all(features[[i]] %in% features_names(object)[[i]]) ))) stop("features not recognised")
    if (any(sapply(features,length)<1)) stop("features not recognised, please read the documentation")
    views <- names(features)
  } else {
    if (paste0(features, collapse="") == "all") { 
      features <- features_names(object)[views]
    } else {
      stop("features not recognised, please read the documentation")
    }
  }

  # Fetch data
  if (denoise) {
    data <- predict(object, views=views, groups=groups)
  } else {
    data <- lapply(object@data[views], function(x) x[groups])
  }
  data <- lapply(views, function(m) lapply(seq_len(length(data[[1]])), function(p) data[[m]][[p]][as.character(features[[m]]),,drop=FALSE]))
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
        tmp <- reshape2::melt(data[[m]][[p]], na.rm=na.rm)
        if(nrow(tmp) >0 & !is.null(tmp)) {
        colnames(tmp) <- c("feature", "sample", "value")
        tmp <- cbind(view = m, group = p, tmp)
        return(tmp) 
        } 
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
#' @details Data is imputed from the generative model of MOFA.
#' @return A list containing the imputed valued or a data.frame if as.data.frame is TRUE
#' @export
#' @examples
#' # Using an existing trained model
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' model <- impute(model)
#' imputed <- get_imputed_data(model)

get_imputed_data <- function(object, views = "all", groups = "all", features = "all", as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (length(object@imputed_data)==0) stop("imputed data not found, did you run: 'object <- impute(object)'?")
  
  # Get views and groups
  views <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)

  # Get features
  if (is(features, "list")) {
    stopifnot(all(sapply(seq_len(length(features)), function(i) all(features[[i]] %in% features_names(object)[[views[i]]]))))
    stopifnot(length(features)==length(views))
    if (is.null(names(features))) names(features) <- views
  } else {
    if (paste0(features, collapse="") == "all") { 
      features <- features_names(object)[views]
    } else {
      stop("features not recognised, please read the documentation")
    }
  }
  
  # Fetch mean
  imputed_data <- lapply(object@imputed_data[views], function(x) x[groups] )
  imputed_data <- lapply(seq_len(length(imputed_data)), function(m) lapply(seq_len(length(imputed_data[[1]])), function(p) imputed_data[[m]][[p]][as.character(features[[m]]),,drop=FALSE]))
  imputed_data <- .name_views_and_groups(imputed_data, views, groups)
  
# Add feature intercepts
# tryCatch( {
#
#   if (add_intercept & length(object@intercepts[[1]])>0) {
#     intercepts <- lapply(object@intercepts[views], function(x) x[groups])
#     intercepts <- .name_views_and_groups(intercepts, views, groups)
#
#     for (m in names(imputed_data)) {
#       for (g in names(imputed_data[[m]])) {
#         imputed_data[[m]][[g]] <- imputed_data[[m]][[g]] + intercepts[[m]][[g]][as.character(features[[m]])]
#       }
#     }
#   } }, error = function(e) { NULL })

  # Convert to long data frame
  if (isTRUE(as.data.frame)) {
    
    imputed_data <- lapply(views, function(m) { 
      lapply(groups, function(g) { 
        tmp <- reshape2::melt(imputed_data[[m]][[g]])
        colnames(tmp) <- c("feature", "sample", "value")
        tmp <- cbind(view = m, group = g, tmp)
        return(tmp) 
      })
    })
    imputed_data <- do.call(rbind, do.call(rbind, imputed_data))
    

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
#' @examples
#' # Using an existing trained model
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' factors <- get_expectations(model, "Z")
#' weights <- get_expectations(model, "W")

get_expectations <- function(object, variable, as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(variable %in% names(object@expectations))
  
  # Get expectations in single matrix or list of matrices (for multi-view nodes)
  exp <- object@expectations[[variable]]

  # unlist single view nodes - single Sigma node across all groups using time warping
  if(variable == "Sigma")
    exp <- exp[[1]]
  
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
      tmp <- reshape2::melt(exp, na.rm=TRUE)
      colnames(tmp) <- c("sample", "factor", "value", "group")
      tmp$sample <- as.character(tmp$sample)
      factor.cols <- c("sample", "factor", "group")
      factor.cols[factor.cols] <- lapply(factor.cols[factor.cols], factor)
    }
    
    # W node
    else if (variable=="W") {
      tmp <- lapply(names(exp), function(m) { 
        tmp <- reshape2::melt(exp[[m]], na.rm=TRUE)
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
          tmp <- reshape2::melt(exp[[m]][[g]], na.rm=TRUE)
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


#' @title Get variance explained values
#' @name get_variance_explained
#' @description Extract the latent factors from the model.
#' @param object a trained \code{\link{MOFA}} object.
#' @param factors character vector with the factor name(s), or numeric vector with the factor index(es).
#' Default is "all".
#' @param groups character vector with the group name(s), or numeric vector with the group index(es).
#' Default is "all".
#' @param views character vector with the view name(s), or numeric vector with the view index(es).
#' Default is "all".
#' @param as.data.frame logical indicating whether to return a long data frame instead of a matrix.
#' Default is \code{FALSE}.
#' @return A list of data matrices with variance explained per group or a \code{data.frame} (if \code{as.data.frame} is TRUE)
#' @export
#'
#' @examples
#' # Using an existing trained model
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#'
#' # Fetch variance explained values (in matrix format)
#' r2 <- get_variance_explained(model)
#'
#' # Fetch variance explained values (in data.frame format)
#' r2 <- get_variance_explained(model, as.data.frame = TRUE)
#'
get_variance_explained <- function(object, groups = "all", views = "all", factors = "all", 
                                   as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get factors and groups
  groups <- .check_and_get_groups(object, groups)
  views <- .check_and_get_views(object, views)
  factors <- .check_and_get_factors(object, factors)
  
  # Fetch R2
  if (.hasSlot(object, "cache") && ("variance_explained" %in% names(object@cache))) {
    r2_list <- object@cache$variance_explained
  } else {
    r2_list <- calculate_variance_explained(object, factors = factors, views = views, groups = groups)
  }
  
  # Convert to data.frame format
  if (as.data.frame) {
    
    # total R2
    r2_total <- reshape2::melt( do.call("rbind",r2_list[["r2_total"]]) )
    colnames(r2_total) <- c("group", "view", "value")
    
    # R2 per factor
    r2_per_factor <- lapply(names(r2_list[["r2_per_factor"]]), function(g) {
      x <- reshape2::melt( r2_list[["r2_per_factor"]][[g]] )
      colnames(x) <- c("factor", "view", "value")
      x$factor <- as.factor(x$factor)
      x$group <- g
      return(x)
    })
    r2_per_factor <- do.call("rbind",r2_per_factor)[,c("group","view","factor","value")]
    r2 <- list("r2_per_factor"=r2_per_factor, "r2_total"=r2_total)
    
  } else {
    r2 <- r2_list
  }
  
  return(r2)
  
}