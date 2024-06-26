
# Function to find "intercept" factors
# .detectInterceptFactors <- function(object, cor_threshold = 0.75) {
#   
#   # Sanity checks
#   if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
#   
#   # Fetch data
#   data <- getTrainData(object)
#   factors <- getfactors_names(object)
#   
#   # Correlate the factors with global means per sample
#   r <- lapply(data, function(x) abs(cor(colSums(x,na.rm=T),factors, use="pairwise.complete.obs")))
#   
#   token <- 0
#   for (i in names(r)) {
#     if (any(r[[i]]>cor_threshold)) {
#       token <- 1
#       message(paste0("Warning: Factor ",which(r[[i]]>cor_threshold)," is strongly correlated with the total expression for each sample in ",i))
#     }
#   }
#   if (token==1)
#     message("Such (strong) factors usually appear when count-based assays are not properly normalised by library size.")
# }

# x: a named vector, where names correspond to sample names
.add_column_to_metadata <- function(object, x, name) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(!is.null(names(x)))
  stopifnot(names(x) %in% unlist(samples_names(object)))
  
  # sort vector to match samples names (fill with NA where applicable)
  vec <- rep(NA,sum(get_dimensions(object)[["N"]]))
  names(vec) <- object@samples_metadata$sample
  vec[names(x)] <- x
  
  # add to metadata
  object@samples_metadata[[name]] <- x
  
  return(object)  
}

.infer_likelihoods <- function(object) {
  
  # Gaussian by default
  likelihood <- rep(x="gaussian", times=object@dimensions$M)
  names(likelihood) <- views_names(object)
  
  for (m in views_names(object)) {
    # data <- get_data(object, views=m)[[1]][[1]]  # take only first group
    data <- object@data[[m]][[1]]
    
    # bernoulli
    if (length(unique(data[!is.na(data)]))==2) {
      likelihood[m] <- "bernoulli"
    # poisson
    } else if (all(data[!is.na(data)]%%1==0)) {
      likelihood[m] <- "poisson"
    }
  }
  
  return(likelihood)
}

# Set view names and group names for nested list objects (e.g. Y)
.name_views_and_groups <- function(nested_list, view_names, group_names) {
  names(nested_list) <- view_names
  for (view in view_names) { names(nested_list[[view]]) <- group_names }
  nested_list
}

#' @importFrom stats sd
.detect_outliers <- function(object, groups = "all", factors = "all") {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Define groups
  groups <- .check_and_get_groups(object, groups)
  H <- length(groups)
  
  # Define factors
  factors <- .check_and_get_factors(object, factors)
  
  for (k in factors) {
    for (g in groups) {
    
      Z <- get_factors(object, groups=g, factors=k)[[1]][,1]
      Z <- Z[!is.na(Z)]
      
      # warning("Outlier detection is independent of the inferred lengthscale currently - might lead to unwanted results")
      cutoff <- 2.5 * 1.96
      tmp <- abs(Z - mean(Z)) / sd(Z)

      outliers <- names(which(tmp>cutoff & abs(Z)>0.5))
      
      if (length(outliers)>0 & length(outliers)<5) {
        object@expectations$Z[[g]][,k][outliers] <- NA
      }
      
    }
  }
  
  # re-compute variance explained
  object@cache[["variance_explained"]] <- calculate_variance_explained(object)
  
  return(object)
}


.flip_factor <- function(model, factor){
  for(g in names(model@expectations$Z)) {
    model@expectations$Z[[g]][,factor] <- - model@expectations$Z[[g]][,factor]
  }
  for(m in names(model@expectations$W)) {
    model@expectations$W[[m]][,factor] <- -model@expectations$W[[m]][,factor]
  }
return(model)
}


.check_and_get_factors <- function(object, factors) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(!any(duplicated(factors)))
  if (is.numeric(factors)) {
    stopifnot(all(factors <= object@dimensions$K))
    factors_names(object)[factors] 
  } else {
    if (paste0(factors, collapse = "") == "all") { 
      factors_names(object)
    } else {
      stopifnot(all(factors %in% factors_names(object)))
      factors
    }
  }
  
}

.check_and_get_covariates <- function(object, covariates) {
  if (!.hasSlot(object, "covariates") || is.null(object@covariates))
    stop("No covariates found in object.")
  stopifnot(!any(duplicated(covariates)))
  if (is.numeric(covariates)) {
    stopifnot(all(covariates <= object@dimensions$C))
    covariates_names(object)[covariates] 
  } else {
    if (paste0(covariates, collapse = "") == "all") { 
      covariates_names(object)
    } else {
      stopifnot(all(covariates %in% covariates_names(object)))
      covariates
    }
  }
}

.check_and_get_views <- function(object, views, non_gaussian=TRUE) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(!any(duplicated(views)))
  if (is.numeric(views)) {
    stopifnot(all(views <= object@dimensions$M))
    views <- views_names(object)[views]
  } else {
    if (paste0(views, sep = "", collapse = "") == "all") { 
      views <- views_names(object)
    } else {
      stopifnot(all(views %in% views_names(object)))
    }
  }
  
  # Ignore non-gaussian views  
  if (isFALSE(non_gaussian)) {
    non_gaussian_views <- names(which(object@model_options$likelihoods!="gaussian"))
    views <- views[!views%in%non_gaussian_views]
  }
  
  return(views)
}


.check_and_get_groups <- function(object, groups) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(!any(duplicated(groups)))
  if (is.numeric(groups)) {
    stopifnot(all(groups <= object@dimensions$G))
    groups_names(object)[groups] 
  } else {
    if (paste0(groups, collapse = "") == "all") { 
      groups_names(object)
    } else {
      stopifnot(all(groups %in% groups_names(object)))
      groups
    }
  }
}


.check_and_get_samples <- function(object, samples) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(!any(duplicated(samples)))
  if (is.numeric(samples)) {
    stopifnot(all(samples <= sum(object@dimensions$N)))
    unlist(samples_names(object))[samples] 
  } else {
    if (paste0(samples, collapse = "") == "all") { 
      unlist(samples_names(object))
    } else {
      stopifnot(all(samples %in% unlist(samples_names(object))))
      samples
    }
  }
}

.check_and_get_features_from_view <- function(object, view, features) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(!any(duplicated(features)))
  if (is.numeric(features)) {
    stopifnot(all(features <= sum(object@dimensions$D[view])))
    unname(unlist(features_names(object)[[view]])[features])
  } else {
    if (paste0(features, collapse = "") == "all") { 
      unlist(features_names(object)[[view]])
    } else {
      stopifnot(all(features %in% unlist(features_names(object)[[view]])))
      features
    }
  }
}

.get_top_features_by_loading <- function(object, view, factors, nfeatures = 10) {
  
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Collect expectations  
  W <- get_weights(object, factors = factors, views = view, as.data.frame=TRUE)
  # Work with absolute values to sort them
  W$value <- abs(W$value)

  # Extract relevant features
  W <- W[with(W, order(-abs(value))), ]

  return(as.character(head(W$feature, nfeatures)))
}


.get_nodes_types <- function() {
  nodes_types <- list(
    multiview_nodes  = c("W", "AlphaW", "ThetaW"),
    multigroup_nodes = c("Z", "AlphaZ", "ThetaZ"),
    twodim_nodes     = c("Y", "Tau"),
    multivariate_singleview_node = "Sigma"
  )
}

setClass("matrix_placeholder", 
         slots=c(rownames = "ANY",
                 colnames = "ANY",
                 nrow     = "integer",
                 ncol     = "integer")
)

setMethod("rownames", "matrix_placeholder", function(x) { x@rownames })
setMethod("colnames", "matrix_placeholder", function(x) { x@colnames })
setMethod("nrow", "matrix_placeholder", function(x) { x@nrow })
setMethod("ncol", "matrix_placeholder", function(x) { x@ncol })

setReplaceMethod("rownames", signature(x = "matrix_placeholder"),
  function(x, value) { 
    x@rownames <- value 
    x@nrow <- length(value)
    x 
    })
setReplaceMethod("colnames", signature(x = "matrix_placeholder"),
  function(x, value) { 
    x@colnames <- value 
    x@ncol <- length(value)
    x 
    })

.create_matrix_placeholder <- function(rownames, colnames) {
  mx <- new("matrix_placeholder")
  mx@rownames <- rownames
  mx@colnames <- colnames
  mx@nrow <- length(rownames)
  mx@ncol <- length(colnames)
  mx
}




# (Hidden) function to define the group
.set_groupby <- function(object, group_by) {
  
  # Option 0: no group
  if (is.null(group_by)) {
    group_by <- rep("1",sum(object@dimensions[["N"]]))
    
    # Option 1: by default group
  } else if (group_by[1] == "group") {
    # group_by = c()
    # for (group in names(samples_names(object))) {
    #   group_by <- c(group_by,rep(group,length(samples_names(object)[[group]])))
    # }
    # group_by = factor(group_by, levels=groups_names(object))
    group_by <- samples_metadata(object)$group
    
    # Option 2: by a metadata column in object@samples$metadata
  } else if ((length(group_by) == 1) && (is.character(group_by)|is.factor(group_by)) & (group_by[1] %in% colnames(samples_metadata(object)))) {
    group_by <- samples_metadata(object)[,group_by]
    # if (is.character(group_by)) group_by <- as.factor( group_by )
    
    # Option 3: input is a data.frame with columns (sample,group)
  } else if (is(group_by,"data.frame")) {
    stopifnot(all(colnames(group_by) %in% c("sample","group")))
    stopifnot(all(unique(group_by$sample) %in% unlist(samples_names(object))))
    
    # Option 4: group_by is a vector of length N
  } else if (length(group_by) > 1) {
    stopifnot(length(group_by) == sum(object@dimensions[["N"]]))
    
    # Option not recognised
  } else {
    stop("'group_by' was specified but it was not recognised, please read the documentation")
  }
  
  # Create data.frame with columns (sample,group)
  if (!is(group_by,"data.frame")) {
    df = data.frame(
      # sample = unlist(samples_names(object)),
      sample = samples_metadata(object)$sample,
      group_by = group_by,
      stringsAsFactors = FALSE
    )
    
  }
  
  return(df)
}

# (Hidden) function to define the color
.set_xax <- function(object, xax) {
  
    # Option 1: by a metadata column in object@samples_metadata
  if ((length(xax) == 1) && (is.character(xax)|is.factor(xax)) & (xax[1] %in% colnames(samples_metadata(object)))) {
    xax <- samples_metadata(object)[,xax]
    
    # Option 2: by a feature present in the training data    
  } else if ((length(xax) == 1) && is.character(xax) && (xax[1] %in% unlist(features_names(object)))) {
    data <- lapply(get_data(object), function(l) Reduce(cbind, l))
    features <- lapply(data, rownames)
    viewidx <- which(sapply(features, function(x) xax %in% x))
    xax <- data[[viewidx]][xax,]
    
    # Option 5: input is a data.frame with columns (sample, value)
  } else if (is(xax, "data.frame")) {
    stopifnot(all(colnames(xax) %in% c("sample", "value")))
    stopifnot(all(unique(xax$sample) %in% unlist(samples_names(object))))
    xax <- dplyr::rename(xax, covariate_value = value)
    # Option 6: color_by is a vector of length N
  } else if (length(xax) > 1) {
    stopifnot(length(xax) == sum(get_dimensions(object)$N))
    
    # Option not recognised
  } else {
    stop("'xax' was specified but it was not recognised, please read the documentation")
  }
  
  # Create data.frame with columns (sample,color)
  if (!is(xax,"data.frame")) {
    xax = data.frame(
      sample = unlist(samples_names(object)),
      covariate_value = xax,
      stringsAsFactors = FALSE
    )
  }
  return(xax)
}

# (Hidden) function to define the color
.set_colorby <- function(object, color_by) {
  
  # Option 0: no color
  if (is.null(color_by)) {
    color_by <- rep("1",sum(object@dimensions[["N"]]))
    
    # Option 1: by default group
  } else if (color_by[1] == "group") {
    color_by <- samples_metadata(object)$group
    
    # Option 2: by a metadata column in object@samples$metadata
  } else if ((length(color_by) == 1) && (is.character(color_by)|is.factor(color_by)) && (color_by[1] %in% colnames(samples_metadata(object)))) {
    color_by <- samples_metadata(object)[,color_by]
    # if (is.character(color_by)) color_by <- as.factor( color_by )
    
    # Option 3: by a feature present in the training data    
  } else if ((length(color_by) == 1) && is.character(color_by) && (color_by[1] %in% unlist(features_names(object)))) {
    viewidx <- which(sapply(features_names(object), function(x) color_by %in% x))
    foo <- list(color_by); names(foo) <- names(viewidx)
    color_by <- lapply(get_data(object, features = foo), function(l) Reduce(cbind, l))[[1]][1,]
    
    # data <- lapply(get_data(object), function(l) Reduce(cbind, l))
    # features <- lapply(data, rownames)
    # viewidx <- which(sapply(features, function(x) color_by %in% x))
    # color_by <- data[[viewidx]][color_by,]
    
    
    # Option 4: by a factor value in object@expectations$Z
  } else if ((length(color_by) == 1) && is.character(color_by) && (color_by[1] %in% colnames(get_factors(object)[[1]]))) {
    color_by <- do.call(rbind, get_factors(object))[,color_by]
    
    # Option 5: input is a data.frame with columns (sample, color)
  } else if (is(color_by, "data.frame")) {
    stopifnot(all(colnames(color_by) %in% c("sample", "color")))
    stopifnot(all(unique(color_by$sample) %in% unlist(samples_names(object))))
    
    # Option 6: color_by is a vector of length N
  } else if (length(color_by) > 1) {
    stopifnot(length(color_by) == sum(get_dimensions(object)$N))
    
    # Option not recognised
  } else {
    stop("'color_by' was specified but it was not recognised, please read the documentation")
  }
  
  # Create data.frame with columns (sample,color)
  if (!is(color_by,"data.frame")) {
    df <- data.frame(
      # sample = unlist(samples_names(object)),
      sample = samples_metadata(object)$sample,
      color_by = color_by,
      stringsAsFactors = FALSE
    )
  }
  if (length(unique(df$color_by)) < 5) df$color_by <- as.factor(df$color_by)
  
  return(df)
}


# (Hidden) function to define the shape
.set_shapeby <- function(object, shape_by) {
  
  # Option 0: no color
  if (is.null(shape_by)) {
    shape_by <- rep("1",sum(object@dimensions[["N"]]))
    
    # Option 1: by default group
  } else if (shape_by[1] == "group") {
    # shape_by = c()
    # for (group in names(samples_names(object))){
    #   shape_by <- c(shape_by,rep(group,length(samples_names(object)[[group]])))
    # }
    shape_by <- samples_metadata(object)$group
    
    
    # Option 2: by a metadata column in object@samples$metadata
  } else if ((length(shape_by) == 1) && is.character(shape_by) & (shape_by %in% colnames(samples_metadata(object)))) {
    shape_by <- samples_metadata(object)[,shape_by]
    
    # Option 3: by a feature present in the training data    
  } else if ((length(shape_by) == 1) && is.character(shape_by) && (shape_by[1] %in% unlist(features_names(object)))) {
    # data <- lapply(get_data(object), function(l) Reduce(cbind, l))
    # features <- lapply(data, rownames)
    # viewidx <- which(sapply(features, function(x) shape_by %in% x))
    # shape_by <- data[[viewidx]][shape_by,]
    
    viewidx <- which(sapply(features_names(object), function(x) shape_by %in% x))
    foo <- list(shape_by); names(foo) <- names(viewidx)
    shape_by <- lapply(get_data(object, features = foo), function(l) Reduce(cbind, l))[[1]][1,]
    
    
    # Option 4: input is a data.frame with columns (sample,color)
  } else if (is(shape_by,"data.frame")) {
    stopifnot(all(colnames(shape_by) %in% c("sample","color")))
    stopifnot(all(unique(shape_by$sample) %in% unlist(samples_names(object))))
    
    # Option 5: shape_by is a vector of length N
  } else if (length(shape_by) > 1) {
    stopifnot(length(shape_by) == sum(object@dimensions[["N"]]))
    
    # Option not recognised
  } else {
    stop("'shape_by' was specified but it was not recognised, please read the documentation")
  }
  
  # Create data.frame with columns (sample,shape)
  if (!is(shape_by,"data.frame")) {
    df = data.frame(
      sample = samples_metadata(object)$sample,
      # sample = unlist(samples_names(object)),
      shape_by = as.factor(shape_by),
      stringsAsFactors = FALSE
    )
  }
  
  return(df)
}



.add_legend <- function(p, df, legend, color_name, shape_name) {
  
  # Add legend for color
  if (is.numeric(df$color_by)) {
    p <- p + 
      # guides(color="none") +
      scale_fill_gradientn(colors=colorRampPalette(rev(brewer.pal(n=5, name="RdYlBu")))(10))  +
      # scale_fill_gradientn(colours = c('lightgrey', 'blue'))
      labs(fill=color_name)
      
  } else {
    if (length(unique(df$color_by))>1) {
      p <- p +
        guides(fill=guide_legend(override.aes = list(shape=21))) +
        labs(fill=color_name)
    } else {
      p <- p + guides(fill="none", color="none") +
        scale_color_manual(values="black") +
        scale_fill_manual(values="gray60")
    }
    
  }
  
  # Add legend for shape
  if (length(unique(df$shape_by))>1) { 
    p <- p + 
      scale_shape_manual(values=c(21,23,24,25)[1:length(unique(df$shape_by))]) +
      guides(shape = guide_legend(override.aes = list(fill = "black"))) +
      labs(shape=shape_name)
  } else { 
    p <- p + 
      scale_shape_manual(values=c(21)) +
      guides(shape="none") 
  }
  
  # Add legend theme
  if (legend) {
    
    p <- p + 
      guides(color=guide_legend(override.aes = list(fill="white"))) +
      theme(
        legend.text = element_text(size=rel(0.8)),
        legend.title = element_text(size=rel(0.8)),
        legend.key = element_rect(fill = "white", color="white")
        # legend.background = element_rect(color = NA, fill=NA),
        # legend.box.background = element_blank()
      )
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}

# Function to define the stroke for each dot
.select_stroke <- function(N) {
  if (N<=1000) { 
    stroke <- 0.5 
  } else if (N>1000 & N<=10000) { 
    stroke <- 0.2
  } else { 
    stroke <- 0.05
  }
}

# # (Hidden) function to define the shape
# .set_shapeby_features <- function(object, shape_by, view) {
#   
#   # Option 1: no color
#   if (is.null(shape_by)) {
#     shape_by <- rep("1",sum(object@dimensions[["D"]][view]))
#     
#     # Option 2: input is a data.frame with columns (feature,color)
#   } else if (is(shape_by,"data.frame")) {
#     stopifnot(all(colnames(shape_by) %in% c("feature","color")))
#     stopifnot(all(unique(shape_by$feature) %in% features(object)[[view]]))
#     
#     # Option 3: by a feature_metadata column
#   } else if ((length(shape_by)==1) && is.character(shape_by) & (shape_by %in% colnames(features_metadata(object)))) {
#     tmp <- features_metadata(object)
#     shape_by <- tmp[tmp$view==view,shape_by]
#     
#     # Option 4: shape_by is a vector of length D
#   } else if (length(shape_by) > 1) {
#     stopifnot(length(shape_by) == object@dimensions[["D"]][[view]])
#     
#     # Option not recognised
#   } else {
#     stop("'shape_by' was specified but it was not recognised, please read the documentation")
#   }
#   
#   # Create data.frame with columns (feature,shape)
#   if (!is(shape_by,"data.frame")) {
#     df = data.frame(
#       feature = features(object)[[view]],
#       shape_by = shape_by,
#       view = view
#     )
#   }
#   
#   return(df)
# }
# 
# 
# # (Hidden) function to define the color
# .set_colorby_features <- function(object, color_by, view) {
#   
#   # Option 1: no color
#   if (is.null(color_by)) {
#     color_by <- rep("1",sum(object@dimensions[["D"]][view]))
#     
#     # Option 2: input is a data.frame with columns (feature,color)
#   } else if (is(color_by,"data.frame")) {
#     stopifnot(all(colnames(color_by) %in% c("feature","color")))
#     stopifnot(all(unique(color_by$feature) %in% features(object)[[view]]))
#     
#     # Option 3: by a feature_metadata column
#   } else if ((length(color_by)==1) && is.character(color_by) & (color_by %in% colnames(features_metadata(object)))) {
#     tmp <- features_metadata(object)
#     color_by <- tmp[tmp$view==view,color_by]
#     
#     # Option 4: color_by is a vector of length D
#   } else if (length(color_by) > 1) {
#     stopifnot(length(color_by) == object@dimensions[["D"]][[view]])
#     
#     # Option not recognised
#   } else {
#     stop("'color_by' was specified but it was not recognised, please read the documentation")
#   }
#   
#   # Create data.frame with columns (feature,color)
#   if (!is(color_by,"data.frame")) {
#     df = data.frame(
#       feature = features(object)[[view]],
#       color_by = color_by,
#       view = view
#     )
#   }
#   
#   return(df)
# }



#' @title Function to add the MOFA representation onto a Seurat object
#' @name add_mofa_factors_to_seurat
#' @description Function to add the MOFA latent representation to a Seurat object
#' @param mofa_object a trained \code{\link{MOFA}} object.
#' @param seurat_object a Seurat object
#' @param views character vector with the view names, or numeric vector with view indexes. Default is 'all'
#' @param factors character vector with the factor names, or numeric vector with the factor indexes. Default is 'all'
#' @details This function calls the \code{CreateDimReducObject} function from Seurat to store the MOFA factors.
#' @return Returns a Seurat object with the 'reductions' slot filled with the MOFA factors. Also adds, if calculated, the UMAP/TSNE obtained with the MOFA factors.
#' @export
#' @examples
#' # Generate a simulated data set
#' MOFAexample <- make_example_data()
add_mofa_factors_to_seurat <- function(mofa_object, seurat_object, views = "all", factors = "all") {
  
  # Sanity checks
  if (!is(mofa_object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package \"Seurat\" is required but is not installed.", call. = FALSE)
  }
  if (!all(colnames(seurat_object)==unlist(samples_names(mofa_object)))) {
    stop("Samples do not match between the MOFA object and the Seurat object")
  }
  
  # Get factors
  factors <- .check_and_get_factors(mofa_object, factors)
  Z <- get_factors(mofa_object, factors=factors)
  Z <- do.call("rbind",Z)
  
  # Get weights (currently not exported)
  views <- .check_and_get_views(mofa_object, views=views)
  W <- get_weights(mofa_object, views=views, factors=factors)
  
  # Collect MOFA options
  mofa_options <- list(
    "data_options" = mofa_object@data_options,
    "model_options" = mofa_object@model_options,
    "training_options" = mofa_object@training_options,
    "dimensions" = mofa_object@dimensions
  )
  
  # Sanity checks
  stopifnot(rownames(Z) %in% colnames(seurat_object))
  stopifnot(views_names(mofa_object) %in% names(seurat_object@assays))
  
  # Add to seurat
  # Add "MOFA" with no view-specific weights to the default assay 
  message("(1) Adding the MOFA factors to the 'reductions' slot of the default Seurat assay with the 'MOFA' key (no feature weights/loadings provided)...")
  seurat_object@reductions[["MOFA"]] <- CreateDimReducObject(
    embeddings = Z, 
    key = "MOFA_", 
    misc = mofa_options
  )
  
  # Add a view-specific "MOFA_" that includes the weights
  # message("(2) Adding the MOFA representation to the 'reductions' slot of each assay, including the feature weights/loadings...")
  # for (m in views_names(mofa_object)) {
  #   seurat_object@reductions[[sprintf("MOFA%s_",m)]] <- CreateDimReducObject(
  #     embeddings = Z, 
  #     loadings = W[[m]], 
  #     assay = m,
  #     key = sprintf("MOFA%s_",m), 
  #     misc = mofa_options
  #   )
  # }
  
  if (length(mofa_object@dim_red)>0) {
    if ("UMAP" %in% names(mofa_object@dim_red)) {
      message("(2) Adding the UMAP representation obtained with the MOFA factors to the 'reductions' slot of the default Seurat assay using the key 'MOFAUMAP'...")
      df <- mofa_object@dim_red$UMAP; mtx <- as.matrix(df[,-1]); rownames(mtx) <- df$sample
      colnames(df) <- paste0("MOFA_UMAP",1:ncol(df))
      seurat_object@reductions[["MOFAUMAP"]] <- CreateDimReducObject(embeddings = mtx, key = "MOFAUMAP_")
    }
    if ("TSNE" %in% names(mofa_object@dim_red)) {
      message("(2) Adding the UMAP representation obtained with the MOFA factors to the 'reductions' slot of the default Seurat assay using the key 'MOFATSNE'...")
      df <- mofa_object@dim_red$UMAP; mtx <- as.matrix(df[,-1]); rownames(mtx) <- df$sample
      seurat_object@reductions[["MOFATSNE"]] <- CreateDimReducObject(embeddings = mtx, key = "MOFATSNE_")
    }
  }
  
  return(seurat_object)
}
