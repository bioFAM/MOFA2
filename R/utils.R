
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

.check_and_get_views <- function(object, views, non_gaussian=TRUE) {
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
    twodim_nodes     = c("Y", "Tau")
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
.set_colorby <- function(object, color_by) {
  
  # Option 0: no color
  if (is.null(color_by)) {
    color_by <- rep("1",sum(object@dimensions[["N"]]))
    
    # Option 1: by default group
  } else if (color_by[1] == "group") {
    color_by <- samples_metadata(object)$group
    
    # Option 2: by a metadata column in object@samples$metadata
  } else if ((length(color_by) == 1) && (is.character(color_by)|is.factor(color_by)) & (color_by[1] %in% colnames(samples_metadata(object)))) {
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
      # guides(color=FALSE) +
      scale_fill_gradientn(colors=colorRampPalette(rev(brewer.pal(n=5, name="RdYlBu")))(10))  +
      # scale_fill_gradientn(colours = c('lightgrey', 'blue'))
      labs(fill=color_name)
      
  } else {
    if (length(unique(df$color_by))>1) {
      p <- p +
        guides(fill=guide_legend(override.aes = list(shape=21))) +
        labs(fill=color_name)
    } else {
      p <- p + guides(fill=FALSE, color=FALSE) +
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
      guides(shape=FALSE) 
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
