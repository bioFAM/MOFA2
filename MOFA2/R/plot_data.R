
###########################################
## Functions to visualise the input data ##
###########################################

#' @title Plot heatmap of relevant features
#' @name plot_data_heatmap
#' @description Function to plot a heatmap of the input data for relevant features, 
#' usually the ones with highest loadings in a given factor.
#' @param object a \code{\link{MOFA}} object.
#' @param view a string vector with the view name, or an integer with the index of the view.
#' @param factor a string with the factor name, or an integer with the index of the factor.
#' @param groups TO-FILL
#' @param features if an integer, the total number of features to plot, based on the absolute value of the loading.
#' If a character vector, a set of manually-defined features. Default is 50.
#' @param transpose logical indicating whether to transpose the output heatmap. 
#' Default corresponds to features as rows and samples as columns.
#' @param imputed logical indicating whether to plot the imputed data instead of the original data. 
#' Default is FALSE.
#' @param denoise logical indicating whether to plot the input data or to plot a potentially denoised version based on the data reconstruction from the MOFA factors.
#' @param annotate_samples TO-FILL
#' @param annotate_features TO-FILL
#' @param ... further arguments that can be passed to \code{\link[pheatmap]{pheatmap}}
#' @details One of the first steps for the annotation of a given factor is to visualise the corresponding loadings, 
#' using for example \code{\link{plot_weights}} or \code{\link{plot_top_weights}}, which show you which are the top features that are driving the heterogeneity. \cr
#' However, one might also be interested in visualising the direct relationship between features and factors, rather than looking at "abstract" weights. \cr
#' This function generates a heatmap for selected features, which should reveal, im the original data space, the underlying pattern that is captured by the latent factor. \cr
#' A similar function for doing scatterplots rather than heatmaps is \code{\link{plot_data_scatter}}.
#' @import pheatmap
#' @export
plot_data_heatmap <- function(object, view, factor, groups = "all", features = 50, transpose = FALSE, imputed = FALSE, denoise = FALSE,
                              annotate_samples = NULL, annotate_features = NULL, ...) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")

  # Get views
  if (is.numeric(view)) view <- views_names(object)[view]
  stopifnot(view %in% views_names(object))

  # Get groups
  groups <- .check_and_get_groups(object, groups)
  
  # Get factors
  if (is.numeric(factor)) {
    factor <- factors_names(object)[factor]
  } else { 
    stopifnot(factor %in% factors_names(object)) 
  }

  # Get weights
  W <- do.call(rbind, get_weights(object, views=view, factors=factor, as.data.frame = F))
  
  # NOTE: By default concatenate all the groups
  Z <- lapply(get_factors(object)[groups], function(z) as.matrix(z[,factor]))
  Z <- do.call(rbind, Z)[,1]
  Z <- Z[!is.na(Z)]
  
  # Get imputed data
  if (imputed) {
    data <- get_imputed_data(object, view, groups)[[1]]
  } else {
    if (denoise) {
      data <- predict(object, views=view, groups=groups, factors="all")[[1]]
      # TO-DO: ADD INTERCEPTS
    } else {
      data <- get_data(object, views=view, groups=groups)[[1]]
    }
  }

  # NOTE: Currently groups are concatenated
  if (class(data) == "list") {
    data <- do.call(cbind, data)
  }

  # Select respective samples
  data <- data[,names(Z)]
  
  # Ignore samples with full missing views
  data <- data[, apply(data, 2, function(x) !all(is.na(x)))]
  
  # Define features
  if (class(features) == "numeric") {
    if (length(features) == 1) {
      features <- rownames(W)[tail(order(abs(W)), n=features)]
    } else {
      features <- rownames(W)[order(-abs(W))[features]]
    }
    stopifnot(all(features %in% features_names(object)[[view]]))  
  } else if (class(features) == "character") {
    stopifnot(all(features %in% features_names(object)[[view]]))
  } else {
    stop("Features need to be either a numeric or character vector")
  }
  data <- data[features,]
  
  # By default, sort features according to the loadings
  W.filt <- W[features,]
  order_features <- names(W.filt)[order(W.filt)]
  data <- data[order_features,]
  
  # By default, sort samples according to the factor values
  order_samples <- names(sort(Z, decreasing=T))
  order_samples <- order_samples[order_samples %in% colnames(data)]
  data <- data[,order_samples]
  
  # Transpose the data
  if (transpose) data <- t(data)

  # Define samples annotations and features annotations
  if (!is.null(annotate_samples) || !is.null(annotate_features)) {
    ann_samples  <- NULL
    ann_features <- NULL

    # Use samples annotation
    if (!is.null(annotate_samples)) {
      stopifnot(annotate_samples %in% colnames(samples_metadata(object)))
      ann_samples <- data.frame(ID = samples_metadata(object)$group_name)
      rownames(ann_samples) <- samples_metadata(object)$sample_name
    }

    # Use features annotation
    if (!is.null(annotate_features)) {
      stopifnot(annotate_features %in% colnames(features_metadata(object)))
      ann_features <- data.frame(ID = features_metadata(object)$view_name)
      rownames(ann_features) <- features_metadata(object)$feature_name
    }

    # Plot heatmap with annotations
    pheatmap::pheatmap(data, annotation_col = ann_samples, annotation_row = ann_features, ...)

  } else {
    
    # Plot heatmap without annotations
    pheatmap::pheatmap(data, ...)
  }
  
}



#' @title Scatterplots of feature values against latent factors
#' @name plot_data_scatter
#' @description Function to do a scatterplot of the feature(s) values against the latent factor values.
#' @param object a \code{\link{MOFA}} object.
#' @param view character vector with a view name, or numeric vector with the index of the view.
#' @param factor character vector with a factor name, or numeric vector with the index of the factor.
#' @param features if an integer, the total number of features to plot (10 by default). If a character vector, a set of manually-defined features.
#' @param groups TO-FILL
#' @param color_by specifies groups or values (either discrete or continuous) used to color the dots (samples). This can be either: 
#' \itemize{
#' \item (default) the string "group", it the dots with respect to their predefined groups.
#' \item a character giving the name of a feature that is present in the input data 
#' \item a character giving the same of a column in the sample metadata slot
#' \item a vector of the same length as the number of samples specifying the value for each sample. 
#' \item a dataframe with two columns: "sample" and "color"
#' }
#' @param shape_by specifies groups or values (only discrete) used to shape the dots (samples). This can be either: 
#' \itemize{
#' \item (default) the string "group": in this case, the plot will shape the dots with respect to their predefined groups.
#' \item a character giving the name of a feature that is present in the input data 
#' \item a character giving the same of a column in the sample metadata slot
#' \item a vector of the same length as the number of samples specifying the value for each sample. 
#' \item a dataframe with two columns: "sample" and "shape"
#' }
#' @param sign can be 'positive', 'negative' or 'all' to show only positive, negative or all weights, respectively. Default is 'all'.
#' @param dot_size numeric indicating dot size.
#' @param text_size numeric indicating text size.
#' @param add_lm logical indicating whether to add a linear regression line for each plot
#' @param imputed logical indicating whether to include imputed measurements
#' @param color_name name for color legend (usually only used if color_by is not a character itself).
#' @param shape_name name for shape legend (usually only used if shape_by is not a character itself).
#' @details One of the first steps for the annotation of factors is to visualise the loadings using \code{\link{plot_weights}} or \code{\link{plot_top_weights}}, 
#' which show you which features drive the heterogeneity of each factor. 
#' However, one might also be interested in visualising the direct relationship between features and factors, rather than looking at "abstract" weights. \cr
#' This function generates scatterplots of features against factors, so that you can observe the association between them. \cr
#' A similar function for doing heatmaps rather than scatterplots is \code{\link{plot_data_heatmap}}.
#' @import ggplot2 dplyr ggpubr
#' @export
plot_data_scatter <- function(object, view, factor, groups = "all", features = 10, sign="all",
                              color_by=NULL, color_name="", color_legend = TRUE,
                              shape_by=NULL, shape_name="", shape_legend = TRUE,
                              dot_size=1, text_size=5, add_lm = TRUE, imputed=FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(length(factor)==1)
  stopifnot(length(view)==1)
  if (is.numeric(view)) view <- views_names(object)[view]
  if (!view %in% views_names(object)) stop(sprintf("The view %s is not present in the object",view))

  groups <- .check_and_get_groups(object, groups)

  # Get factors
  if (is.numeric(factor)) {
    factor <- factors_names(object)[factor]
  } else { 
    stopifnot(factor %in% factors_names(object)) 
  }
      
  # Collect relevant data
  N <- get_dimensions(object)[["N"]]
  W <- get_weights(object)[[view]][,factor]
  
  if (imputed) {
    Y <- do.call(cbind, object@imputed_data[[view]][groups])
  } else {
    Y <- do.call(cbind, object@data[[view]][groups])
  }
  
  # NOTE: By default concatenate all the groups
  # Z <- lapply(get_factors(object)[groups], function(z) as.matrix(z[,factor]))
  # Z <- do.call(rbind, Z)[,1]
  # Z <- Z[!is.na(Z)]
  Z <- get_factors(object, factors=factor, groups = groups, as.data.frame=T)
  Z <- Z[,c("sample","value")]
  colnames(Z) <- c("sample","x")
  
  # Get features
  if (sign=="all") {
    W <- abs(W)
  } else if (sign=="positive") {
    W <- W[W>0]
  } else if (sign=="negative") {
    W <- W[W<0]
  }
  
  if (class(features) == "numeric") {
    if (length(features) == 1) {
      features <- names(tail(sort(abs(W)), n=features))
    } else {
      features <- names(sort(-abs(W))[features])
    }
    stopifnot(all(features %in% features_names(object)[[view]]))  
  } else if (class(features)=="character") {
    stopifnot(all(features %in% features_names(object)[[view]]))
  } else {
    stop("Features need to be either a numeric or character vector")
  }
  W <- W[features]
  Y <- Y[features,,drop=F]
  
  # Set group/color/shape
  if (length(color_by)==1 & is.character(color_by)) color_name <- color_by
  if (length(shape_by)==1 & is.character(shape_by)) shape_name <- shape_by
  color_by <- .set_colorby(object, color_by)
  shape_by <- .set_shapeby(object, shape_by)
  
  # Merge factor values with color and shape information
  df1 <- merge(Z, color_by, by="sample")
  df1 <- merge(df1, shape_by, by="sample")
  
  # Create data frame 
  # df1 <- data.frame(sample = names(Z), x = Z, shape_by = shape_by, color_by = color_by, stringsAsFactors = F)
  df2 <- get_data(object, views = view, groups = groups, features = list(features), as.data.frame = T)
  df <- left_join(df1, df2, by = "sample")
  
  #remove values missing color or shape annotation
  # if(!showMissing) df <- df[!(is.nan(df$shape_by) & !(is.nan(df$color_by))]
  
  # Generate plot
  p <- ggplot(df, aes_string(x = "x", y = "value", color = "color_by", shape = "shape_by")) + 
    geom_point(size=dot_size) +
    scale_shape_manual(values=c(19,1,2:18)[1:length(unique(shape_by))]) +
    labs(x="Factor values", y="") +
    facet_wrap(~feature, scales="free_y") +
    theme_classic() + theme(
      axis.text = element_text(size = rel(1), color = "black"), 
      axis.title = element_text(size = rel(1.0), color="black"), 
      legend.key = element_rect(fill = "white")
    )
  
  # Add linear regression line
  if (add_lm) {
    p <- p +
      stat_smooth(method="lm", color="grey", fill="grey", alpha=0.5) +
      stat_cor(method = "pearson", label.sep="\n", output.type = "latex", label.y = max(df$value,na.rm=T), size=text_size, color="black")
  }
  
  # Add legend for color
  if (length(unique(df$color_by))>1) {
    if (color_legend) { 
      p <- p + labs(color = color_name) 
    } else { 
      p <- p + guides(color = FALSE)
    }
  } else {
    p <- p + guides(color = FALSE) + scale_colour_manual(values="black")
  }
  
  # Add legend for shape
  if (length(unique(df$shape))>1 & shape_legend) { 
    p <- p + labs(shape=shape_name)
  } else { 
    p <- p + guides(shape=F) 
  }
  
  return(p)
}



#' @title Overview of the input data
#' @name plot_data_overview
#' @description Function to do a tile plot showing the missing value structure of the input data
#' @param object a \code{\link{MOFA}} object.
#' @param colors a vector specifying the colors per view (see example for details).
#' @details This function is helpful to get an overview of the structure of the data. 
#' It shows the number of samples, groups, views and features and it indicates which measurements are missing.
#' @import ggplot2 dplyr reshape2
#' @export
plot_data_overview <- function(object, colors = NULL) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Collect relevant data
  data <- object@data
  
  M <- get_dimensions(object)[["M"]]
  G <- get_dimensions(object)[["G"]]
  if (M==1 & G==1) stop("This function is not useful when there is just one view and one group")
  if (is.null(dim(data[[1]][[1]]))) stop("Data not found")
  
  # Define colors  
  if (is.null(colors)) {
    # colors <- rep("#5CACEE", M)
    palette <- c("#D95F02", "#377EB8", "#E6AB02", "#31A354", "#7570B3", "#E7298A", "#66A61E",
                 "#A6761D", "#666666", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
                 "#A65628", "#F781BF", "#1B9E77")
    if (M < 17) colors <- palette[1:M] else colors <- rainbow(M)
    names(colors) <- views_names(object)
  } else {
      stopifnot(sort(names(colors))==sort(views_names(object)))
  }
  if (length(colors) != M) stop("Length of 'colors' does not match the number of views")

  # Define availability binary matrix to indicate whether assay j is profiled in sample i
  tmp <- lapply(data, function(m) sapply(m, function(g) apply(g, 2, function(x) !all(is.na(x)))))
  ovw <- do.call(cbind, lapply(1:M, function(m) {
    do.call(rbind, lapply(tmp[[m]], as.data.frame))
  }))
  rownames(ovw) <- object@samples_metadata$sample_name
  colnames(ovw) <- views_names(object)

  ovw$sample <- object@samples_metadata$sample_name
  ovw$group <- object@samples_metadata$group_name

  # Melt to data.frame
  molten_ovw <- reshape2::melt(ovw, id.vars = c("sample", "group"), var=c("view"))
  molten_ovw$sample <- factor(molten_ovw$sample, levels = rownames(ovw))

  n <- length(unique(molten_ovw$sample))
  
  # Add number of samples and features per view/group
  molten_ovw$combi  <- ifelse(molten_ovw$value, as.character(molten_ovw$view), "missing")
  molten_ovw$ntotal <- paste("N=", sapply(data[[1]], function(e) ncol(e))[ as.character(molten_ovw$group) ], sep="")
  molten_ovw$ptotal <- paste("D=", sapply(data, function(e) nrow(e[[1]]))[ as.character(molten_ovw$view) ], sep="")
    
  # Define y-axis label
  molten_ovw <- mutate(molten_ovw, view_label = paste(view, ptotal, sep="\n"), group_label = paste(group, ntotal, sep="\n"))
  
  # Plot
  p <- ggplot(molten_ovw, aes_string(x="sample", y="view_label", fill="combi")) +
    geom_tile() +
    scale_fill_manual(values = c("missing"="grey", colors)) +
    # xlab(paste0("Samples (N=", n, ")")) + ylab("") +
    guides(fill=F) + 
    facet_wrap(~group_label, scales="free_x", nrow=length(unique(molten_ovw$view_label))) +
    theme(
      panel.background = element_rect(fill="white"),
      text = element_text(size=14),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color="black"),
      strip.background = element_blank(),
      panel.grid = element_blank()
    )
  
  return(p)
}

#' @title Visualize the structure of the data in the terminal
#' @name plot_ascii_data
#' @description A Fancy printing method
#' @param object a \code{\link{MOFA}} object
#' @param header a logical value specifying whether to show the MOFA header.
#' @details This function is helpful to get an overview of the structure of the data as a text output
#' @export
plot_ascii_data <- function(object, header = FALSE) {
  stopifnot(class(object) == "MOFA")

  if (!.hasSlot(object, "dimensions") | length(object@dimensions) == 0)
    stop("Error: dimensions not defined")
  if (!.hasSlot(object, "status") | length(object@status) == 0)
    stop("Error: status not defined")

  vis_lines <- ""

  lpad <- max(sapply(views_names(object), function(v) nchar(v)))
  wlim <- max(sapply(groups_names(object), function(v) nchar(v)))
  igr_sp <- .rep_string(5, " ")
  s <- 8             # extra lpadding shift
  w <- max(8, wlim)  # width of one block (minus 2 walls)
  hat    <- paste0(" ", .rep_string(w, "_"), " ")
  walls  <- paste0("|", .rep_string(w, " "), "|")
  ground <- paste0("|", .rep_string(w, "_"), "|")

  if (header) {
    cat("
         \U2588︎\U2588︎\U2588︎\U2588︎\U2588︎     \U2588︎\U2588︎   \U2588\U2588︎\U2588︎\U2588︎\U2588︎
mofa   \U2588︎\U2588︎\U2588︎\U2588︎\U2588︎  =  \U2588︎\U2588︎ x \U2588︎\U2588︎\U2588︎\U2588︎\U2588︎
         \U2588︎\U2588︎\U2588︎\U2588︎\U2588︎     \U2588︎\U2588︎   
    ")
  }

  groups_line    <- .pad_left(lpad + s, .cpaste(groups_names(object), w+2, collapse = igr_sp))
  nsamples_line  <- .pad_left(lpad + s, .cpaste(get_dimensions(object)$N, w+2, collapse = igr_sp))
  vis_lines      <- c(vis_lines, groups_line, nsamples_line) 

  # Calculate percentage of missing values in every view and every group
  content_pct <- lapply(object@data, function(view) sapply(view, function(group) sum(is.na(group))))
  if (length(content_pct) == 0) {
    content_pct <- lapply(object@data, function(view) sapply(view, function(group) sum(is.na(group))))
  }
  content_pct <- lapply(1:length(content_pct), function(m) {
    paste0(as.character(100 - content_pct[[m]] / object@dimensions$N / object@dimensions$D[m] * 100), sep = "%")
  })

  for (m in 1:length(views_names(object))) {
    # browser()
    toprect_line   <- .pad_left(lpad + s, paste(.rep_string(get_dimensions(object)$G, hat, collapse = igr_sp)))
    midrect_line   <- .pad_left(lpad + s, paste(.rep_string(get_dimensions(object)$G, walls, collapse = igr_sp)))
    dfeatures_line <- .pad_left_with(lpad + s, 
                                     paste(.insert_inside(content_pct[[m]], rep(walls, get_dimensions(object)$G)), collapse = igr_sp), 
                                     with = paste(c(views_names(object)[m], .cpaste(get_dimensions(object)$D[m], s)), collapse = ""))
    botrect_line   <- .pad_left(lpad + s, paste(.rep_string(get_dimensions(object)$G, ground, collapse = igr_sp)))

    vis_lines      <- c(vis_lines, toprect_line, midrect_line, dfeatures_line, botrect_line)  
  }

  cat(paste(vis_lines, collapse = "\n"))

  cat("\n\n")  
}

.rep_string <- function(times, string, collapse = "") {
  paste(replicate(times, string), collapse = collapse)
}

.pad_left_with <- function(len, string, with = "") {
  wlen <- nchar(with)
  len  <- max(len - wlen, 0)
  paste0(with, paste(replicate(len, " "), collapse = ""), string)
}

.pad_left <- function(len, string) {
  .pad_left_with(len, string, with = "")
}

.insert_inside <- function(values, boxes) {
  sapply(1:length(boxes), function(i) {
    box <- boxes[i]
    v <- values[i]
    paste0(substr(box, 1, 1), .cpaste(v, nchar(box) - 2), substr(box, length(box), length(box)))
  })
}

# Center and paste
.cpaste <- function(vals, cwidth, collapse = "") {
  vals <- sapply(vals, function(e) {
    e <- toString(e)
    lendiff <- cwidth - nchar(e)
    if (lendiff > 1) {
      paste0(.rep_string(ceiling(lendiff / 2), " "),
             e,
             .rep_string(floor(lendiff / 2), " "))
    } else {
      e
    }
  })
  paste(vals, collapse = collapse)
}


# (Hidden) function to define the color
.set_colorby <- function(object, color_by) {
  
  # Option 0: no color
  if (is.null(color_by)) {
    color_by <- rep("1",sum(object@dimensions[["N"]]))
    
    # Option 1: by default group
  } else if (color_by[1] == "group") {
    color_by <- groups(object)$group_name
    
    # Option 2: by a feature present in the training data    
  } else if ((length(color_by) == 1) && is.character(color_by) && (color_by[1] %in% unlist(features_names(object)))) {
    data <- lapply(get_data(object), function(l) Reduce(cbind, l))
    features_names <- lapply(data, rownames)
    viewidx <- which(sapply(features_names, function(x) color_by %in% x))
    color_by <- data[[viewidx]][color_by,]
    
    # Option 3: by a metadata column in object@samples$metadata
  } else if ((length(color_by) == 1) && is.character(color_by) & (color_by[1] %in% colnames(samples_metadata(object)))) {
    color_by <- samples_metadata(object)[,color_by]
    
    # Option 4: input is a data.frame with columns (sample, color)
  } else if (is(color_by, "data.frame")) {
    stopifnot(all(colnames(color_by) %in% c("sample", "color")))
    stopifnot(all(unique(color_by$sample) %in% unlist(samples_names(object))))
    
    # Option 5: color_by is a vector of length N
  } else if (length(color_by) > 1) {
    stopifnot(length(color_by) == sum(get_dimensions(object)$N))
    
    # Option not recognised
  } else {
    stop("'color_by' was specified but it was not recognised, please read the documentation")
  }
  
  # Create data.frame with columns (sample,color)
  if (!is(color_by,"data.frame")) {
    df = data.frame(
      sample = unlist(samples_names(object)),
      color_by = color_by,
      stringsAsFactors = F
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
    shape_by = c()
    for (group in names(samples_names(object))){
      shape_by <- c(shape_by,rep(group,length(samples_names(object)[[group]])))
    }
    
    # Option 2: by a feature present in the training data    
  } else if ((length(shape_by) == 1) && is.character(shape_by) && (shape_by[1] %in% unlist(features_names(object)))) {
    data <- lapply(get_data(object), function(l) Reduce(cbind, l))
    features_names <- lapply(data, rownames)
    viewidx <- which(sapply(features_names, function(x) shape_by %in% x))
    shape_by <- data[[viewidx]][shape_by,]
    
    # Option 3: input is a data.frame with columns (sample,color)
  } else if (is(shape_by,"data.frame")) {
    stopifnot(all(colnames(shape_by) %in% c("sample","color")))
    stopifnot(all(unique(shape_by$sample) %in% unlist(samples_names(object))))
    
    # Option 4: by a metadata column in object@samples$metadata
  } else if ((length(shape_by) == 1) && is.character(shape_by) & (shape_by %in% colnames(samples_metadata(object)))) {
    shape_by <- samples_metadata(object)[,shape_by]
    
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
      sample = unlist(samples_names(object)),
      shape_by = as.factor(shape_by),
      stringsAsFactors = FALSE
    )
  }
  
  return(df)
}
