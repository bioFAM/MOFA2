
###########################################
## Functions to visualise the input data ##
###########################################



#' @title Plot heatmap of relevant features
#' @name plot_data_heatmap
#' @description Function to plot a heatmap of the data for relevant features, typically the ones with high weights.
#' @param object a \code{\link{MOFA}} object.
#' @param factor a string with the factor name, or an integer with the index of the factor.
#' @param view a string with the view name, or an integer with the index of the view. Default is the first view.
#' @param groups groups to plot. Default is "all".
#' @param features if an integer (default), the total number of features to plot based on the absolute value of the weights.
#' If a character vector, a set of manually defined features.
#' @param annotation_samples annotation metadata for samples (columns). 
#' Either a character vector specifying columns in the sample metadata, or a data.frame that will be passed to \code{\link[pheatmap]{pheatmap}} as \code{annotation_row}
#' @param annotation_features annotation metadata for features (rows). 
#' Either a character vector specifying columns in the feature metadata, or a data.frame that will be passed to \code{\link[pheatmap]{pheatmap}} as \code{annotation_col}
#' @param transpose logical indicating whether to transpose the heatmap. 
#' Default corresponds to features as rows and samples as columns.
#' @param imputed logical indicating whether to plot the imputed data instead of the original data. Default is FALSE.
#' @param denoise logical indicating whether to plot a denoised version of the data reconstructed using the MOFA factors. 
#' @param max.value numeric indicating the maximum value to display in the heatmap (i.e. the matrix values will be capped at \code{max.value} ).
#' @param min.value numeric indicating the minimum value to display in the heatmap (i.e. the matrix values will be capped at \code{min.value} ).
#' See \code{\link{predict}}. Default is FALSE.
#' @param ... further arguments that can be passed to \code{\link[pheatmap]{pheatmap}}
#' @details One of the first steps for the annotation of a given factor is to visualise the corresponding weights, 
#' using for example \code{\link{plot_weights}} or \code{\link{plot_top_weights}}. \cr
#' However, one might also be interested in visualising the direct relationship between features and factors, rather than looking at "abstract" weights. \cr
#' This function generates a heatmap for selected features, which should reveal the underlying pattern that is captured by the latent factor. \cr
#' A similar function for doing scatterplots rather than heatmaps is \code{\link{plot_data_scatter}}.
#' @return A  \code{\link[pheatmap]{pheatmap}} object
#' @importFrom pheatmap pheatmap
#' @importFrom utils tail
#' @export
#' @examples
#' # Using an existing trained model
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' plot_data_heatmap(model, factor = 1, show_rownames = FALSE, show_colnames = FALSE)

plot_data_heatmap <- function(object, factor, view = 1, groups = "all", features = 50, 
    annotation_features = NULL, annotation_samples = NULL, transpose = FALSE, 
    imputed = FALSE, denoise = FALSE, max.value = NULL, min.value = NULL, ...) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(length(factor)==1)
  stopifnot(length(view)==1)
  
  # Define views, factors and groups
  groups <- .check_and_get_groups(object, groups)
  factor <- .check_and_get_factors(object, factor)
  view <- .check_and_get_views(object, view)
  
  # Get weights
  W <- do.call(rbind, get_weights(object, views=view, factors=factor, as.data.frame = FALSE))
  
  # NOTE: By default concatenate all the groups
  Z <- lapply(get_factors(object)[groups], function(z) as.matrix(z[,factor]))
  Z <- do.call(rbind, Z)[,1]
  Z <- Z[!is.na(Z)]
  

  # Get data
  if (isTRUE(denoise)) {
    data <- predict(object, views=view, groups=groups)[[1]]
  } else {
    if (isTRUE(imputed)) {
      data <- get_imputed_data(object, view, groups)[[1]]
    } else {
      data <- get_data(object, views=view, groups=groups)[[1]]
    }
  }

  # Concatenate groups
  if (is(data, "list")) {
    data <- do.call(cbind, data)
  }
  
  # Subset features
  if (is(features, "numeric")) {
    if (length(features)==1) {
      features <- rownames(W)[tail(order(abs(W)), n=features)]
    } else {
      features <- rownames(W)[order(-abs(W))[features]]
    }
    # Sort features according to the weights
    features <- names(W[features,])[order(W[features,])]
  } else if (is(features, "character")) {
    stopifnot(all(features %in% features_names(object)[[view]]))
  } else {
    stop("Features need to be either a numeric or character vector")
  }
  data <- data[features,]
  

  # Select respective samples
  data <- data[,names(Z)]
  
  # Ignore samples with full missing views
  data <- data[, apply(data, 2, function(x) !all(is.na(x)))]
  
  # By default, sort samples according to the factor values
  order_samples <- names(sort(Z, decreasing = TRUE))
  order_samples <- order_samples[order_samples %in% colnames(data)]
  data <- data[,order_samples]
    
  # Add sample annotations
  if (!is.null(annotation_samples)) {
    
    # Predefined data.frame
    if (is.data.frame(annotation_samples)) {
      message("'annotation_samples' provided as a data.frame, please make sure that the rownames match the sample names")
      if (any(!colnames(data)%in%rownames(annotation_samples))) {
        stop("There are rownames in annotation_samples that do not correspond to sample names in the model")
      }
      annotation_samples <- annotation_samples[colnames(data), , drop = FALSE]
      
    # Extract metadata from the sample metadata  
    } else if (is.character(annotation_samples)) {
      stopifnot(annotation_samples%in%colnames(object@samples_metadata))
      # tmp <- tibble::column_to_rownames(object@samples_metadata,"sample")[order_samples,,drop=F]
      tmp <- object@samples_metadata
      rownames(tmp) <- tmp$sample
      tmp$sample <- NULL
      tmp <- tmp[order_samples,,drop=FALSE]
      annotation_samples <- tmp[,annotation_samples, drop=FALSE]
      rownames(annotation_samples) <- rownames(tmp)
    } else {
      stop("Input format for 'annotation_samples' not recognised ")
    }
    
    # Convert character columns to factors
    foo <- sapply(annotation_samples, function(x) is.logical(x) || is.character(x))
    if (any(foo)) annotation_samples[,which(foo)] <- lapply(annotation_samples[,which(foo),drop=FALSE], as.factor)
  }

  
  # Add feature annotations
  if (!is.null(annotation_features)) {
    stop("'annotation_features' is currently not implemented")
  }
  
  # Transpose the data
  if (transpose) {
    data <- t(data)
    if (!is.null(annotation_samples)) {
      annotation_features <- annotation_samples
      annotation_samples <- NULL
    }
    if (!is.null(annotation_features)) {
      annotation_samples <- annotation_features
      annotation_features <- NULL
    }
  }
  
  # Cap values
  if (!is.null(max.value)) data[data>=max.value] <- max.value
  if (!is.null(min.value)) data[data<=min.value] <- min.value
  
  # Plot heatmap
  pheatmap(data, 
    annotation_row = annotation_features, 
    annotation_col = annotation_samples, 
    ...
  )
  
}



#' @title Scatterplots of feature values against latent factors
#' @name plot_data_scatter
#' @description Function to do a scatterplot of features against factor values.
#' @param object a \code{\link{MOFA}} object.
#' @param factor string with the factor name, or an integer with the index of the factor.
#' @param view string with the view name, or an integer with the index of the view. Default is the first view.
#' @param groups groups to plot. Default is "all".
#' @param features if an integer (default), the total number of features to plot. If a character vector, a set of manually-defined features.
#' @param sign can be 'positive', 'negative' or 'all' (default) to show only positive, negative or all weights, respectively.
#' @param color_by specifies groups or values (either discrete or continuous) used to color the dots (samples). This can be either: 
#' \itemize{
#' \item the string "group": dots are coloured with respect to their predefined groups.
#' \item a character giving the name of a feature that is present in the input data 
#' \item a character giving the same of a column in the sample metadata slot
#' \item a vector of the same length as the number of samples specifying the value for each sample. 
#' \item a dataframe with two columns: "sample" and "color"
#' }
#' @param shape_by specifies groups or values (only discrete) used to shape the dots (samples). This can be either: 
#' \itemize{
#' \item the string "group": dots are shaped with respect to their predefined groups.
#' \item a character giving the name of a feature that is present in the input data 
#' \item a character giving the same of a column in the sample metadata slot
#' \item a vector of the same length as the number of samples specifying the value for each sample. 
#' \item a dataframe with two columns: "sample" and "shape"
#' }
#' @param legend logical indicating whether to add a legend
#' @param dot_size numeric indicating dot size (default is 5).
#' @param text_size numeric indicating text size (default is 5).
#' @param stroke numeric indicating the stroke size (the black border around the dots, default is NULL, infered automatically).
#' @param alpha numeric indicating dot transparency (default is 1).
#' @param add_lm logical indicating whether to add a linear regression line for each plot
#' @param lm_per_group logical indicating whether to add a linear regression line separately for each group
#' @param imputed logical indicating whether to include imputed measurements
#' @details One of the first steps for the annotation of factors is to visualise the weights using \code{\link{plot_weights}} or \code{\link{plot_top_weights}}.
#' However, one might also be interested in visualising the direct relationship between features and factors, rather than looking at "abstract" weights. \cr
#' A similar function for doing heatmaps rather than scatterplots is \code{\link{plot_data_heatmap}}.
#' @import ggplot2
# #' @importFrom ggpubr stat_cor
#' @importFrom dplyr left_join
#' @importFrom utils tail
#' @importFrom stats quantile
#' @return A \code{\link{ggplot}} object
#' @export
#' @examples
#' # Using an existing trained model
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' plot_data_scatter(model)

plot_data_scatter <- function(object, factor = 1, view = 1, groups = "all", features = 10, sign = "all",
                              color_by = "group", legend = TRUE, alpha = 1, shape_by = NULL, stroke = NULL,
                              dot_size = 2.5, text_size = NULL, add_lm = TRUE, lm_per_group = TRUE, imputed = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(length(factor)==1)
  stopifnot(length(view)==1)
  if (lm_per_group) add_lm = TRUE
  
  # Define views, factors and groups
  groups <- .check_and_get_groups(object, groups)
  factor <- .check_and_get_factors(object, factor)
  view <- .check_and_get_views(object, view)

  # Collect relevant data
  N <- get_dimensions(object)[["N"]]
  W <- get_weights(object)[[view]][,factor]
  
  if (imputed) {
    Y <- do.call(cbind, object@imputed_data[[view]][groups])
  } else {
    Y <- do.call(cbind, object@data[[view]][groups])
  }
  
  # Feetch factors
  Z <- get_factors(object, factors = factor, groups = groups, as.data.frame = TRUE)
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
  
  if (is(features, "numeric")) {
    if (length(features) == 1) {
      features <- names(tail(sort(abs(W)), n=features))
    } else {
      features <- names(sort(-abs(W))[features])
    }
    stopifnot(all(features %in% features_names(object)[[view]]))  
  } else if (is(features, "character")) {
    stopifnot(all(features %in% features_names(object)[[view]]))
  } else {
    stop("Features need to be either a numeric or character vector")
  }
  W <- W[features]
  Y <- Y[features,,drop = FALSE]
  
  # Set group/color/shape
  if (length(color_by)==1 & is.character(color_by)) color_name <- color_by
  if (length(shape_by)==1 & is.character(shape_by)) shape_name <- shape_by
  color_by <- .set_colorby(object, color_by)
  shape_by <- .set_shapeby(object, shape_by)
  
  # Merge factor values with color and shape information
  df1 <- merge(Z, color_by, by="sample")
  df1 <- merge(df1, shape_by, by="sample")
  
  # Create data frame 
  foo <- list(features); names(foo) <- view
  df2 <- get_data(object, groups = groups, features = foo, as.data.frame = TRUE)
  df2$sample <- as.character(df2$sample)
  df <- dplyr::left_join(df1, df2, by = "sample")
  
  # (Q) Remove samples with missing values in Factor values
  df <- df[!is.na(df$value),]
  
  # Set stroke
  if (is.null(stroke)) {
    stroke <- .select_stroke(N=length(unique(df$sample)))
  }
  
  # Set Pearson text size
  if (add_lm && is.null(text_size)) {
    text_size <- .select_pearson_text_size(N=length(unique(df$feature)))
  }
  
  # Set axis text size
  axis.text.size <- .select_axis.text.size(N=length(unique(df$feature)))
  
  # Generate plot
  p <- ggplot(df, aes_string(x = "x", y = "value")) + 
    geom_point(aes_string(fill = "color_by", shape = "shape_by"), colour = "black", size = dot_size, stroke = stroke, alpha = alpha) +
    labs(x="Factor values", y="") +
    facet_wrap(~feature, scales="free_y") +
    theme_classic() + 
    theme(
      axis.text = element_text(size = rel(axis.text.size), color = "black"), 
      axis.title = element_text(size = rel(1.0), color="black")
    )

  # Add linear regression line
  if (add_lm) {
    if (lm_per_group && length(groups)>1) {
      p <- p +
        stat_smooth(formula=y~x, aes_string(color="group"), method="lm", alpha=0.4) +
        ggpubr::stat_cor(aes_string(color="group", label = "..r.label.."), method = "pearson", label.sep="\n", output.type = "latex", size = text_size)# +
        # guides(color = FALSE)
    } else {
      p <- p +
        stat_smooth(formula=y~x, method="lm", color="grey", fill="grey", alpha=0.4) +
        ggpubr::stat_cor(method = "pearson", label.sep="\n", output.type = "latex", size = text_size, color = "black")
    }
  }
  
  # Add legend
  p <- .add_legend(p, df, legend, color_name, shape_name)
  
  return(p)
}



#' @title Overview of the input data
#' @name plot_data_overview
#' @description Function to do a tile plot showing the missing value structure of the input data
#' @param object a \code{\link{MOFA}} object.
#' @param colors a vector specifying the colors per view (see example for details).
#' @param show_dimensions logical indicating whether to plot the dimensions of the data (default is TRUE).
#' @details This function is helpful to get an overview of the structure of the data. 
#' It shows the model dimensionalities (number of samples, groups, views and features) 
#' and it indicates which measurements are missing.
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom dplyr mutate
#' @return A \code{\link{ggplot}} object
#' @export
#' @examples
#' # Using an existing trained model
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' plot_data_overview(model)

plot_data_overview <- function(object, colors = NULL, show_dimensions = TRUE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Collect relevant data
  data <- object@data
  
  M <- get_dimensions(object)[["M"]]
  G <- get_dimensions(object)[["G"]]
  if (M==1 & G==1) warning("This function is not useful when there is just one view and one group")
  if (is.null(dim(data[[1]][[1]]))) stop("Data not found")
  
  # Define colors  
  if (is.null(colors)) {
    # colors <- rep("#5CACEE", M)
    palette <- c("#FF7F50", "#D95F02", "#377EB8", "#E6AB02", "#31A354", "#7570B3", "#E7298A", 
                 "#66A61E", "#A6761D", "#666666", "#E41A1C", "#4DAF4A", "#984EA3", "#FF7F00", 
                 "#FFFF33", "#A65628", "#F781BF", "#1B9E77")
    if (M < 18) colors <- palette[seq_len(M)] else colors <- rainbow(M)
    names(colors) <- views_names(object)
  } else {
    if (length(colors) != M) stop("Length of 'colors' does not match the number of views")
    if(is.null(names(colors))) {
      names(colors) <- views_names(object)
    } else {
      stopifnot(sort(names(colors))==sort(views_names(object)))
    }
  }

  # Define availability binary matrix to indicate whether assay j is profiled in sample i
  tmp <- lapply(data, function(m) sapply(m, function(g) apply(g, 2, function(x) !all(is.na(x)))))
  ovw <- do.call(cbind, lapply(seq_len(M), function(m) {
    do.call(rbind, lapply(tmp[[m]], as.data.frame))
  }))
  rownames(ovw) <- object@samples_metadata$sample
  colnames(ovw) <- views_names(object)

  ovw$sample <- object@samples_metadata$sample
  ovw$group <- object@samples_metadata$group

  # Melt to data.frame
  to.plot <- melt(ovw, id.vars = c("sample", "group"), var=c("view"))
  to.plot$sample <- factor(to.plot$sample, levels = rownames(ovw))

  n <- length(unique(to.plot$sample))
  
  # Add number of samples and features per view/group
  to.plot$combi  <- ifelse(to.plot$value, as.character(to.plot$view), "missing")
  if (show_dimensions) {
    to.plot$ntotal <- paste("N=", sapply(data[[1]], function(e) ncol(e))[ as.character(to.plot$group) ], sep="")
    to.plot$ptotal <- paste("D=", sapply(data, function(e) nrow(e[[1]]))[ as.character(to.plot$view) ], sep="")
    if (length(unique(to.plot$group))==1) { 
      to.plot <- mutate(to.plot, view_label = paste(view, ptotal, sep="\n"), group_label = ntotal)
    } else {
      to.plot <- mutate(to.plot, view_label = paste(view, ptotal, sep="\n"), group_label = paste(group, ntotal, sep="\n"))
    }
  } else {
    to.plot <- mutate(to.plot, view_label = view, group_label = group)
  }
    
  # Plot
  p <- ggplot(to.plot, aes_string(x="sample", y="view_label", fill="combi")) +
    geom_tile() +
    scale_fill_manual(values = c("missing"="grey", colors)) +
    # xlab(paste0("Samples (N=", n, ")")) + ylab("") +
    guides(fill = FALSE) + 
    facet_wrap(~group_label, scales="free_x", nrow=length(unique(to.plot$view_label))) +
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
#' @param nonzero a logical value specifying whether to calculate the fraction of non-zero values (non-NA values by default)
#' @details This function is helpful to get an overview of the structure of the data as a text output
#' @return None
#' @export
#' @examples
#' # Using an existing trained model
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' plot_ascii_data(model)

plot_ascii_data <- function(object, nonzero = FALSE) {
  stopifnot(is(object, "MOFA"))

  if (!.hasSlot(object, "dimensions") || length(object@dimensions) == 0)
    stop("Error: dimensions not defined")
  if (!.hasSlot(object, "status") || length(object@status) == 0)
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

  groups_line    <- .pad_left(lpad + s, .cpaste(groups_names(object), w+2, collapse = igr_sp))
  nsamples_line  <- .pad_left(lpad + s, .cpaste(get_dimensions(object)$N, w+2, collapse = igr_sp))
  vis_lines      <- c(vis_lines, groups_line, nsamples_line) 

  # Calculate percentage of missing values in every view and every group
  if (nonzero) {
    content_pct <- lapply(object@data, function(view) sapply(view, function(group) sum(group == 0)))
  } else {
    content_pct <- lapply(object@data, function(view) sapply(view, function(group) sum(is.na(group))))
  }
  content_pct <- lapply(seq_len(length(content_pct)), function(m) {
    paste0(as.character(round(100 - content_pct[[m]] / object@dimensions$N / object@dimensions$D[m] * 100)), sep = "%")
  })

  for (m in seq_len(length(views_names(object)))) {
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
  sapply(seq_len(length(boxes)), function(i) {
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


# Function to define the axis text size for plot_data_scatter
.select_axis.text.size <- function(N) {
  if (N>=4) {
    return(0.5)
  } else if (N>=2 & N<4) {
    return(0.6)
  } else if (N==1) {
    return(0.8)
  }
}

# Function to define the text size for the pearson correlation coefficient
.select_pearson_text_size <- function(N) {
  if (N>=4) {
    return(3)
  } else if (N>=2 & N<4) {
    return(4)
  } else if (N==1) {
    return(5)
  }
}
