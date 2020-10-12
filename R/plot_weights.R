########################################
## Functions to visualise the weights ##
########################################

#' @title Plot heatmap of the weights
#' @name plot_weights_heatmap
#' @description Function to visualize the weights for a given set of factors in a given view. \cr 
#' This is useful to visualize the overall pattern of the weights but not to individually characterise the factors. \cr
#' To inspect the weights of individual factors, use the functions \code{\link{plot_weights}} and \code{\link{plot_top_weights}}
#' @param object a trained \code{\link{MOFA}} object.
#' @param view character vector with the view name(s), or numeric vector with the index of the view(s) to use. 
#' Default is the first view.
#' @param features character vector with the feature name(s), or numeric vector with the index of the feature(s) to use. 
#' Default is 'all'.
#' @param factors character vector with the factor name(s), or numeric vector with the index of the factor(s) to use. 
#' Default is 'all'.
#' @param threshold threshold on absolute weight values, so that weights with a magnitude below this threshold (in all factors) are removed
#' @param ... extra arguments passed to \code{\link[pheatmap]{pheatmap}}.
#' @importFrom pheatmap pheatmap
#' @return A \code{\link{pheatmap}} object
#' @export
#' @examples 
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' plot_weights_heatmap(model)

plot_weights_heatmap <- function(object, view = 1, features = "all", factors = "all", threshold = 0, ...) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")

  if (is.numeric(view)) view <- views_names(object)[view]
  stopifnot(all(view %in% views_names(object)))  
  
  # Define factors
  factors <- .check_and_get_factors(object, factors)
  
  # Define features
  if (paste(features, collapse="") =="all") { 
    features <- features_names(object)[[view]]
  } else if (is.numeric(features)) {
    features <- features_names(object)[[view]][features]
  } else {
    stopifnot(all(features %in% features_names(object)[[view]]))  
  }

  # Get relevant data
  W <- get_weights(object, views=view, factors=factors)[[1]][features,]
  
  # apply thresholding of weights
  W <- W[!apply(W,1,function(r) all(abs(r)<threshold)),]
  W <- W[,!apply(W,2,function(r) all(abs(r)<threshold))]

  # Plot heatmap
  pheatmap(t(W), ...)
}


#' @title Scatterplots of weights
#' @name plot_weights_scatter
#' @description Scatterplot of the weights values for two factors
#' @param object a trained \code{\link{MOFA}} object.
#' @param factors a vector of length two with the factors to plot. Factors can be specified either as a characters
#' using the factor names, or as numeric with the index of the factors
#' @param view character vector with the voiew name, or numeric vector with the index of the view to use. Default is the first view.
#' @param color_by specifies groups or values used to color the features. This can be either 
#' \itemize{
#' \item a character giving the same of a column in the feature metadata slot
#' \item a vector specifying the value for each feature. 
#' \item a dataframe with two columns: "feature" and "color"
#'}
#' @param shape_by specifies groups or values used to shape the features. This can be either 
#' \itemize{
#' \item a character giving the same of a column in the feature metadata slot
#' \item a vector specifying the value for each feature. 
#' \item a dataframe with two columns: "feature" and "shape"
#'}
#' @param name_color name for color legend (usually only used if color_by is not a character itself)
#' @param name_shape name for shape legend (usually only used if shape_by is not a character itself)
#' @param show_missing logical indicating whether to include dots for which \code{shape_by} or \code{color_by} is missing
#' @param dot_size numeric indicating dot size.
#' @param abs logical indicating whether to take the absolute value of the weights.
#' @param scale logical indicating whether to scale all weights from -1 to 1 (or from 0 to 1 if \code{abs=TRUE}).
#' @param legend logical indicating whether to add a legend to the plot (default is TRUE).
#' @details One of the first steps for the annotation of factors is to visualise and group/color them using known covariates such as phenotypic or clinical data.
#' This method generates a single scatterplot for the combination of two latent factors.
#' @return Returns a \code{ggplot2} object
#' @import ggplot2
#' @export
#' @examples 
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' plot_weights_scatter(model, factors = 1:2)

plot_weights_scatter <- function (object, factors, view = 1, color_by = NULL, shape_by = NULL, dot_size = 1,  
                                 name_color="", name_shape="", show_missing = TRUE, abs = FALSE, scale = TRUE, legend = TRUE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(length(factors)==2)
  
  # Get views  
  if (is.numeric(view)) view <- views_names(object)[view]
  stopifnot(all(view %in% views_names(object))) 

  # Get factor
  if(is.numeric(factors)) {
    factors <- factors_names(object)[factors]
  } else { 
    stopifnot(all(factors %in% factors_names(object)))
  }
  
  # Collect relevant data  
  D <- object@dimensions[["D"]][view]
  W <- get_weights(object, views=view, factors=factors, as.data.frame = FALSE)
  W <- as.data.frame(W); colnames(W) <- c("x","y")
  W$view <- view
  W$feature <- features_names(object)[[view]]

  # Set color and shape
  if (length(color_by)==1 & is.character(color_by)) color_name <- color_by
  if (length(shape_by)==1 & is.character(shape_by)) shape_name <- shape_by
  color_by <- .set_colorby_features(object, color_by, view)
  shape_by <- .set_shapeby_features(object, shape_by, view)
  
  # Merge factor values with group/color/shape information
  df <- merge(W, color_by, by=c("feature","view"))
  df <- merge(df, shape_by, by=c("feature","view"))
  
  # Remove values missing color or shape annotation
  if (!show_missing) df <- df[!is.na(df$shape_by) & !is.na(df$color_by),]

  # turn into factors
  df$shape_by[is.na(df$shape_by)] <- "NA"
  df$shape_by <- as.factor(df$shape_by)
  if(length(unique(df$color_by)) < 5) df$color_by <- as.factor(df$color_by)
 
  # Calculate absolute value
  if (abs) {
    df$x <- abs(df$x)
    df$y <- abs(df$y)
  }
  
  # Scale values
  if (scale) {
    df$x <- df$x/max(abs(df$x))
    df$y <- df$y/max(abs(df$y))
  }
  
  # Create plot
  p <- ggplot(df, aes_string(x="x", y="y")) + 
    geom_point(aes_string(color = "color_by", shape = "shape_by"), size=dot_size) + 
    labs(x=factors[1], y=factors[2]) +
    geom_segment(x=min(df$x,na.rm=TRUE), xend=max(df$x,na.rm=TRUE), y=0, yend=0, size=0.25, color="orange") +
    geom_segment(y=min(df$y,na.rm=TRUE), yend=max(df$y,na.rm=TRUE), x=0, xend=0, size=0.25, color="orange") +
    theme_classic() +
    theme(
      axis.text = element_text(size=rel(1), color="black"), 
      axis.title = element_text(size=rel(1.3), color="black"), 
      axis.ticks = element_line(color="black")
    )
  
  if (scale) {
    if (abs) {
      p <- p + coord_cartesian(xlim=c(0,1), ylim=c(0,1))
    } else {
      p <- p + coord_cartesian(xlim=c(-1,1), ylim=c(-1,1))
    }
  }
  
  if (length(unique(df$color_by))==1) 
    p <- p + guides(color=FALSE) + scale_color_manual(values="black")
  
  # Add legend
  if ( (length(unique(df$color_by))>1 || length(unique(df$shape_by))>1) && legend) {
    p <- p + labs(color = name_color, shape = name_shape) + 
      theme(
      legend.key = element_rect(fill = "white"),
      legend.text = element_text(size=16),
      legend.title = element_text(size=16)
    )
  } else {
    p <- p + theme(
        legend.position = "none"
      )
  }
  
  return(p)
}


#' @title Plot distribution of feature weights (weights)
#' @name plot_weights
#' @description An important step to annotate factors is to visualise the corresponding feature weights. \cr
#' This function plots all weights for a given latent factor and view, labeling the top ones. \cr
#' In contrast, the function \code{\link{plot_top_weights}} displays only the top features with highest loading.
#' @param object a \code{\link{MOFA}} object.
#' @param view a string with the view name, or an integer with the index of the view.
#' @param factors character vector with the factor name(s), or numeric vector with the index of the factor(s).
#' @param nfeatures number of top features to label.
#' @param color_by specifies groups or values (either discrete or continuous) used to color the dots (features). This can be either: 
#' \itemize{
#' \item (default) the string "group": in this case, the plot will color the dots with respect to their predefined groups.
#' \item a character giving the name of a feature that is present in the input data 
#' \item a character giving the same of a column in the features metadata slot
#' \item a vector of the same length as the number of features specifying the value for each feature 
#' \item a dataframe with two columns: "feature" and "color"
#' }
#' @param shape_by specifies groups or values (only discrete) used to shape the dots (features). This can be either: 
#' \itemize{
#' \item (default) the string "group": in this case, the plot will shape the dots with respect to their predefined groups.
#' \item a character giving the name of a feature that is present in the input data 
#' \item a character giving the same of a column in the features metadata slot
#' \item a vector of the same length as the number of features specifying the value for each feature 
#' \item a dataframe with two columns: "feature" and "shape"
#' }
#' @param abs logical indicating whether to take the absolute value of the weights.
#' @param manual A nested list of character vectors with features to be manually labelled (see the example for details).
#' @param color_manual a character vector with colors, one for each element of 'manual'
#' @param scale logical indicating whether to scale all weights from -1 to 1 (or from 0 to 1 if abs=TRUE).
#' @param dot_size numeric indicating the dot size.
#' @param text_size numeric indicating the text size.
#' @param legend logical indicating whether to add legend.
#' @param return_data logical indicating whether to return the data frame to plot instead of plotting
#' @import ggplot2 dplyr tidyr
#' @importFrom magrittr %>%
#' @importFrom ggrepel geom_text_repel
#' @return A \code{\link{ggplot}} object or a \code{data.frame} if return_data is TRUE
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Plot distribution of weights for Factor 1 and View 1
#' plot_weights(model, view = 1, factors = 1)
#' 
#' # Plot distribution of weights for Factors 1 to 3 and View 1
#' plot_weights(model, view = 1, factors = 1:3)
#' 
#' # Take the absolute value and highlight the top 10 features
#' plot_weights(model, view = 1, factors = 1, nfeatures = 10, abs = TRUE)
#' 
#' # Change size of dots and text
#' plot_weights(model, view = 1, factors = 1, text_size = 5, dot_size = 1)
#' 
plot_weights <- function(object, view = 1, factors = 1, nfeatures = 10, 
                         color_by = NULL, shape_by = NULL,
                         abs = FALSE, manual = NULL, color_manual = NULL, scale = TRUE, 
                         dot_size = 1, text_size = 5, legend = TRUE, return_data = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(length(view)==1)
  
  # Get views
  view <- .check_and_get_views(object, view)
  
  # Get factor names
  factors <- .check_and_get_factors(object, factors)
  
  # Collect expectations  
  W <- get_weights(object, views = view, factors = factors, as.data.frame = TRUE)
  
  # Convert factor names to a factor to preserve order
  W$factor <- factor(W$factor, levels = unique(factors))


  ################
  ## Parse data ##
  ################
  
  # Scale values
  if (scale && sum(W$value>0)>0) W$value <- W$value / max(abs(W$value))
  
  # Take the absolute value
  if (abs) W$value <- abs(W$value)
    
  # Define groups for labelling
  W$labelling_group <- "0"
  
  # Define group of features to color according to the loading
  if (is.null(manual) & nfeatures>0) {
    for (f in factors) {
      features <- W[W$factor==f,] %>% group_by(view) %>% top_n(n=nfeatures, abs(value)) %>% .$feature
      W[(W$feature %in% features) & (W$factor==f), "labelling_group"] <- "1"
    }
  }
  
  # Define group of features to label manually
  if(!is.null(manual)) {
    if (is.null(color_manual)) {
      if (length(manual)>1) {
        # color_manual <- hcl(h = seq(15, 375, length=length(manual)+1), l=65, c=100)[seq_len(length(manual))]
        color_manual <- RColorBrewer::brewer.pal(n=length(manual)+1, "Dark2")
      } else {
        color_manual <- "black"
      }
    } else {
      stopifnot(length(color_manual)==length(manual)) 
    }
    
    # Add labelling group (0 for non-labeled, >= 1 for labeled)
    for (m in seq_len(length(manual)))
      W$labelling_group[W$feature %in% manual[[m]]] <- as.character(m+1)
  }
  
  # Make features names unique
  W$feature_id <- W$feature
  if ((length(unique(W$view)) > 1) && (nfeatures > 0) && (any(duplicated(W[W$factor == factors[1],]$feature_id)))) {
    message("Duplicated feature names across views, we will add the view name as a prefix")
    W$feature_id <- paste(W$view, W$feature, sep="_")
  }
  
  # labelling_indicator is TRUE for labeled, FALSE for non-labeled
  W$labelling_indicator <- as.factor(W$labelling_group != "0")

  # Set color and shape
  if (length(color_by)==1 & is.character(color_by)) color_name <- color_by
  if (length(shape_by)==1 & is.character(shape_by)) shape_name <- shape_by
  obj_color_by <- .set_colorby_features(object, color_by, view)
  obj_shape_by <- .set_shapeby_features(object, shape_by, view)
  
  # Merge factor values with group/color/shape information
  W <- merge(W, obj_color_by, by=c("feature", "view"))
  W <- merge(W, obj_shape_by, by=c("feature", "view"))

  # Sort features by weight
  W <- by(W, list(W$factor), function(x) x[order(x$value),])
  W <- do.call(rbind, W)

  # In order to re-order features across multiple factors, make them unique for different factors
  W$feature_id <- paste(W$feature_id, W$factor, sep="_")
  W$feature_id <- factor(W$feature_id, levels = unique(W$feature_id))

  # Return data if requested instead of plotting
  if (return_data) return(W)
  
  # Generate plot
  p <- ggplot(W, aes_string(x = "value", y = "feature_id", col = "labelling_group")) +
    scale_y_discrete(expand = c(0.03,0.03)) +
    geom_point(aes_string(shape = "shape_by", size="labelling_indicator")) + 
    labs(x="Weight", y="Rank", size=dot_size)
  
  # Add labels to the top features
  if (nfeatures>0 || length(unique(W$labelling_group))>0) {
    p <- p + geom_text_repel(
      force = 10,
      data = W[W$labelling_group != "0",], aes_string(label = "feature", col = "labelling_group"),
      size=text_size, segment.alpha=0.25, segment.color="black", segment.size=0.3, 
      box.padding = unit(0.5,"lines"), show.legend = FALSE)
  }
  
  # Configure axis 
  if (scale) {
    if (abs) {
      p <- p + 
        coord_cartesian(xlim=c(0,1)) +
        scale_x_continuous(breaks=c(0,1)) +
        expand_limits(x=c(0,1))
    } else {
      p <- p + 
        coord_cartesian(xlim=c(-1,1)) +
        scale_x_continuous(breaks=c(-1,0,1)) +
        expand_limits(x=c(-1,1))
    }
  }
  
  # Define dot size
  p <- p + scale_size_manual(values=c(dot_size/2,dot_size*2)) + guides(size = FALSE)
  
  # Define dot colours and legend for colours
  if (!is.null(color_by)) { 
    p <- p + labs(color=color_name)
  } else {
    foo <- c("grey","black",color_manual); names(foo) <- as.character(0:(length(foo)-1))
    p <- p + guides(color=FALSE) + scale_color_manual(values=foo)
  }
  
  # Add legend for shape
  if (!is.null(shape_by)) { 
    p <- p + labs(shape=shape_name)
  } else { 
    p <- p + guides(shape=FALSE) 
  }
  
  # Facet if multiple factors
  if (length(unique(W$factor)) > 1) {
    p <- p + facet_wrap(~factor, nrow=1, scales="free")
  }
  
  # Add Theme  
  p <- p +
    theme_bw() + 
    theme(
      plot.title = element_text(size=rel(1.3), hjust=0.5),
      axis.title = element_text(size=rel(1.3), color="black"),
      axis.text.x = element_text(size=rel(1.3), color="black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      
      # facets
      strip.text = element_text(size=rel(1.2)),
      panel.spacing = unit(1,"lines"),

      # gridlines
      panel.grid.major.y = element_blank(),
    )


  # Configure the legend
  if (legend) {
    p <- p + theme(
      legend.text = element_text(size=rel(1.2)),
      legend.title = element_text(size=rel(1.2))
    )
  } else {
      p <- p + theme(legend.position = "none")
  }
  
  return(p)
}


#' @title Plot top weights
#' @name plot_top_weights
#' @description Plot top weights for a given factor and view.
#' @param object a trained \code{\link{MOFA}} object.
#' @param view a string with the view name, or an integer with the index of the view.
#' @param factors a character string with factors names, or an integer vector with factors indices.
#' @param nfeatures number of top features to display.
#' Default is 10
#' @param abs logical indicating whether to use the absolute value of the weights (Default is FALSE).
#' @param sign can be 'positive', 'negative' or 'all' to show only positive, negative or all weights, respectively. Default is 'all'.
#' @param scale logical indicating whether to scale all weights from -1 to 1 (or from 0 to 1 if abs=TRUE). Default is TRUE.
#' @details An important step to annotate factors is to visualise the corresponding feature weights. \cr
#' This function displays the top features with highest loading whereas the function \code{\link{plot_top_weights}} plots all weights for a given latent factor and view. \cr
#' Importantly, the weights of the features within a view have relative values and they should not be interpreted in an absolute scale.
#' Therefore, for interpretability purposes we always recommend to scale the weights with \code{scale=TRUE}.
#' @import ggplot2
#' @importFrom dplyr group_by top_n desc
#' @return Returns a \code{ggplot2} object
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Plot top weights for Factors 1 and 2 and View 1
#' plot_top_weights(model, view = 1, factors = c(1,2))
#' 
#' # Do not take absolute value
#' plot_weights(model, abs = FALSE)
#' 
plot_top_weights <- function(object, view = 1, factors = 1,
                             nfeatures = 10, abs = TRUE, scale = TRUE, sign = "all") {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (nfeatures <= 0) stop("'nfeatures' has to be greater than 0")
  if (sign=="all") { abs <- TRUE}
  if (is.numeric(view)) view <- views_names(object)[view]
  stopifnot(view %in% views_names(object))
  
  # Get views
  view <- .check_and_get_views(object, view)
  
  # Get factor names
  factors <- .check_and_get_factors(object, factors)
  
  # Collect expectations  
  W <- get_weights(object, factors = factors, views = view, as.data.frame=TRUE)

  # Scale values by weight with highest (absolute) value
  if (scale) W$value <- W$value/max(abs(W$value))

  # Store sign
  W <- W[W$value!=0,]
  W$sign <- ifelse(W$value>0, "+", "-")

  # Select subset of only positive or negative weights
  if (sign=="positive") { W <- W[W$value>0,] } else if (sign=="negative") { W <- W[W$value<0,] }

  # Absolute value
  if (abs) W$value <- abs(W$value)
  
  # Extract relevant features
  W <- W[with(W, order(-abs(value))), ]

  # Sort according to weights for each factor
  W <- as.data.frame(top_n(group_by(W, factor), n = nfeatures, wt = value))
  #
  
  # Make features names unique
  W$feature_id <- W$feature
  if ((length(unique(W$view)) > 1) && (nfeatures > 0) && (any(duplicated(W[W$factor == factors[1],]$feature_id)))) {
    message("Duplicated feature names across views, we will add the view name as a prefix")
    W$feature_id <- paste(W$view, W$feature, sep="_")
  }

  # In order to re-order features across multiple factors, 
  # make them unique for different factors
  W$feature_id <- factor(W$feature_id, levels = rev(unique(W$feature_id)))
  
  p <- ggplot(W, aes_string(x="feature_id", y="value")) +
    geom_point(size=2) +
    geom_segment(aes_string(xend="feature_id"), size=0.75, yend=0) +
    scale_colour_gradient(low="grey", high="black") +
    coord_flip() +
    labs(y="Weight") +

    # Theme
    theme_bw() +
    theme(
      axis.title.x = element_text(color='black'),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.1), hjust=1, color='black'),
      axis.text.x = element_text(color='black'),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(),
      legend.position = 'top',
      legend.title = element_blank(),
      legend.text = element_text(color="black"),
      legend.key = element_rect(fill='transparent'),
      
      # facets
      strip.text = element_text(size=rel(1.2)),
      panel.background = element_blank(),
      panel.spacing = unit(1,"lines"),

      # gridlines
      panel.grid.major.y = element_blank(),
    ) +
    facet_wrap(~factor, nrow=1, scales="free")
  
  if (sign=="negative") p <- p + scale_x_discrete(position = "top")

  # If absolute values are used, add the corresponding signs to the plot
  if (abs) {
    p <- p + 
      ylim(0,max(W$value)+0.1) + 
      geom_text(label=W$sign,y=max(W$value)+0.1, size=10)
  }

  return(p)
  
}



# (Hidden) function to define the shape
.set_shapeby_features <- function(object, shape_by, view) {
  
  # Option 1: no color
  if (is.null(shape_by)) {
    shape_by <- rep("1",sum(object@dimensions[["D"]][view]))
    
  # Option 2: input is a data.frame with columns (feature,color)
  } else if (is(shape_by,"data.frame")) {
    stopifnot(all(colnames(shape_by) %in% c("feature","color")))
    stopifnot(all(unique(shape_by$feature) %in% features_names(object)[[view]]))
    
  # Option 3: by a feature_metadata column
  } else if ((length(shape_by)==1) && is.character(shape_by) & (shape_by %in% colnames(features_metadata(object)))) {
    tmp <- features_metadata(object)
    shape_by <- tmp[tmp$view==view,shape_by]
    
  # Option 4: shape_by is a vector of length D
  } else if (length(shape_by) > 1) {
    stopifnot(length(shape_by) == object@dimensions[["D"]][[view]])
    
  # Option not recognised
  } else {
    stop("'shape_by' was specified but it was not recognised, please read the documentation")
  }
  
  # Create data.frame with columns (feature,shape)
  if (!is(shape_by,"data.frame")) {
    df = data.frame(
      feature = features_names(object)[[view]],
      shape_by = shape_by,
      view = view
    )
  }
  
  return(df)
}


# (Hidden) function to define the color
.set_colorby_features <- function(object, color_by, view) {
  
  # Option 1: no color
  if (is.null(color_by)) {
    color_by <- rep("1",sum(object@dimensions[["D"]][view]))
    
    # Option 2: input is a data.frame with columns (feature,color)
  } else if (is(color_by,"data.frame")) {
    stopifnot(all(colnames(color_by) %in% c("feature","color")))
    stopifnot(all(unique(color_by$feature) %in% features_names(object)[[view]]))
    
    # Option 3: by a feature_metadata column
  } else if ((length(color_by)==1) && is.character(color_by) & (color_by %in% colnames(features_metadata(object)))) {
    tmp <- features_metadata(object)
    color_by <- tmp[tmp$view==view,color_by]
    
    # Option 4: color_by is a vector of length D
  } else if (length(color_by) > 1) {
    stopifnot(length(color_by) == object@dimensions[["D"]][[view]])
    
    # Option not recognised
  } else {
    stop("'color_by' was specified but it was not recognised, please read the documentation")
  }
  
  # Create data.frame with columns (feature,color)
  if (!is(color_by,"data.frame")) {
    df = data.frame(
      feature = features_names(object)[[view]],
      color_by = color_by,
      view = view
    )
  }
  
  return(df)
}
