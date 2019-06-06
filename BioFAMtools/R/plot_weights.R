
########################################
## Functions to visualise the weights ##
########################################

#' @title Plot heatmap of the weights
#' @name plot_weights_heatmap
#' @description Function to visualize the loadings for a given set of factors in a given view. \cr 
#' This is useful to visualize the overall pattern of the weights but not to individually characterise the factors. \cr
#' To inspect the loadings of individual factors, use the functions \code{\link{plot_weights}} and \code{\link{plot_top_weights}}
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param view character vector with the view name(s), or numeric vector with the index of the view(s) to use. 
#' Default is 'all'
#' @param features character vector with the feature name(s), or numeric vector with the index of the feature(s) to use. 
#' Default is 'all'
#' @param factors character vector with the factor name(s), or numeric vector with the index of the factor(s) to use. 
#' Default is 'all'
#' @param threshold threshold on absolute weight values, so that loadings with a magnitude below this threshold (in all factors) are removed
#' @param ... extra arguments passed to \code{\link[pheatmap]{pheatmap}}.
#' @importFrom pheatmap pheatmap
#' @export
plot_weights_heatmap <- function(object, view, features = "all", factors = "all", threshold = 0, ...) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")

  if (is.numeric(view)) view <- views_names(object)[view]
  stopifnot(all(view %in% views_names(object)))  
  
  # Get factors
  if (paste0(factors, collapse="") == "all") { 
  	factors <- factors_names(object)
  } else if (is.numeric(factors)) {
  	factors <- factors_names(object)[factors]
  } else { 
  	stopifnot(all(factors %in% factors_names(object)))
  }
  
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
  

  # Set title
  # if (is.null(main)) { main <- paste("Loadings of Latent Factors on", view) }
  
  # set colors and breaks if not specified
  # if (is.null(color) & is.null(breaks)) {
  #   palLength <- 100
  #   minV <- min(W)
  #   maxV <- max(W)
  #   color <- colorRampPalette(colors=c("black", "blue", "white", "orange","red"))(palLength)
  #   breaks <- c(seq(minV, 0, length.out=ceiling(palLength/2) + 1), 
  #               seq(maxV/palLength,maxV, length.out=floor(palLength/2)))
  # }
  
  # apply thresholding of loadings
  W <- W[!apply(W,1,function(r) all(abs(r)<threshold)),]
  W <- W[,!apply(W,2,function(r) all(abs(r)<threshold))]

  # Plot heatmap
  pheatmap::pheatmap(t(W), ...)
}


#' @title Scatterplot of weights for two factors
#' @name plot_weight_scatter
#' @description Scatterplot of the weights values for two factors
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param view character vector with the voiew name, or numeric vector with the index of the view to use.
#' @param factors a vector of length two with the factors to plot. Factors can be specified either as a characters
#' using the factor names, or as numeric with the index of the factors
#' @param color_by specifies groups or values used to color the samples. 
#' This can be either 
#' a character giving the name of a feature present in the training data, 
#' a character giving the same of a covariate (only if using \code{\link{MultiAssayExperiment}} as input), 
#' or a vector of the same length as the number of samples specifying discrete groups or continuous numeric values.
#' @param shape_by specifies groups or values used to shape the samples. 
#' This can be either
#' a character giving the name of a feature present in the training data, 
#' a character giving the same of a covariate (only if using \code{\link{MultiAssayExperiment}} as input), 
#' or a vector of the same length as the number of samples specifying discrete groups.
#' @param name_color name for color legend (usually only used if color_by is not a character itself)
#' @param name_shape name for shape legend (usually only used if shape_by is not a character itself)
#' @param showMissing logical indicating whether to include samples for which \code{shape_by} or \code{color_by} is missing
#' @details One of the first steps for the annotation of factors is to visualise and group/color them using known covariates such as phenotypic or clinical data.
#' This method generates a single scatterplot for the combination of two latent factors.
#' @return Returns a \code{ggplot2} object
#' @import ggplot2
#' @export
plot_weight_scatter <- function (object, view, factors, color_by = NULL, shape_by = NULL, 
                                 name_color="", name_shape="", showMissing = TRUE, abs=T) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  stopifnot(length(factors)==2)

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
  W <- get_weights(object, views=view, factors=factors, as.data.frame = F)
  W <- W[[view]][,factors]
  
  # Set color
  colorLegend <- T
  if (!is.null(color_by)) {
    
    # It is the name of a covariate or a feature in the TrainData
    if (length(color_by) == 1 & is.character(color_by)) {
      if(name_color=="") name_color <- color_by
      training_data <- get_training_data(object)
      features_names <- features_names(object)
      if(color_by %in% Reduce(union, features_names)) {
        viewidx <- which(sapply(features_names, function(vnm) color_by %in% vnm))
        color_by <- training_data[[viewidx]][[1]][color_by,]
      }
      
    # It is a vector of length N
    } else if (length(color_by) > 1) {
      stopifnot(length(color_by) == N)
      # color_by <- as.factor(color_by)
    } else {
      stop("'color_by' was specified but it was not recognised, please read the documentation")
    }
  } else {
    color_by <- rep(TRUE, D)
    colorLegend <- F
  }

  # Set shape
  shapeLegend <- T
  if (!is.null(shape_by)) {
    
    # It is the name of a covariate 
    if (length(shape_by) == 1 & is.character(shape_by)) {
      if(name_shape=="") name_shape <- shape_by
      features_names <- features_names(object)
      training_data <- get_training_data(object)
      if (shape_by %in% Reduce(union,features_names)) {
        viewidx <- which(sapply(features_names, function(vnm) shape_by %in% vnm))
        shape_by <- training_data[[viewidx]][[1]][shape_by,]
      }

    # It is a vector of length N
    } else if (length(shape_by) > 1) {
      stopifnot(length(shape_by) == N)
    } else {
      stop("'shape_by' was specified but it was not recognised, please read the documentation")
    }
  } else {
    shape_by <- rep(TRUE, D)
    shapeLegend <- F
  }
  
  # Create data frame to plot
  df <- data.frame(x = W[, factors[1]], y = W[, factors[2]], shape_by = shape_by, color_by = color_by)
  
  # remove values missing color or shape annotation
  if (!showMissing) df <- df[!is.na(df$shape_by) & !is.na(df$color_by),]

  #turn into factors
  df$shape_by[is.na(df$shape_by)] <- "NA"
  df$shape_by <- as.factor(df$shape_by)
  if(length(unique(df$color_by)) < 5) df$color_by <- as.factor(df$color_by)
 
  if (abs) {
    df$x <- abs(df$x)
    df$y <- abs(df$y)
  }
  
  xlabel <- paste("Latent factor", factors[1])
  ylabel <- paste("Latent factor", factors[2])
                                
  p <- ggplot(df, aes_string(x = "x", y = "y")) + 
      geom_point(aes_string(color = "color_by", shape = "shape_by")) + xlab(xlabel) + ylab(ylabel) +
      # scale_shape_manual(values=c(19,1,2:18)[1:length(unique(shape_by))]) +
      theme(plot.margin = margin(20, 20, 10, 10), 
            axis.text = element_text(size = rel(1), color = "black"), 
            axis.title = element_text(size = 16), 
            axis.title.y = element_text(size = rel(1.1), margin = margin(0, 10, 0, 0)), 
            axis.title.x = element_text(size = rel(1.1), margin = margin(10, 0, 0, 0)), 
            axis.line = element_line(color = "black", size = 0.5), 
            axis.ticks = element_line(color = "black", size = 0.5),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(),
            legend.key = element_rect(fill = "white"),
            legend.text = element_text(size = 16),
            legend.title = element_text(size =16)
            )
  if (colorLegend) { p <- p + labs(color = name_color) } else { p <- p + guides(color = FALSE) }
  if (shapeLegend) { p <- p + labs(shape = name_shape) }  else { p <- p + guides(shape = FALSE) }
  return(p)
}


#' @title Plot Weights
#' @name plot_weights
#' @description An important step to annotate factors is to visualise the corresponding feature loadings. \cr
#' This function plots all loadings for a given latent factor and view, labeling the top ones. \cr
#' In contrast, the function \code{\link{plot_top_weights}} displays only the top features with highest loading.
#' @param object a \code{\link{BioFAModel}} object.
#' @param views character vector with the voiew name, or numeric vector with the index of the view to use.
#' @param factors character vector with the factor name, or numeric vector with the index of the factor to use.
#' @param nfeatures number of top features to label.
#' @param abs logical indicating whether to use the absolute value of the weights.
#' @param manual A nested list of character vectors with features to be manually labelled.
#' @param color_manual a character vector with colors, one for each element of 'manual'
#' @param scale logical indicating whether to scale all loadings from 0 to 1.
#' @param sort_by_factor name or numeric value for the factor to sort all loadings by, or all (default)
#' when "all", features between factors do not match
#' @details The weights of the features within a view are relative andthey should not be interpreted in an absolute scale.
#' Therefore, for interpretability purposes we always recommend to scale the weights with \code{scale=TRUE}.
#' @import ggplot2 ggrepel
#' @export
plot_weights <- function(object, views = "all", factors = "all", nfeatures = 10, 
                         abs = FALSE, manual = NULL, color_manual = NULL, scale = TRUE, 
                         sort_by_factor = "all",
                         dot_size = 1, text_size = 5) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  
  views <- .check_and_get_views(object, views)

  # Get factor
  if (factors[1] == "all") {
    factors <- factors_names(object)
  } else if (is.numeric(factors)) {
    factors <- factors_names(object)[factors]
  } else {
    stopifnot(factors %in% factors_names(object)) 
  }

  # Get factor to sort by
  if (!is.null(sort_by_factor)) {
    if (is.numeric(sort_by_factor)) {
      sort_by_factor <- factors_names(object)[sort_by_factor]
    } else if (sort_by_factor == "all") {
      # 
    } else { 
      stopifnot(sort_by_factor %in% factors_names(object)) 
    }
    stopifnot(length(sort_by_factor) == 1)
  }
  
  # Collect expectations  
  W <- get_weights(object, views = views, factors = factors, as.data.frame = T)
  W <- W[(W$factor %in% factors) & (W$view %in% views),]
  # Convert factor names to a factor to preserve order
  W$factor <- factor(W$factor, levels = unique(W$factor))
  
  # Scale values
  if (scale) W$value <- W$value / max(abs(W$value))
  
  # Parse the weights
  if (abs) W$value <- abs(W$value)
    
  # Define groups for labelling
  W$group <- "0"
  
  # Define group of features to color according to the loading
  if (nfeatures > 0) {
    for (f in factors) {
      # This uses dplyr
      features <- W %>% filter(factor == f) %>% group_by(view) %>% top_n(n=nfeatures, abs(value)) %>% .$feature
      W[(W$feature %in% features) & (W$factor == f),"group"] <- "1"
    }
  }
  
  # Define group of features to label manually
  if(!is.null(manual)) {
    if (is.null(color_manual)) {
      color_manual <- hcl(h = seq(15, 375, length = length(manual) + 1), l = 65, c = 100)[1:length(manual)]
    } else {
      stopifnot(length(color_manual)==length(manual)) 
    }
    for (m in 1:length(manual))
      W$group[W$feature %in% manual[[m]]] <- as.character(m+1)
  }
  
  # Make features names unique
  W$feature_id <- W$feature
  if ((length(unique(W$view)) > 1) && (nfeatures > 0) && (any(duplicated(W[W$factor == factors[1],]$feature_id)))) {
    message("Duplicated feature names across views, we will add the view name as a prefix")
    W$feature_id <- paste(W$view, W$feature, sep="_")
  }

  # Sort by loading
  if (!is.null(sort_by_factor)) {
    W <- by(W, list(W$factor), function(x) {
      x[order(x$value),]
    })
    if (sort_by_factor == "all") {
      W <- do.call(rbind, W)
      W$feature_id <- paste(W$view, W$feature, W$factor, sep="_")
    } else {
      W <- do.call(rbind, c(W[sort_by_factor], W[names(W) != sort_by_factor]))
      W$feature_id <- W$feature
    }
  }
  W$feature_id <- factor(W$feature_id, levels = unique(W$feature_id))
  
  # Define plot title
  # if(is.null(main)) main <- paste("Distribution of weigths of LF", factor, "in", view, "view")
  
  # Generate plot

  # Convert plotting group
  W$tmp <- as.factor(W$group != "0")

  gg_W <- ggplot(W, aes(x=feature_id, y=value, col=group)) +
    # scale_y_continuous(expand = c(0.01,0.01)) + scale_x_discrete(expand = c(0.01,0.01)) +
    # scale_y_discrete(expand = c(0.01, 0.01)) +
    scale_x_discrete(expand = c(0.01, 0.01)) +
    geom_point(aes(size=tmp)) + labs(x="Rank position", y="Loading", size=dot_size) +
    ggrepel::geom_text_repel(data = W[W$group!="0",], aes(label = feature, col = group),
                             size=text_size, segment.alpha=0.1, segment.color="black", segment.size=0.3, box.padding = unit(0.5, "lines"), show.legend=F)
  
  if (scale) {
    gg_W <- gg_W + coord_cartesian(ylim=c(-1,1))
  }
  
  # Define size
  gg_W <- gg_W + scale_size_manual(values=c(dot_size,dot_size*2)) + guides(size=F)
  
  # Define colors
  cols <- c("grey", "black", color_manual)
  gg_W <- gg_W + scale_color_manual(values=cols) + guides(col=F)
  
  # Facet if multiple views and for multiple factors
  if ((length(unique(W$view)) > 1) && (length(unique(W$factor)) > 1)) {
    gg_W <- gg_W + facet_wrap(view ~ factor, scales="free")
  }
  else if (length(unique(W$factor)) > 1) {
    gg_W <- gg_W + facet_wrap(~factor, nrow=1, scales="free")
  }
  else if (length(unique(W$view)) > 1) {
    gg_W <- gg_W + facet_wrap(~view, ncol=1, scales="free")
  }
  
  
  # Add Theme  
  gg_W <- gg_W +
    theme_minimal() +
    theme_bw() + 
    theme(
      # panel.spacing = margin(5,5,5,5),
      # panel.border = element_rect(colour = "black", fill=NA, size=0.75),
      plot.title = element_text(size=rel(1.3), hjust=0.5),
      # axis.title.x = element_blank(),
      axis.text.x = element_text(size=rel(1.3), color="black"),
      axis.text.y = element_blank(),  # loadings names are hidden by default
      axis.title.y = element_text(size=rel(1.5), color="black"),
      axis.ticks.y = element_blank(),

      # white background and dark border
      # panel.background = element_rect(fill = "white", colour = NA),
      # panel.border     = element_rect(fill = NA, colour = "grey20"),
      strip.background = element_blank(),
      panel.border = element_blank(),

      # make gridlines dark, same contrast with white as in theme_grey
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    ) + 
    coord_flip() +
    expand_limits(y = c(-1, 1))  # loadings axis
  
  return(gg_W)
}


#' @title Plot top weights
#' @name plot_top_weights
#' @description Plot top weights for a given latent in a given view.
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param view character vector with the view name, or numeric vector with the index of the view to use.
#' @param factor character vector with the factor name, or numeric vector with the index of the factor to use.
#' @param nfeatures number of top features to display.
#' Default is 10
#' @param abs logical indicating whether to use the absolute value of the weights.
#' Default is TRUE
#' @param sign can be 'positive', 'negative' or 'both' to show only positive, negative or all weigths, respectively.
#' Default is 'both'.
#' @param scale logical indicating whether to scale all loadings from 0 to 1.
#' Default is TRUE.
#' @details An important step to annotate factors is to visualise the corresponding feature loadings. \cr
#' This function displays the top features with highest loading whereas the function \code{\link{plot_top_weights}} plots all loadings for a given latent factor and view. \cr
#' Importantly, the weights of the features within a view have relative values and they should not be interpreted in an absolute scale.
#' Therefore, for interpretability purposes we always recommend to scale the weights with \code{scale=TRUE}.
#' @import ggplot2
#' @return Returns a \code{ggplot2} object
#' @export
plot_top_weights <- function(object, view, factor, nfeatures = 10, abs = TRUE, scale = TRUE, sign = "both") {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  
  if (is.numeric(view)) view <- views_names(object)[view]
  stopifnot(view %in% views_names(object))
  # if(!is.null(manual_features)) { stopifnot(class(manual_features)=="list"); stopifnot(all(Reduce(intersect,manual_features) %in% features_names(object)[[view]]))  }
  
  # Collect expectations  
  W <- get_weights(object, factors=factor, views=view, as.data.frame=T)

  # Scale values by loading with highest (absolute) value
  if(scale) W$value <- W$value/max(abs(W$value))

  # store sign
  W <- W[W$value!=0,]
  W$sign <- ifelse(W$value>0, "+", "-")

 # select subset of only positive or negative loadings
  if (sign=="positive") { W <- W[W$value>0,] } else if (sign=="negative") { W <- W[W$value<0,] }

   # Absolute value
  if (abs) W$value <- abs(W$value)

  
  # Extract relevant features
  W <- W[with(W, order(-abs(value))), ]
  if (nfeatures>0) features <- head(W$feature,nfeatures) # Extract top hits
  # if (!is.null(manual_features)) features <- W$feature[W$feature %in% manual_features] # Extract manual hits
  W <- W[W$feature %in% features,]
  
  # Sort according to loadings
  W <- W[with(W, order(-value, decreasing = T)), ]
  W$feature <- factor(W$feature, levels=W$feature)
  
  p <- ggplot(W, aes(x=feature, y=value)) +
    geom_point(size=2) +
    geom_segment(aes(xend=feature, yend=0), size=0.75) +
    scale_colour_gradient(low="grey", high="black") +
    # scale_colour_manual(values=c("#F8766D","#00BFC4")) +
    # guides(colour = guide_legend(title.position="top", title.hjust = 0.5)) +
    coord_flip() +
    theme(
      axis.title.x = element_text(size=rel(1.5), color='black'),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.2), hjust=1, color='black'),
      axis.text.x = element_text(size=rel(1.5), color='black'),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(),
      legend.position='top',
      # legend.title=element_text(size=rel(1.5), color="black"),
      legend.title=element_blank(),
      legend.text=element_text(size=rel(1.3), color="black"),
      legend.key=element_rect(fill='transparent'),
      panel.background = element_blank(),
      aspect.ratio = .7
      )
  
  if (sign=="negative") p <- p + scale_x_discrete(position = "top")

  # If absolute values are used, add the corresponding signs to the plot
  if (abs) p <- p +  ylim(0,max(W$value)+0.1)+ geom_text(label=W$sign,y=max(W$value)+0.1, size=10)

  if(abs & scale) p <-  p + ylab(paste("Absolute loading on factor", factor))  
  else if(abs & !scale) p <- p + ylab(paste("Absolute loading on factor", factor))
  else if(!abs & scale) p <- p + ylab(paste("Loading on factor", factor))
  else p <- p + ylab(paste("Loading on factor", factor))
  return(p)
  
}


#' @title Plot correlation matrix between the columns of the weights (1 column = 1 factor) for a given view
#' @name plot_weight_cor
#' @description Function to plot the correlation matrix between the columns of the weights
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param method a character indicating the type of correlation coefficient to be computed: pearson (default), kendall, or spearman.
#' @param ... arguments passed to \code{\link[corrplot]{corrplot}}
#' @details This method plots the correlation matrix between the weights columns. \cr 
#' The model encourages the factors to be uncorrelated, so this function usually yields a diagonal correlation matrix. \cr 
#' However, it is not a hard constraint such as in Principal Component Analysis and correlations between factors can occur, 
#' particularly with large number factors. \cr
#' Generally, correlated factors are redundant and should be avoided, as they make interpretation harder. Therefore, 
#' if you have too many correlated factors we suggest you try reducing the number of factors.
#' @return Returns a symmetric matrix with the correlation coefficient between every pair of factors.
#' @importFrom corrplot corrplot
#' @export
plot_weight_cor <- function(object, view, method = "pearson", ...) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  
  # Fetch factors
  W <- get_weights(object)[[view]]
  
  # Compute and plot correlation
  r <- abs(cor(x=W, y=W, method=method, use = "complete.obs"))
  p <- corrplot(r, tl.col="black", ...)
  
  return(r)
}
