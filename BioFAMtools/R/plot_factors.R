
###########################################
## Functions to visualise latent factors ##
###########################################

#' @title Beeswarm plots of factor values
#' @name plot_factor
#' @description Beeswarm plot of the latent factor values.
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param factors character vector with the factor name(s), or numeric vector with the index of the factor(s) to use. 
#' Default is 'all'
#' @param color_by specifies groups or values (discrete or continuous) used to color the samples
#' This can be either: 
#' the string "group" : in this case, the plot will color samples with respect to the groups they belong to
#' a character giving the name of a feature, 
#' a character giving the same of a covariate (only if using \code{\link{MultiAssayExperiment}} as input), 
#' or a vector of the same length as the number of samples specifying discrete groups or continuous numeric values.
#' or a dataframe with two columns : "sample" (names of the samples) and "color" (values for each sample)
#' @param group_by specifies groups used to separate the samples : one plot per group 
#' This can be either: 
#' the string "group" : in this case, the plot will separate samples with respect to the groups they belong to
#' a character giving the name of a feature, 
#' a character giving the same of a covariate (only if using \code{\link{MultiAssayExperiment}} as input), 
#' or a vector of the same length as the number of samples specifying discrete groups.
#' or a dataframe with two columns : "sample" (names of the samples) and "group" (groups for each sample)
#' @param color_name name for color legend (usually only used if color_by is not a character itself)
#' @param show_missing logical indicating whether to remove samples for which \code{shape_by} or \code{color_by} is missing.
#' @details One of the main steps for the annotation of factors is to visualise and color them using known covariates or phenotypic data. \cr
#' This function generates a Beeswarm plot of the sample values in a given latent factor. \cr
#' Similar functions are \code{\link{plot_factor_scatter}} for doing scatter plots.
#' @return Returns a \code{ggplot2} object
#' @import ggplot2
#' @import ggbeeswarm
#' @import grDevices
#' @import RColorBrewer
#' @export
plot_factors <- function(object, factors = "all", group_by = "group", add_dots = TRUE, add_violin = TRUE, show_missing = TRUE, dot_size = 1,
                                 color_by = "group", color_name = "", shape_by = NULL, shape_name = "", 
                                 dots_alpha = 1.0, legend = TRUE,
                                 violin_alpha = 0.5, color_violin=TRUE,
                                 rasterize = FALSE, dodge = FALSE) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  
  # Get factor values
  if ((length(factors)==1) && (factors == "all")) {
    factors <- factors_names(object)
  } else if (is.numeric(factors)) {
    factors <- factors_names(object)[unique(factors)]
  } else { 
    stopifnot(all(factors %in% factors_names(object)))
  }
  Z <- get_factors(object, factors=factors, as.data.frame=T)
  Z$factor <- factor(Z$factor, levels=factors)
  
  # Set group/color/shape
  group_by <- .set_groupby(object, group_by)
  color_by <- .set_colorby(object, color_by)
  shape_by <- .set_shapeby(object, shape_by)
  
  # Remove samples with missing values
  Z <- Z[complete.cases(Z),]
  
  # Merge factor values with group/color/shape information
  df <- merge(Z, group_by, by="sample")
  df <- merge(df, color_by, by="sample")
  df <- merge(df, shape_by, by="sample")
  
  df$color_by <- as.character(df$color_by)
  df$group_by <- as.character(df$group_by)
  df$shape_by <- as.character(df$shape_by)
  
  # Remove samples with no sample metadata
  if (!show_missing) df <- filter(df, !is.na(color_by) & !is.na(shape_by))
  
  # Generate plot
  p <- ggplot(df, aes(x=group_by, y=value, color=color_by, shape=shape_by)) +
    facet_wrap(~factor, scales="free") +
    labs(x="", y="Factor value")

  # Add dots
  if (add_dots) {
    if (rasterize) {
      warning("geom_jitter is not available with rasterise==TRUE. We use instead ggrastr::geom_quasirandom_rast()")
      if (dodge) {
        p <- p + ggrastr::geom_quasirandom_rast(size=dot_size, position="dodge", dodge.width=1)
      } else {
        p <- p + ggrastr::geom_quasirandom_rast(size=dot_size)
      }
    } else {
      if (dodge) {
        p <- p + geom_jitter(size = dot_size, alpha = dots_alpha, position=position_jitterdodge(dodge.width=1))
      } else {
        p <- p + geom_jitter(size = dot_size, alpha = dots_alpha)
      }
    }
  }
  
  # Add violin plot
  if (add_violin) {
    if (color_violin) {
      tmp <- summarise(group_by(df, factor, color_by), n=n())
      if (min(tmp$n)==1) {
        warning("Warning: some 'color_by' groups have only one observation, violin plots are not displayed")
      } else {
        # violin_color <- ifelse(is.na(violin_color), color_by, violin_color)
        p <- p + geom_violin(aes(fill=color_by), alpha=violin_alpha, trim=T, scale="width", position=position_dodge(width=1))
        # if (add_dots) p <- p + scale_color_discrete(guide = FALSE)
      }
    } else {
      p <- p + geom_violin(color="black", fill="grey", alpha=violin_alpha, trim=T, scale="width")
    }
  }

  # If 'color_by' is numeric, define the default gradient
  if (is.numeric(df$color))
    p <- p + scale_color_gradientn(colors=colorRampPalette(rev(brewer.pal(n = 5, name = "RdYlBu")))(10)) 
  
  # Add legend for color
  if (length(unique(df$color_by))>1) { 
    p <- p + labs(color_by=color_name)
  } else { 
    p <- p + guides(fill=F, color=F) + 
      scale_color_manual(values="black") +
      scale_fill_manual(values="gray60")
  }
  
  # Add legend for shape
  if (length(unique(df$shape))>1) { 
    p <- p + labs(shape=shape_name)
  } else { 
    p <- p + guides(shape=F) 
  }

  # Use unified theme across the plots
  p <- p +
    theme_bw() + 
    theme(
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size=.1),
        panel.spacing = unit(1, "lines"),
        # axis.text.x = element_text(size=rel(1.0), color="black", angle=30, hjust=0, vjust=0.5),
        axis.text.x = element_text(size=rel(1.0), color="black"),
        axis.text.y = element_text(size=rel(1.0)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=rel(0.9), color="black"),
        # axis.line = element_line(color="black", size=1.0),
        axis.ticks = element_line(color = "black"),
        axis.title = element_text(size=rel(1.2)),
    )
  
  if (length(unique(df$factor))>1) {
    p <- p + scale_y_continuous(breaks=NULL)
  }
  
  # Add legend
  if (legend) {
    p <- p + theme(
      legend.title = element_text(size=rel(1.1), hjust=0.5, color="black"),
      legend.text = element_text(size=rel(1.1), hjust=0, color="black"),
      legend.position = "right", 
      legend.direction = "vertical",
      legend.key = element_blank()
    )
  } else {
    p <- p + theme(
      legend.position = "none"
    )
  }
  
  return(p)
}

#' @title Scatterplots of two factor values
#' @name plot_embeddings
#' @description Scatterplot of the values of two latent factors.
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param factors a vector of length two with the factors to plot. Factors can be specified either as a characters
#' using the factor names, or as numeric with the index of the factors
#' @param color_by specifies groups or values used to color the samples. 
#' This can be either 
#' the string "group" : in this case, the plot will color samples with respect to the groups they belong to
#' a character giving the name of a feature present in the training data, 
#' a character giving the same of a covariate (only if using \code{\link{MultiAssayExperiment}} as input), 
#' or a vector of the same length as the number of samples specifying discrete groups or continuous numeric values.
#' @param shape_by specifies groups or values used to shape the samples. 
#' This can be either
#' a character giving the name of a feature present in the training data, 
#' a character giving the same of a covariate (only if using \code{\link{MultiAssayExperiment}} as input), 
#' or a vector of the same length as the number of samples specifying discrete groups.
#' @param color_name name for color legend (usually only used if color_by is not a character itself)
#' @param shape_name name for shape legend (usually only used if shape_by is not a character itself)
#' @param show_missing logical indicating whether to include samples for which \code{shape_by} or \code{color_by} is missing
#' @param return_data logical indicating whether to return the data frame to plot instead of plotting
#' @details One of the first steps for the annotation of factors is to visualise and group/color them using known covariates such as phenotypic or clinical data.
#' This method generates a single scatterplot for the combination of two latent factors.
#' Similar function is
#' \code{\link{plot_factors}} for doing Beeswarm plots for factors.
#' @return Returns a \code{ggplot2} object
#' @import ggplot2
#' @export
plot_embeddings <- function(object, factors = c(1,2), show_missing = TRUE,
                            color_by = NULL, shape_by = NULL, color_name = NULL, shape_name = NULL,
                            dot_size = 1.5, alpha = 1, legend = TRUE, return_data = FALSE) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  
  # If plotting one or multiple factors, re-direct to other functions 
  if (length(unique(factors)) == 1) {
    .args <- as.list(match.call()[-1])
    return(do.call(plot_factors, .args))
  } else if (length(factors) > 2) {
    .args <- as.list(match.call()[-1])
    p <- do.call(.plot_multiple_factors, .args)
    return(p)
  }

  # Remember color_name and shape_name if not provided
  if (!is.null(color_by) && (length(color_by) == 1) && is.null(color_name))
    color_name <- color_by
  if (!is.null(shape_by) && (length(shape_by) == 1) && is.null(shape_name))
    shape_name <- shape_by
  
  # Get factors
  if ((length(factors) == 1) && (factors[1] == "all")) {
    factors <- factors_names(object)
  } else if (is.numeric(factors)) {
    factors <- factors_names(object)[unique(factors)]
  } else { 
    stopifnot(all(factors %in% factors_names(object)))
  }
  Z <- get_factors(object, factors=factors, as.data.frame=TRUE)
  # Z$factor <- factor(Z$factor, levels=factors)
  
  # Set color and shape
  color_by <- .set_colorby(object, color_by)
  shape_by <- .set_shapeby(object, shape_by)
  
  # Remove samples with missing values
  Z <- Z[complete.cases(Z),]
  
  # Merge factor values with color and shape information
  df <- merge(Z, color_by, by="sample")
  df <- merge(df, shape_by, by="sample")
  
  # Remove missing values
  if(!show_missing) df <- filter(df, !is.na(color_by) & !is.na(shape_by))
  
  # spread over factors
  df <- tidyr::spread(df, key="factor", value="value")
  df <- df[,c(colnames(df)[1:4], factors)]
  df <- magrittr::set_colnames(df, c(colnames(df)[1:4], "x", "y"))

  # Return data if requested instead of plotting
  if (return_data) return(df)
  
  # Generate plot
  p <- ggplot(df, aes(x=x, y=y)) + 
    geom_point(aes(color = color_by, shape = shape_by), size=dot_size, alpha=alpha) +
    # ggrastr::geom_point_rast(aes(color = color_by, shape = shape_by)) +
    labs(x=factors[1], y=factors[2]) +
    theme_classic() +
    theme(
      axis.text = element_text(size = rel(1.0), color = "black"), 
      axis.title = element_text(size = rel(1.3), color = "black"), 
      axis.line = element_line(color = "black", size = 0.5), 
      axis.ticks = element_line(color = "black", size = 0.5)
    )
  
  # If color is numeric, define the default gradient
  if (is.numeric(df$color))
    p <- p + scale_color_gradientn(colors=colorRampPalette(rev(brewer.pal(n = 5, name = "RdYlBu")))(10)) 
  
  # Add legend for color
  if (length(unique(df$color))>1) { 
    p <- p + labs(color=color_name)
  } else { 
    p <- p + guides(color=F) + scale_color_manual(values="black")
  }
  
  # Add legend for shape
  if (length(unique(df$shape))>1) { 
    p <- p + labs(shape=shape_name)
  } else { 
    p <- p + guides(shape=F) 
  }
  
  if (legend) {
    p <- p + theme(
      legend.key = element_rect(fill = "white"),
      legend.text = element_text(size=rel(1.2)),
      legend.title = element_text(size=rel(1.2))
    )
  } else {
      p <- p + theme(legend.position = "none")
  }
  
  return(p)
}


  
# Plot multiple factors as pairwise scatterplots
.plot_multiple_factors <- function(object, factors = "all", show_missing = TRUE, dot_size = 1,
                                   color_by = NULL, color_name = "", shape_by = NULL, shape_name = "", legend = TRUE) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  
  # Get factors
  if (is.numeric(factors)) {
    factors <- factors_names(object)[factors] 
  } else if (paste0(factors, collapse="") == "all") { 
    factors <- factors_names(object) 
  } else {
    stopifnot(all(factors %in% factors_names(object)))  
  }
  
  # Collect relevant data
  Z <- get_factors(object, factors=factors, as.data.frame=TRUE)

  # Set color and shape
  color_by <- .set_colorby(object, color_by)
  shape_by <- .set_shapeby(object, shape_by)
  
  # Remove samples with missing factor values
  Z <- Z[complete.cases(Z),]
  
  # Merge factor values with color and shape information
  df <- merge(Z, color_by, by="sample")
  df <- merge(df, shape_by, by="sample")
  
  # Remove missing values
  if(!show_missing) df <- filter(df, !is.na(color_by) & !is.na(shape_by))

  # Spread over factors
  df <- tidyr::spread(df, key="factor", value="value")
  
  # Prepare the legend
  p <- ggplot(df, aes_string(x=factors[1], y=factors[2], color="color_by", shape="shape_by")) +
    geom_point() +
    theme(
      legend.key = element_rect(fill = "white"),
      legend.text = element_text(size=rel(1.2)),
      legend.title = element_text(size=rel(1.2))
    )
  if (length(unique(df$color))>1) { p <- p + labs(color=color_name) } else { p <- p + guides(color=F) + scale_color_manual(values="black") }
  if (is.numeric(df$color)) p <- p + scale_color_gradientn(colors=colorRampPalette(rev(brewer.pal(n=5, name="RdYlBu")))(10)) 
  if (length(unique(df$shape))>1) { p <- p + labs(shape=shape_name) } else { p <- p + guides(shape = FALSE) }
  if (length(unique(df$color))>1 | length(unique(df$shape))>1) { legend <- GGally::grab_legend(p) } else { legend <- NULL }
  
  
  # Generate plot
  p <- GGally::ggpairs(df, 
    columns = factors,
    lower = list(continuous=GGally::wrap("points", size=dot_size)), 
    diag = list(continuous='densityDiag'), 
    upper = list(continuous=GGally::wrap("points", size=dot_size)), 
    mapping = aes(color=color_by, shape=shape_by), 
    title = "", 
    legend = legend
    )
  p <- p + theme_bw() + theme(
    # axis.line = element_line(color="black", size=rel(1.0)),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )
  
  if (!legend)
    p <- p + theme(legend.position = "none")
  
  return(p)
}
  


#' @title Plot correlation matrix between latent factors
#' @name plot_factor_cor
#' @description Function to plot the correlation matrix between the latent factors.
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param method a character indicating the type of correlation coefficient to be computed: pearson (default), kendall, or spearman.
#' @param ... arguments passed to \code{\link[corrplot]{corrplot}}
#' @details This method plots the correlation matrix between the latent factors. \cr 
#' The model encourages the factors to be uncorrelated, so this function usually yields a diagonal correlation matrix. \cr 
#' However, it is not a hard constraint such as in Principal Component Analysis and correlations between factors can occur, 
#' particularly with large number factors. \cr
#' Generally, correlated factors are redundant and should be avoided, as they make interpretation harder. Therefore, 
#' if you have too many correlated factors we suggest you try reducing the number of factors.
#' @return Returns a symmetric matrix with the correlation coefficient between every pair of factors.
#' @import corrplot
#' @export
plot_factor_cor <- function(object, method = "pearson", ...) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  
  # Fetch factors
  Z <- get_factors(object)
  
  # Compute and plot correlation
  r <- abs(cor(x=do.call(rbind, Z), y=do.call(rbind, Z), method=method, use = "complete.obs"))
  p <- corrplot(r, tl.col="black", ...)
  
  return(r)
}


# (Hidden) function to define the group
.set_groupby <- function(object, group_by) {
  
  # Option 0: no group
  if (is.null(group_by)) {
    group_by <- rep("1",sum(object@dimensions[["N"]]))
    
    # Option 1: by default group
  } else if (group_by[1] == "group") {
    group_by = c()
    for (group in names(samples_names(object))) {
      group_by <- c(group_by,rep(group,length(samples_names(object)[[group]])))
    }
    group_by = factor(group_by, levels=groups_names(object))
    
    # Option 2: by a metadata column in object@samples$metadata
  } else if ((length(group_by) == 1) && is.character(group_by) & (group_by[1] %in% colnames(samples_metadata(object)))) {
      group_by <- samples_metadata(object)[,group_by]

    # Option 3: input is a data.frame with columns (sample,group)
  } else if (is(group_by,"data.frame")) {
    stopifnot(all(colnames(group_by) %in% c("sample","group")))
    stopifnot(all(unique(group_by$sample) %in% unlist(samples_names(model))))
    
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
      sample = unlist(samples_names(object)),
      group_by = group_by
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
    color_by <- samples_groups(object)$group_name
    
  # Option 2: by a feature present in the training data    
  } else if ((length(color_by) == 1) && is.character(color_by) && (color_by[1] %in% unlist(features_names(object)))) {
      training_data <- lapply(get_training_data(object), function(l) Reduce(cbind, l))
      features_names <- lapply(training_data, rownames)
      viewidx <- which(sapply(features_names, function(vnm) color_by %in% vnm))
      color_by <- training_data[[viewidx]][color_by,]
    
  # Option 3: by a metadata column in object@samples$metadata
  } else if ((length(color_by) == 1) && is.character(color_by) & (color_by[1] %in% colnames(samples_metadata(object)))) {
      color_by <- samples_metadata(object)[,color_by]
        
  # Option 4: input is a data.frame with columns (sample, color)
  } else if (is(color_by, "data.frame")) {
    stopifnot(all(colnames(color_by) %in% c("sample", "color")))
    stopifnot(all(unique(color_by$sample) %in% unlist(samples_names(model))))
    
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
      color_by = color_by
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
    
  # Option 2: input is a data.frame with columns (sample,color)
  } else if (is(shape_by,"data.frame")) {
    stopifnot(all(colnames(shape_by) %in% c("sample","color")))
    stopifnot(all(unique(shape_by$sample) %in% unlist(samples_names(model))))

    # Option 3: by a metadata column in object@samples$metadata
  } else if ((length(shape_by) == 1) && is.character(shape_by) & (shape_by %in% colnames(samples_metadata(object)))) {
      shape_by <- samples_metadata(object)[,shape_by]
    
  # Option 4: shape_by is a vector of length N
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
      shape_by = as.factor(shape_by)
    )
  }

  return(df)
}
