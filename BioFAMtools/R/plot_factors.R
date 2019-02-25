
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
#' @export
plot_factor <- function(object, factors = "all", group_by = "group", add_dots=TRUE, add_violin=TRUE, show_missing = TRUE, dot_size = 1,
                                 color_by = NULL, color_name = "", shape_by = NULL, alpha=0.25, shape_name = "", rasterize = FALSE, dodge=FALSE) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")
  
  # Get factor values
  Z <- get_factors(object, factors=factors, as.data.frame=T)
  Z$factor <- as.factor(Z$factor)
  
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
  
  # Remove samples with no sample metadata
  if (!show_missing) df <- filter(df, !is.na(color_by) & !is.na(shape_by))
  
  # Generate plot
  p <- ggplot(df, aes(x=group_by, y=value, color=color_by, shape=shape_by)) +
    facet_wrap(~factor, scales="free") +
    ylab("Factor value") + xlab("") +
    theme(
      axis.text.y = element_text(size = rel(1.5), color = "black"),
      axis.title.y = element_text(size = rel(1.5), color = "black"),
      axis.line = element_line(color = "black", size = 0.4),
      axis.ticks.length = unit(0.25,"cm"),
      axis.ticks = element_line(color = "black"),
      panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      panel.background = element_blank(),
      legend.title=element_text(size=20, hjust=0.5, color="black"),
      legend.text=element_text(size=18, hjust=0.5, color="black"),
      legend.position = "right", 
      legend.direction = "vertical",
      legend.key = element_blank()
    ) 
  
  # Add dots
  if (add_dots) {
    if (rasterize) {
      if (dodge) {
        p <- p + ggrastr::geom_quasirandom_rast(size=dot_size, position="dodge", dodge.width=1)
      } else {
        p <- p + ggrastr::geom_quasirandom_rast(size=dot_size)
      }
    } else {
      if (dodge) {
        p <- p + ggbeeswarm::geom_quasirandom(size=dot_size, position="dodge", dodge.width=1)
      } else {
        p <- p + ggbeeswarm::geom_quasirandom(size=dot_size)
        
      }
    }
  }
  
  # Add violin plot
  if (add_violin) {
    if (dodge) {
      if (min(table(df$color_by))==1) {
        warning("Warning: some 'color_by' groups have only one observation, violin plots are not displayed")
      } else {
        p <- p + geom_violin(aes(fill=color_by), color="black", alpha=alpha, trim=F, scale="width", position=position_dodge(width = 1))
        if (add_dots) p <- p + scale_fill_discrete(guide = FALSE)
      }
    } else {
      p <- p + geom_violin(color="black", fill="grey", alpha=alpha, trim=F, scale="width")
    }
  }
  
  # If 'color_by' is numeric, define the default gradient
  if (is.numeric(df$color))
    p <- p + scale_color_gradientn(colors=colorRampPalette(rev(brewer.pal(n = 5, name = "RdYlBu")))(10)) 
  
  # Add legend for color
  if (length(unique(df$color))>1) { 
    p <- p + labs(color=color_name) + theme(legend.text.align = 0)
  } else { 
    p <- p + guides(color = FALSE) 
  }
  
  # Add legend for shape
  if (length(unique(df$shape))>1) { 
    p <- p + labs(shape=shape_name) + theme(legend.text.align = 0)
  } else { 
    p <- p + guides(shape = FALSE) 
  }
  
  return(p)
}

#' @title Scatterplots of two factor values
#' @name plot_factors
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
#' @details One of the first steps for the annotation of factors is to visualise and group/color them using known covariates such as phenotypic or clinical data.
#' This method generates a single scatterplot for the combination of two latent factors.
#' Similar functions are \code{\link{plot_factor_scatters}} for doing multiple scatter plots and 
#' \code{\link{plot_factor_beeswarm}} for doing Beeswarm plots for single factors.
#' @return Returns a \code{ggplot2} object
#' @import ggplot2
#' @export
plot_factors <- function(object, factors, show_missing = TRUE, dot_size=1,
                         color_by = NULL, shape_by = NULL, color_name="", shape_name="") {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  
  # If plotting one or multiple factors, re-direct to other functions 
  if (length(factors)==1) {
    .args <- as.list(match.call()[-1])
    do.call(plot_factor_beeswarm, .args)   
  } else if (length(factors)>2) {
    .args <- as.list(match.call()[-1])
    do.call(.plot_multiple_factors, .args)
  }
  
  # Get factors
  if (is.numeric(factors)) {
    factors <- factors_names(object)[factors]
  } else { 
    stopifnot(all(factors %in% factors_names(object)))
  }
  Z <- get_factors(object, factors=factors, as.data.frame=TRUE)
  Z$factor <- as.factor(Z$factor)
  
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
  df <- magrittr::set_colnames(df,c(colnames(df)[1:4],"x","y"))
  
  # Generate plot  
  p <- ggplot(df, aes(x=x, y=y)) + 
    geom_point(aes(color = color_by, shape = shape_by)) +
    # ggrastr::geom_point_rast(aes(color = color_by, shape = shape_by)) +
    xlab(factors[1]) + ylab(factors[2]) +
    theme(
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
  
  # If color is numeric, define the default gradient
  if (is.numeric(df$color))
    p <- p + scale_color_gradientn(colors=colorRampPalette(rev(brewer.pal(n = 5, name = "RdYlBu")))(10)) 
  
  # Add legend for color
  if (length(unique(df$color))>1) { 
    p <- p + labs(color=color_name) + theme(legend.text.align = 0)
  } else { 
    p <- p + guides(color = FALSE) 
  }
  
  # Add legend for shape
  if (length(unique(df$shape))>1) { 
    p <- p + labs(shape=shape_name) + theme(legend.text.align = 0)
  } else { 
    p <- p + guides(shape = FALSE) 
  }
  
  return(p)
}
  
  

  
# Plot multiple factors as pairwise scatterplots
.plot_multiple_factors <- function(object, factors = "all", show_missing=TRUE, dot_size=1,
                                   color_by=NULL, color_name="", shape_by=NULL, shape_name="") {
  
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
    geom_point(size = dot_size)
  if (length(unique(df$color))>1) { p <- p + labs(color=color_name) } else { p <- p + guides(color = FALSE) }
  if (is.numeric(df$color)) p <- p + scale_color_gradientn(colors=colorRampPalette(rev(brewer.pal(n = 5, name = "RdYlBu")))(10)) 
  if (length(unique(df$shape))>1) { p <- p + labs(shape=shape_name) } else { p <- p + guides(shape = FALSE) }
  if (length(unique(df$color))>1 | length(unique(df$shape))>1) { legend <- GGally::grab_legend(p) } else { legend <- NULL }
  
  
  # Generate plot
  p <- GGally::ggpairs(df, 
    columns = factors,
    lower = list(continuous="points"), diag=list(continuous='blankDiag'), upper=list(continuous='points'),
    mapping = aes(color=color_by, shape=shape_by), 
    title = "", 
    legend = legend
    )

  p <- p + theme_bw() +
    theme(
      plot.title = element_text(size = 16, hjust=0.5, color="black"),
      axis.title = element_text(size = 10, color="black"),
      axis.text = element_text(size = 9, color="black"),
      legend.position = "right",
      legend.direction = "vertical"
    )
  
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
    group_by <- rep(TRUE,sum(object@dimensions[["N"]]))
    
    # Option 1: by default group
  } else if (group_by[1] == "group") {
    group_by = c()
    for (group in names(samples_names(object))) {
      group_by <- c(group_by,rep(group,length(samples_names(object)[[group]])))
    }
    group_by = factor(group_by, levels=groups_names(object))
    
    # Option 2: input is a data.frame with columns (sample,group)
  } else if (is(group_by,"data.frame")) {
    stopifnot(all(colnames(group_by) %in% c("sample","group")))
    stopifnot(all(unique(group_by$sample) %in% unlist(samples_names(model))))
    
    # Option 3: group_by is a vector of length N
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
    color_by <- rep(TRUE,sum(object@dimensions[["N"]]))
    
  # Option 1: by default group
  } else if (color_by[1] == "group") {
    color_by = c()
    for (group in names(samples_names(object))){
      color_by <- c(color_by,rep(group,length(samples_names(object)[[group]])))
    }
    
  # Option 2: by a feature present in the training data    
  } else if (length(color_by) == 1 & is.character(color_by) & color_by[1]%in%unlist(features_names(object))) {
      training_data <- lapply(get_training_data(object), function(l) Reduce(cbind, l))
      features_names <- lapply(training_data, rownames)
      viewidx <- which(sapply(features_names, function(vnm) color_by %in% vnm))
      color_by <- training_data[[viewidx]][color_by,]
        
  # Option 3: input is a data.frame with columns (sample,color)
  } else if (is(color_by,"data.frame")) {
    stopifnot(all(colnames(color_by) %in% c("sample","color")))
    stopifnot(all(unique(color_by$sample) %in% unlist(samples_names(model))))
    
  # Option 4: color_by is a vector of length N
  } else if (length(color_by) > 1) {
    stopifnot(length(color_by) == sum(object@dimensions[["N"]]))
    
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
    shape_by <- rep(TRUE,sum(object@dimensions[["N"]]))
    
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
    
  # Option 3: shape_by is a vector of length N
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
