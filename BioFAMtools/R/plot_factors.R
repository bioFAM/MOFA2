
###########################################
## Functions to visualise latent factors ##
###########################################


#' @title Plot histogram of latent factor values
#' @name plot_factor_hist
#' @description Plot a histogram of latent factor values.
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param factor character vector with the factor name or numeric vector with the index of the factor.
#' @param group_by specifies groups used to color the samples of the histogram. 
#' This can be either: 
#' a character giving the name of a feature,
#' the name of a covariate (only if using a \code{\link{MultiAssayExperiment}} as input), 
#' or a vector of the same length as the number of samples.
#' @param group_names names for the groups.
#' @param alpha transparency parameter. 
#' Default is 0.5
#' @param binwidth binwidth for histogram. Default is \code{NULL}, which uses \code{ggplot2} default calculation.
#' @param show_missing boolean indicating whether to remove sample for which \code{group_by} is missing (default is FALSE)
#' @details One of the first steps for the annotation of factors is to visualise and color them using known covariates such as phenotypic or clinical data. \cr
#' This method generates a histogram of the sample values in a given latent factor. \cr
#' Similar functions are \code{\link{plot_factor_scatter}} for doing scatter plots between pairs of factors 
#' and \code{\link{plot_factor_beeswarm}} for doing Beeswarm plots of single factors.
#' @return Returns a \code{ggplot2} object
#' @import ggplot2
#' @export
plot_factor_hist <- function(object, factor, group_by = NULL, group_names = "", alpha = 0.5, binwidth = NULL, show_missing = FALSE, ...) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  if(!factor %in% factors_names(object)) { stop("factor not recognised") }
  
  # Collect relevant data
  N <- sum(object@dimensions[["N"]])
  Z <- get_factors(object, factors = factor, as.data.frame = TRUE)
  
  # get groups
  group_legend <- T
  if (!is.null(group_by)) {
    
    # It is the name of a covariate or a feature in the training_data
    if (length(group_by) == 1 & is.character(group_by)) {
      if(group_names=="") group_names <- group_by
      training_data <- get_training_data(object)
      features_names <- lapply(training_data(object), rownames)
      if(group_by %in% Reduce(union,features_names)) {
        viewidx <- which(sapply(features_names, function(vnm) group_by %in% vnm))
        group_by <- training_data[[viewidx]][group_by,]
      } else if(class(object@input_data) == "MultiAssayExperiment"){
        group_by <- get_covariates(object, group_by)
      }
      else stop("'group_by' was specified but it was not recognised, please read the documentation")
      
    # It is a vector of length N
    } else if (length(group_by) > 1) {
      stopifnot(length(group_by) == N)
      
    # It is not recognised
    } else {
      stop("'group_by' was specified but it was not recognised, please read the documentation")
    }
    
  } else {
    group_by <- rep(TRUE, N)
    group_legend <- F
  }
  
  names(group_by) <- unname(do.call(c, samples_names(object)))
  Z$group_by <- group_by[Z$sample]

  # Remove missing samples
  if(!show_missing) Z <- Z[!is.na(Z$group_by),]
  Z$group_by <- as.factor(Z$group_by)
  
  # Generate plot
  p <- ggplot(Z, aes_string(x="value", group="group_by")) + 
    geom_histogram(aes(fill=group_by), alpha=alpha, binwidth=binwidth, position="identity", ...) + 
    scale_y_continuous(expand=c(0,0)) +
    guides(fill=guide_legend(title=group_names)) +
    theme(plot.margin = margin(40,40,20,20), 
          axis.text = element_text(size=rel(1.3), color = "black"), 
          axis.title.y = element_text(size=rel(1.5), margin=margin(0,15,0,0)), 
          axis.title.x = element_text(size=rel(1.5), margin=margin(15,0,0,0)), 
          axis.line = element_line(color="black", size=rel(1.0)),
          # axis.ticks = element_line(color="black", size=0.5),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.key = element_rect(fill = "white")
    )
  
  if (!group_legend) { p <- p + guides(fill = FALSE) }
  
  return(p)
}


#' @title Beeswarm plot of latent factors
#' @name plot_factor_beeswarm
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
#' @param name_color name for color legend (usually only used if color_by is not a character itself)
#' @param show_missing logical indicating whether to remove samples for which \code{shape_by} or \code{color_by} is missing.
#' @details One of the main steps for the annotation of factors is to visualise and color them using known covariates or phenotypic data. \cr
#' This function generates a Beeswarm plot of the sample values in a given latent factor. \cr
#' Similar functions are \code{\link{plot_factor_scatter}} for doing scatter plots and \code{\link{plot_factor_hist}} for doing histogram plots
#' @return Returns a \code{ggplot2} object
#' @import ggplot2
#' @import ggbeeswarm
#' @import grDevices
#' @export
plot_factor_beeswarm <- function(object, factors = "all", group_by = NULL, color_by = NULL, name_color = "", show_missing = FALSE) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")

  # Collect relevant data
  N <- sum(object@dimensions[["N"]])
  Z <- get_factors(object, factors=factors, include_intercept=FALSE, as.data.frame=T)
  Z$factor <- as.factor(Z$factor)
  
  if (!is.null(group_by)){
    if (group_by == "group"){
      group_by = c()
      for (group in names(samples_names(object))){
        group_by <- c(group_by,rep(group,length(samples_names(object)[[group]])))
      }
    }
  }
  
  # Set color
  if (!is.null(group_by)) {
    if (class(group_by)!="data.frame"){
      # It is the name of a covariate or a feature in the training_data
      if (length(group_by) == 1 & is.character(group_by)) {
        if(name_color=="") name_color <- group_by
        training_data <- lapply(get_training_data(object), function(l) Reduce(cbind, l))
        features_names <- lapply(training_data, rownames)
        if(group_by %in% Reduce(union,features_names)) {
          viewidx <- which(sapply(features_names, function(vnm) group_by %in% vnm))
          group_by <- training_data[[viewidx]][group_by,]
        } else if(class(object@input_data) == "MultiAssayExperiment"){
          group_by <- get_covariates(object, group_by)
        }
        else stop("'group_by' was specified but it was not recognised, please read the documentation")
        # It is a vector of length N
      } else if (length(group_by) > 1) {
        stopifnot(length(group_by) == N)
        # group_by <- as.factor(group_by)
      } else {
        stop("'group_by' was specified but it was not recognised, please read the documentation")
      }
    }
  } else {
    group_by <- rep(TRUE,N)
  }
  
  if (!is.null(color_by)){
    if (color_by == "group"){
      color_by = c()
      for (group in names(samples_names(object))){
        color_by <- c(color_by,rep(group,length(samples_names(object)[[group]])))
      }
    }
  }
  
  # Set color
  color_legend <- T
  if (!is.null(color_by)) {
    if (class(color_by)!="data.frame"){
      # It is the name of a covariate or a feature in the training_data
      if (length(color_by) == 1 & is.character(color_by)) {
        if(name_color=="") name_color <- color_by
        training_data <- lapply(get_training_data(object), function(l) Reduce(cbind, l))
        features_names <- lapply(training_data, rownames)
        if(color_by %in% Reduce(union,features_names)) {
          viewidx <- which(sapply(features_names, function(vnm) color_by %in% vnm))
          color_by <- training_data[[viewidx]][color_by,]
        } else if(class(object@input_data) == "MultiAssayExperiment"){
          color_by <- get_covariates(object, color_by)
        }
        else stop("'color_by' was specified but it was not recognised, please read the documentation")
        # It is a vector of length N
      } else if (length(color_by) > 1) {
        stopifnot(length(color_by) == N)
        # color_by <- as.factor(color_by)
      } else {
        stop("'color_by' was specified but it was not recognised, please read the documentation")
      }
    }
  } else {
    color_by <- rep(TRUE,N)
    color_legend <- F
  }
  
  # Remove samples with missing values
  if (!show_missing) {
    Z <- Z[complete.cases(Z),]
    # Z <- Z[!is.na(Z$value),]
  }
  
  if (class(color_by)!="data.frame"){
    df <- data.frame("sample"=Reduce(c,samples_names(object)),"color_name"= color_by)
    df <- left_join(Z, df, by="sample")
  }
  else{
    df <- merge(Z,color_by,by="sample")
  }
  
  if (class(group_by)!="data.frame"){
    df2 <- data.frame("sample"=Reduce(c,samples_names(object)),"group_name"= group_by)
    df <- left_join(df, df2, by="sample")
  }
  else{
    df <- merge(df,group_by,by="sample")
  }
  
  p <- ggplot(df, aes(x=group_name, y=value)) +
     #geom_violin(aes(color=grouped_by)) +  
     ggbeeswarm::geom_quasirandom(aes(color=color_name)) 
     #geom_boxplot(width=0.1) +
     #scale_x_continuous(breaks=NULL)

  # Generate plot
  
  p <- p +
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
    ) + facet_wrap(~factor, scales="free")
  
  # If color_by is numeric, define the default gradient
  # if (is.numeric(color_by)) { p <- p + scale_color_gradientn(colors=terrain.colors(10)) }
  if (is.numeric(color_by)) { p <- p + scale_color_gradientn(colors=colorRampPalette(rev(brewer.pal(n = 5, name = "RdYlBu")))(10)) }
  
  # Add legend
  if (color_legend) { p <- p + labs(color=name_color) } else { p <- p + guides(color = FALSE) }
  
  return(p)
}

#' @title Scatterplot of two latent factors
#' @name plot_factor_scatter
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
#' @param name_color name for color legend (usually only used if color_by is not a character itself)
#' @param name_shape name for shape legend (usually only used if shape_by is not a character itself)
#' @param show_missing logical indicating whether to include samples for which \code{shape_by} or \code{color_by} is missing
#' @details One of the first steps for the annotation of factors is to visualise and group/color them using known covariates such as phenotypic or clinical data.
#' This method generates a single scatterplot for the combination of two latent factors.
#' Similar functions are \code{\link{plot_factor_scatters}} for doing multiple scatter plots and 
#' \code{\link{plot_factor_beeswarm}} for doing Beeswarm plots for single factors.
#' @return Returns a \code{ggplot2} object
#' @import ggplot2
#' @export
plot_factor_scatter <- function(object, factors, groups = "all", color_by = NULL, shape_by = NULL, name_color="",
                         name_shape="", show_missing = TRUE) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  stopifnot(length(factors)==2)
  stopifnot(all(factors %in% factors_names(object)))
  
  groups <- .check_and_get_groups(object, groups)
  
  # Collect relevant data  
  N <- sum(object@dimensions[["N"]])
  Z <- get_factors(object, factors=factors, include_intercept=FALSE, as.data.frame=TRUE)
  
  if (!is.null(color_by)){
    if (color_by == "group"){
      color_by = c()
      for (group in names(samples_names(object))){
        color_by <- c(color_by,rep(group,length(samples_names(object)[[group]])))
      }
    }
  }
  
  # Set color
  color_legend <- T
  if (!is.null(color_by)) {
    # It is the name of a covariate or a feature in the training_data
    if (length(color_by) == 1 & is.character(color_by)) {
      if(name_color=="") name_color <- color_by
      training_data <- lapply(get_training_data(object), function(l) Reduce(cbind, l))
      features_names <- lapply(training_data, rownames)
      if(color_by %in% Reduce(union,features_names)) {
        viewidx <- which(sapply(features_names, function(vnm) color_by %in% vnm))
        color_by <- training_data[[viewidx]][color_by,]
      } else if(class(object@input_data) == "MultiAssayExperiment"){
        color_by <- get_covariates(object, color_by)
    }
    else stop("'color_by' was specified but it was not recognised, please read the documentation")
    # It is a vector of length N
    } else if (length(color_by) > 1) {
      stopifnot(length(color_by) == N)
      # color_by <- as.factor(color_by)
    } else {
      stop("'color_by' was specified but it was not recognised, please read the documentation")
    }
  } else {
    color_by <- rep(TRUE,N)
    color_legend <- F
  }

  # Set shape
  shape_legend <- T
  if (!is.null(shape_by)) {
    # It is the name of a covariate 
    if (length(shape_by) == 1 & is.character(shape_by)) {
      if(name_shape=="") name_shape <- shape_by
      training_data <- lapply(get_training_data(object), function(l) Reduce(cbind, l))
      features_names <- lapply(training_data, rownames)
      if(shape_by %in% Reduce(union,features_names)) {
        viewidx <- which(sapply(features_names, function(vnm) shape_by %in% vnm))
        shape_by <- training_data[[viewidx]][shape_by,]
      } else if(class(object@input_data) == "MultiAssayExperiment"){
        shape_by <- get_covariates(object, shape_by)
    }
    else stop("'shape_by' was specified but it was not recognised, please read the documentation")
    # It is a vector of length N
    # It is a vector of length N
    } else if (length(shape_by) > 1) {
      stopifnot(length(shape_by) == N)
    } else {
      stop("'shape_by' was specified but it was not recognised, please read the documentation")
    }
  } else {
    shape_by <- rep(TRUE,N)
    shape_legend <- F
  }
  
  df_aes <- data.frame(sample=Reduce(c,samples_names(object)), color_by = color_by, shape_by = shape_by) 
  df <- left_join(Z, df_aes, by="sample")
  
  # Remove missing values
  if(!show_missing) df <- filter(df, !is.na(color_by) & !is.na(shape_by))
  
  # Subset to relevant groups
  df <- filter(df, group %in% groups)
  
  #turn into factors
  df$shape_by[is.na(df$shape_by)] <- "NA"
  df$shape_by <- as.factor(df$shape_by)
  if(length(unique(df$color_by)) < 5) df$color_by <- as.factor(df$color_by)
  
  xlabel <- paste("Latent factor", factors[1])
  ylabel <- paste("Latent factor", factors[2])
  
  df_mat <- df %>% spread(key="factor", value="value")
  df_mat <-  set_colnames(df_mat,c(colnames(df_mat)[1:4], "x", "y"))
  
  p <- ggplot(df_mat, aes(x = x, y = y)) + 
      geom_point(aes(color = color_by, shape = shape_by)) + xlab(xlabel) + ylab(ylabel) +
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
  if (color_legend) { p <- p + labs(color = name_color) } else { p <- p + guides(color = FALSE) }
  if (shape_legend) { p <- p + labs(shape = name_shape) }  else { p <- p + guides(shape = FALSE) }
  return(p)
}
  
  
#' @title Pairwise scatterplots of multiple latent factors
#' @name plot_factor_scatters
#' @description Scatterplots of the sample values for pair-wise combinations of multiple latent factors.
#' @param object a \code{\link{BioFAModel}} object.
#' @param factors character vector with the factor name(s), or numeric vector with the index of the factor(s) to use. 
#' Default is 'all'
#' @param color_by specifies groups or values used to color the samples. 
#' This can be either: 
#' the string "group" : in this case, the plot will color samples with respect to the groups they belong to
#' a character giving the name of a feature present in the training data, 
#' a character giving the same of a covariate (only if using \code{\link{MultiAssayExperiment}} as input), 
#' or a vector of the same length as the number of samples specifying discrete groups or continuous numeric values.
#' @param shape_by specifies groups or values used to shape the samples. 
#' This can be either: 
#' a character giving the name of a feature present in the training data, 
#' a character giving the same of a covariate (only if using \code{\link{MultiAssayExperiment}} as input), 
#' or a vector of the same length as the number of samples specifying discrete groups.
#' @param name_color name for color legend (usually only used if color_by is not a character itself)
#' @param name_shape name for shape legend (usually only used if shape_by is not a character itself)
#' @param show_missing logical indicating whether to include samples for which \code{shape_by} or \code{color_by} is missing
#' @details One of the first steps for the annotation of factors is to visualise and group/color them using known covariates such as phenotypic or clinical data.
#' This method generates multiple scatterplots for pairwise combinations of several latent factors.
#' Similar functions are \code{\link{plot_factor_scatter}} for doing single scatter plots and 
#' \code{\link{plot_factor_beeswarm}} for doing Beeswarm plots for single factors.
#' @return \code{ggplot2} object
#' @import ggplot2
#' @export
plot_factor_scatters <- function(object, factors = "all", groups = "all", 
                                 show_missing=TRUE, color_by=NULL, name_color="",  
                                 shape_by=NULL, name_shape="") {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")

  groups <- .check_and_get_groups(object, groups)

  # Get factors
  if (paste0(factors, collapse="") == "all") { 
    factors <- factors_names(object) 
    # if(is.null(factors)) factors <- 1:ncol(Z) # old object are not compatible with factro names
  } else {
    stopifnot(all(factors %in% factors_names(object)))  
  }
  
  # Collect relevant data
  N <- sum(object@dimensions[["N"]])
  Z <- get_factors(object, factors=factors, include_intercept=FALSE, as.data.frame=TRUE)

  # Remove constant factors 
  tmp <- group_by(Z, factor) %>% mutate(var=var(value, na.rm = TRUE)) %>% ungroup()
  if (any(tmp$var==0)) {
    Z <- filter(Z, var>0)
    factors <-  unqiue(Z$factor)
  }
  
  if (!is.null(color_by)){
    if (color_by == "group"){
      color_by = c()
      for (group in names(samples_names(object))){
        color_by <- c(color_by,rep(group,length(samples_names(object)[[group]])))
      }
    }
  }
    
  # Set color
  color_legend <- T
  if (!is.null(color_by)) {
    # It is the name of a covariate or a feature in the training_data
    if (length(color_by) == 1 & is.character(color_by)) {
      if(name_color=="") name_color <- color_by
      training_data <- lapply(get_training_data(object), function(l) Reduce(cbind, l))
      features_names <- lapply(training_data, rownames)
      if(color_by %in% Reduce(union,features_names)) {
        viewidx <- which(sapply(features_names, function(vnm) color_by %in% vnm))
        color_by <- training_data[[viewidx]][color_by,]
      } else if(class(object@input_data) == "MultiAssayExperiment"){
        color_by <- get_covariates(object, color_by)
    }
    else stop("'color_by' was specified but it was not recognised, please read the documentation")
    # It is a vector of length N
    } else if (length(color_by) > 1) {
      stopifnot(length(color_by) == N)
      # color_by <- as.factor(color_by)
    } else {
      stop("'color_by' was specified but it was not recognised, please read the documentation")
    }
  } else {
    color_by <- rep(TRUE,N)
    color_legend <- F
  }

  # Set shape
  shape_legend <- T
  if (!is.null(shape_by)) {
    # It is the name of a covariate 
    if (length(shape_by) == 1 & is.character(shape_by)) {
      if(name_shape=="") name_shape <- shape_by
      training_data <- lapply(get_training_data(object), function(l) Reduce(cbind, l))
      features_names <- lapply(training_data, rownames)
      if (shape_by %in% Reduce(union,features_names)) {
        viewidx <- which(sapply(features_names, function(vnm) shape_by %in% vnm))
        shape_by <- training_data[[viewidx]][shape_by,]
      } else if(class(object@input_data) == "MultiAssayExperiment"){
        shape_by <- get_covariates(object, shape_by)
    }
    else stop("'shape_by' was specified but it was not recognised, please read the documentation")
    # It is a vector of length N
    # It is a vector of length N
    } else if (length(shape_by) > 1) {
      stopifnot(length(shape_by) == N)
    } else {
      stop("'shape_by' was specified but it was not recognised, please read the documentation")
    }
  } else {
    shape_by <- rep(TRUE,N)
    shape_legend <- F
  }
  
  df_aes <- data.frame(sample=Reduce(c,samples_names(object)), color_by = color_by, shape_by = shape_by) 
  df <- left_join(Z, df_aes, by="sample")

  # Remove missing values
  if(!show_missing) df <- filter(df, !is.na(color_by) & !is.na(shape_by))

  # Subset to relevant groups
  df <- filter(df, group %in% groups)
  
  #turn into factors
  df$shape_by[is.na(df$shape_by)] <- "NA"
  df$shape_by <- as.factor(df$shape_by)
  if(length(unique(df$color_by)) < 5) df$color_by <- as.factor(df$color_by)
  
  # Define title and legend of the plot
  main <- "" 
  p <- ggplot(df, aes_string(x=colnames(df)[1], y=colnames(df)[2], color="color_by", shape="shape_by")) + geom_point()
  if (color_legend | shape_legend) { 
    p <- p +
      theme(
        legend.title=element_text(size=15, hjust=0.5, color="black"),
        legend.position = "right", 
        legend.direction = "vertical",
        legend.key = element_blank()
      )
    
    # If color_by is numeric, define the default gradient
    if (is.numeric(df$color_by)) { p <- p + scale_color_gradientn(colors=terrain.colors(10)) }
    
    if (color_legend) { p <- p + labs(color = name_color) } else { p <- p + guides(color = FALSE) }
    if (shape_legend) { p <- p + labs(shape = name_shape) }  else { p <- p + guides(shape = FALSE) }
    # Extract the legend
    legend <- GGally::grab_legend(p)
  } else {
    legend <- NULL
  }
  
  # Generate plot
  df_mat <- df %>% spread(key="factor", value="value")
  p <- GGally::ggpairs(df_mat, columns = colnames(df_mat[,!colnames(df_mat) %in% c("sample","group","color_by","shape_by")]), 
                  lower=list(continuous="points"), diag=list(continuous='blankDiag'), upper=list(continuous='points'),
          mapping=aes(color=color_by, shape=shape_by), title=main, legend=legend) +
    theme_bw() +
    theme(plot.title = element_text(size = 16, hjust=0.5, color="black"), 
          axis.title = element_text(size = 10, color="black"), 
          axis.text = element_text(size = 9, color="black"),
          legend.position = "right", 
          legend.direction = "vertical"
          )
  
  # If color_by is numeric, define the default gradient
  if (is.numeric(df$color_by)) { 
    for(i in 1:p$nrow) {
      for(j in 1:p$ncol){
        p[i,j] <- p[i,j] + scale_color_gradientn(colors=terrain.colors(10)) 
      }
    }
  }
  
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
#' @importFrom corrplot corrplot
#' @export
plot_factor_cor <- function(object, method = "pearson", ...) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  
  # Fetch factors
  Z <- get_factors(object)
  
  # Remove intercept
  if(object@model_options$learn_intercept) Z <- lapply(Z, function(z) z[,-1])
  
  # Compute and plot correlation
  r <- abs(cor(x=do.call(rbind, Z), y=do.call(rbind, Z), method=method, use = "complete.obs"))
  p <- corrplot(r, tl.col="black", ...)
  
  return(r)
}

#' @title Plot the sparsity of factors
#' @name plot_sparsity_factors
#' @description Plot the sparsity of the factors (computed as 1-<Theta_k> for factor k)
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param threshold_variance_explained : do not plot the factors explaining less of the variance explained than the threshold for the corresponding views
#' @param show_variance_explained : indicate the variance explained by each of the factor in each view 
#' @details This method enables to visualize simultaneously the factors that are active in the different views (with the use of the threshold)
#' and their sparsity level. A sparse factor could be biologically meaningful, while a dense one is probably a technical source of variability.
#' @return Returns a \code{ggplot2} object
#' @import ggplot2
#' @export
plot_sparsity_factors <- function(object, threshold_variance_explained = 0.01, show_variance_explained = TRUE) {
  
  factors    <- factors_names(object)
  views      <- views_names(object)
  factor_val <- rep(factors_names(object), length(views))

  view_val  <- c()
  theta_val <- c()
  var_val   <- c()
  
  for (view in views) {
    view_val  <- c(view_val, rep(view, length(factors)))
    theta_val <- c(theta_val, object@expectations$ThetaW[[view]])
    var_val   <- c(var_val, calculate_variance_explained(object,only="views")$r2_per_factor[,view])
  }

  d <- data.frame(factor_val, view_val, 1-theta_val, var_val)
  colnames(d) <- c("factor", "view", "sparsity", "variance_explained")
  
  #sort factors by names (work only if factors are int encoded by strings)
  d$factor <- factor(d$factor, levels=seq(1,max(as.integer(levels(d$factor)))))
  
  d <- d[d$variance_explained > threshold_variance_explained,]
  if (show_variance_explained) {
    p = ggplot(d, aes(x=factor, y=sparsity, fill=view, alpha=variance_explained)) +
      geom_bar(position="dodge", stat="identity", aes(alpha=variance_explained))
  }
  else{
    p = ggplot(d, aes(x=factor, y=sparsity, fill=view)) +
      geom_bar(position="dodge", stat="identity") 
  }
  
  return(p)
}



