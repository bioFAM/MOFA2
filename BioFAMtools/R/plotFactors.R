
###########################################
## Functions to visualise latent factors ##
###########################################


#' @title Plot histogram of latent factor values
#' @name plotFactorHist
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
#' @param showMissing boolean indicating whether to remove sample for which \code{group_by} is missing (default is FALSE)
#' @details One of the first steps for the annotation of factors is to visualise and color them using known covariates such as phenotypic or clinical data. \cr
#' This method generates a histogram of the sample values in a given latent factor. \cr
#' Similar functions are \code{\link{plotFactorScatter}} for doing scatter plots between pairs of factors 
#' and \code{\link{plotFactorBeeswarm}} for doing Beeswarm plots of single factors.
#' @return Returns a \code{ggplot2} object
#' @import ggplot2
#' @export
plotFactorHist <- function(object, factor, group_by = NULL, group_names = "", alpha = 0.5, binwidth = NULL, showMissing = FALSE, ...) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  if(!factor %in% factorNames(object)) { stop("factor not recognised") }
  
  # Collect relevant data
  N <- sum(object@Dimensions[["N"]])
  Z <- getFactors(object, factors = factor, as.data.frame = TRUE)
  
  # get groups
  groupLegend <- T
  if (!is.null(group_by)) {
    
    # It is the name of a covariate or a feature in the TrainData
    if (length(group_by) == 1 & is.character(group_by)) {
      if(group_names=="") group_names <- group_by
      TrainData <- getTrainData(object)
      featureNames <- lapply(TrainData(object), rownames)
      if(group_by %in% Reduce(union,featureNames)) {
        viewidx <- which(sapply(featureNames, function(vnm) group_by %in% vnm))
        group_by <- TrainData[[viewidx]][group_by,]
      } else if(class(object@InputData) == "MultiAssayExperiment"){
        group_by <- getCovariates(object, group_by)
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
    groupLegend <- F
  }
  
  names(group_by) <- unname(do.call(c, sampleNames(object)))
  Z$group_by <- group_by[Z$sample]

  # Remove missing samples
  if(!showMissing) Z <- Z[!is.na(Z$group_by),]
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
  
  if (!groupLegend) { p <- p + guides(fill = FALSE) }
  
  return(p)
}


#' @title Beeswarm plot of latent factors
#' @name plotFactorBeeswarm
#' @description Beeswarm plot of the latent factor values.
#' @param object a trained \code{\link{BioFAModel}} object.
#' @param factors character vector with the factor name(s), or numeric vector with the index of the factor(s) to use. 
#' Default is 'all'
#' @param color_by specifies groups or values used to color the samples. 
#' This can be either: 
#' a character giving the name of a feature, 
#' a character giving the same of a covariate (only if using \code{\link{MultiAssayExperiment}} as input), 
#' or a vector of the same length as the number of samples specifying discrete groups or continuous numeric values.
#' @param name_color name for color legend (usually only used if color_by is not a character itself)
#' @param showMissing logical indicating whether to remove samples for which \code{shape_by} or \code{color_by} is missing.
#' @details One of the main steps for the annotation of factors is to visualise and color them using known covariates or phenotypic data. \cr
#' This function generates a Beeswarm plot of the sample values in a given latent factor. \cr
#' Similar functions are \code{\link{plotFactorScatter}} for doing scatter plots and \code{\link{plotFactorHist}} for doing histogram plots
#' @return Returns a \code{ggplot2} object
#' @import ggplot2
#' @import ggbeeswarm
#' @import grDevices
#' @export
plotFactorBeeswarm <- function(object, factors = "all", color_by = NULL, name_color = "", showMissing = FALSE) {
  
  # Sanity checks
  if (!is(object, "BioFAModel")) stop("'object' has to be an instance of BioFAModel")

  # Collect relevant data
  N <- sum(object@Dimensions[["N"]])
  Z <- getFactors(object, factors=factors, include_intercept=FALSE, as.data.frame=T)
  Z$factor <- as.factor(Z$factor)
  
  # Set color
  colorLegend <- T
  if (!is.null(color_by)) {
    # It is the name of a covariate or a feature in the TrainData
    if (length(color_by) == 1 & is.character(color_by)) {
      if(name_color=="") name_color <- color_by
      TrainData <- getTrainData(object)
      featureNames <- lapply(TrainData(object), rownames)
      if(color_by %in% Reduce(union,featureNames)) {
        viewidx <- which(sapply(featureNames, function(vnm) color_by %in% vnm))
        color_by <- TrainData[[viewidx]][color_by,]
      } else if(class(object@InputData) == "MultiAssayExperiment") {
        color_by <- getCovariates(object, color_by)
      } else stop("'color_by' was specified but it was not recognised, please read the documentation")
    # It is a vector of length N
    } else if (length(color_by) > 1) {
      stopifnot(length(color_by) == N)
      # color_by <- as.factor(color_by)
    } else {
      stop("'color_by' was specified but it was not recognised, please read the documentation")
    }
  } else {
    color_by <- rep(TRUE, N)
    colorLegend <- F
  }
  names(color_by) <- unlist(sampleNames(object))
  if(length(unique(color_by)) < 5) color_by <- as.factor(color_by)

  Z$color_by <- color_by[Z$sample]
  
  # Remove samples with missing values
  if (!showMissing) {
    Z <- Z[complete.cases(Z),]
    # Z <- Z[!is.na(Z$value),]
  }
  
  # Generate plot
  p <- ggplot(Z, aes_string(x=0, y="value")) + 
    ggbeeswarm::geom_quasirandom(aes(color=color_by)) +
    ylab("Factor value") + xlab("") +
    scale_x_continuous(breaks=NULL) +
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
  if (colorLegend) { p <- p + labs(color=name_color) } else { p <- p + guides(color = FALSE) }
  
  return(p)
}

#' @title Scatterplot of two latent factors
#' @name plotFactorScatter
#' @description Scatterplot of the values of two latent factors.
#' @param object a trained \code{\link{BioFAModel}} object.
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
#' Similar functions are \code{\link{plotFactorScatters}} for doing multiple scatter plots and 
#' \code{\link{plotFactorBeeswarm}} for doing Beeswarm plots for single factors.
#' @return Returns a \code{ggplot2} object
#' @import ggplot2
#' @export
plotFactorScatter <- function (object, factors, color_by = NULL, shape_by = NULL, name_color="",
                         name_shape="", showMissing = TRUE) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  stopifnot(length(factors)==2)
  stopifnot(all(factors %in% factorNames(object)))
  
  # Collect relevant data  
  N <- sum(object@Dimensions[["N"]])
  Z <- getFactors(object, factors = factors, include_intercept = FALSE, as.data.frame = T)
  Z$factor <- as.factor(Z$factor)
  
  # Set color
  colorLegend <- T
  if (!is.null(color_by)) {
    # It is the name of a covariate or a feature in the TrainData
    if (length(color_by) == 1 & is.character(color_by)) {
      if(name_color=="") name_color <- color_by
      TrainData <- getTrainData(object)
      featureNames <- lapply(TrainData(object), rownames)
      if(color_by %in% Reduce(union,featureNames)) {
        viewidx <- which(sapply(featureNames, function(vnm) color_by %in% vnm))
        color_by <- TrainData[[viewidx]][color_by,]
      } else if(class(object@InputData) == "MultiAssayExperiment"){
        color_by <- getCovariates(object, color_by)
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
    colorLegend <- F
  }

  # Set shape
  shapeLegend <- T
  if (!is.null(shape_by)) {
    # It is the name of a covariate 
    if (length(shape_by) == 1 & is.character(shape_by)) {
      if(name_shape=="") name_shape <- shape_by
      TrainData <- getTrainData(object)
      featureNames <- lapply(TrainData(object), rownames)
      if(shape_by %in% Reduce(union,featureNames)) {
        viewidx <- which(sapply(featureNames, function(vnm) shape_by %in% vnm))
        shape_by <- TrainData[[viewidx]][shape_by,]
      } else if(class(object@InputData) == "MultiAssayExperiment"){
        shape_by <- getCovariates(object, shape_by)
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
    shapeLegend <- F
  }
  
  # Create data frame to plot
  # NOTE: factor are concatenated across groupes
  # tmp <- lapply(Z, function(z) data.frame(x = z[, factors[1]], y = z[, factors[2]], shape_by = shape_by, color_by = color_by))
  # df <- do.call(rbind, tmp)
  df <- data.frame(x = Z[Z$factor == factors[[1]],'value'], y = Z[Z$factor == factors[[2]],'value'], shape_by = shape_by, color_by = color_by)
  
  # remove values missing color or shape annotation
  if (!showMissing) df <- df[!is.na(df$shape_by) & !is.na(df$color_by),]

   #turn into factors
   df$shape_by[is.na(df$shape_by)] <- "NA"
   df$shape_by <- as.factor(df$shape_by)
   if(length(unique(df$color_by)) < 5) df$color_by <- as.factor(df$color_by)
 
  
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
  
  
#' @title Pairwise scatterplots of multiple latent factors
#' @name plotFactorScatters
#' @description Scatterplots of the sample values for pair-wise combinations of multiple latent factors.
#' @param object a \code{\link{BioFAModel}} object.
#' @param factors character vector with the factor name(s), or numeric vector with the index of the factor(s) to use. 
#' Default is 'all'
#' @param color_by specifies groups or values used to color the samples. 
#' This can be either: 
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
#' @param showMissing logical indicating whether to include samples for which \code{shape_by} or \code{color_by} is missing
#' @details One of the first steps for the annotation of factors is to visualise and group/color them using known covariates such as phenotypic or clinical data.
#' This method generates multiple scatterplots for pairwise combinations of several latent factors.
#' Similar functions are \code{\link{plotFactorScatter}} for doing single scatter plots and 
#' \code{\link{plotFactorBeeswarm}} for doing Beeswarm plots for single factors.
#' @return \code{ggplot2} object
#' @import ggplot2
#' @export
plotFactorScatters <- function(object, factors = "all", showMissing=TRUE, 
                         color_by=NULL, name_color="",  
                         shape_by=NULL, name_shape="") {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")

  # Get factors
  if (paste0(factors, collapse="") == "all") { 
    factors <- factorNames(object) 
    # if(is.null(factors)) factors <- 1:ncol(Z) # old object are not compatible with factro names
  } else {
    stopifnot(all(factors %in% factorNames(object)))  
  }

  # Collect relevant data
  N <- object@Dimensions[["N"]]
  Z <- getFactors(object, factors=factors, include_intercept=FALSE, as.data.frame=T)
  Z$factor <- as.factor(Z$factor)

  # Remove constant factors 
  tmp <- apply(Z, 2, var,na.rm=T)
  if (any(tmp==0)) {
    # message(paste0("Removing constant factors: ", paste(which(tmp==0), collapse="")))
    Z <- Z[,!tmp==0]
    factors <- factors[!tmp==0]
  }
  
  # Set color
  colorLegend <- T
  if (!is.null(color_by)) {
    # It is the name of a covariate or a feature in the TrainData
    if (length(color_by) == 1 & is.character(color_by)) {
      if(name_color=="") name_color <- color_by
      TrainData <- getTrainData(object)
      featureNames <- lapply(TrainData(object), rownames)
      if(color_by %in% Reduce(union,featureNames)) {
        viewidx <- which(sapply(featureNames, function(vnm) color_by %in% vnm))
        color_by <- TrainData[[viewidx]][color_by,]
      } else if(class(object@InputData) == "MultiAssayExperiment"){
        color_by <- getCovariates(object, color_by)
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
    colorLegend <- F
  }

  # Set shape
  shapeLegend <- T
  if (!is.null(shape_by)) {
    # It is the name of a covariate 
    if (length(shape_by) == 1 & is.character(shape_by)) {
      if(name_shape=="") name_shape <- shape_by
      TrainData <- getTrainData(object)
      featureNames <- lapply(TrainData(object), rownames)
      if (shape_by %in% Reduce(union,featureNames)) {
        viewidx <- which(sapply(featureNames, function(vnm) shape_by %in% vnm))
        shape_by <- TrainData[[viewidx]][shape_by,]
      } else if(class(object@InputData) == "MultiAssayExperiment"){
        shape_by <- getCovariates(object, shape_by)
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
    shapeLegend <- F
  }

  # Remove missing values
  if(!showMissing) {
    Z <- Z[!is.na(color_by),]
    color_by <- color_by[!is.na(color_by)]
    shape_by <- shape_by[!is.na(shape_by)]
  }

  # Crete data.frame
  df <- as.data.frame(Z); colnames(df) <- paste0("LF_",colnames(df))
  df <- cbind(df, color_by=color_by, shape_by=shape_by)

  #turn into factors
  df$shape_by[is.na(df$shape_by)] <- "NA"
  df$shape_by <- as.factor(df$shape_by)
  if(length(unique(df$color_by)) < 5) df$color_by <- as.factor(df$color_by)
  
  # Define title and legend of the plot
  main <- "" 
  p <- ggplot(df, aes_string(x=colnames(df)[1], y=colnames(df)[2], color="color_by", shape="shape_by")) + geom_point()
  if (colorLegend | shapeLegend) { 
    p <- p +
      theme(
        legend.title=element_text(size=15, hjust=0.5, color="black"),
        legend.position = "right", 
        legend.direction = "vertical",
        legend.key = element_blank()
      )
    
    # If color_by is numeric, define the default gradient
    if (is.numeric(df$color_by)) { p <- p + scale_color_gradientn(colors=terrain.colors(10)) }
    
    if (colorLegend) { p <- p + labs(color = name_color) } else { p <- p + guides(color = FALSE) }
    if (shapeLegend) { p <- p + labs(shape = name_shape) }  else { p <- p + guides(shape = FALSE) }
    # Extract the legend
    legend <- GGally::grab_legend(p)
  } else {
    legend <- NULL
  }
  
  # Generate plot
  p <- GGally::ggpairs(df, columns = colnames(df[,!colnames(df) %in% c("color_by","shape_by")]), 
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
#' @name plotFactorCor
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
plotFactorCor <- function(object, method = "pearson", ...) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  
  # Fetch factors
  Z <- getFactors(object)
  
  # Remove intercept
  if(object@ModelOptions$LearnIntercept) Z <- lapply(Z, function(z) z[,-1])
  
  # Compute and plot correlation
  r <- abs(cor(x=do.call(rbind, Z), y=do.call(rbind, Z), method=method, use = "complete.obs"))
  p <- corrplot(r, tl.col="black", ...)
  
  return(r)
}


