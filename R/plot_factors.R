
###########################################
## Functions to visualise latent factors ##
###########################################

#' @title Beeswarm plot of factor values
#' @name plot_factor
#' @description Beeswarm plot of the latent factor values.
#' @param object a trained \code{\link{MOFA}} object.
#' @param factors character vector with the factor names, or numeric vector with the indices of the factors to use, or "all" to plot all factors.
#' @param groups character vector with the groups names, or numeric vector with the indices of the groups of samples to use, or "all" to use samples from all groups.
#' @param scale logical indicating whether to scale factor values.
#' @param group_by specifies grouping of samples:
#' \itemize{
#' \item (default) the string "group": in this case, the plot will color samples with respect to their predefined groups.
#' \item a character giving the name of a feature that is present in the input data 
#' \item a character giving the name of a column in the sample metadata slot
#' \item a vector of the same length as the number of samples specifying the value for each sample. 
#'}
#' @param color_by specifies color of samples. This can be either: 
#' \itemize{
#' \item (default) the string "group": in this case, the plot will color the dots with respect to their predefined groups.
#' \item a character giving the name of a feature that is present in the input data 
#' \item a character giving the name of a column in the sample metadata slot
#' \item a vector of the same length as the number of samples specifying the value for each sample. 
#' }
#' @param shape_by specifies shape of samples. This can be either: 
#' \itemize{
#' \item (default) the string "group": in this case, the plot will shape the dots with respect to their predefined groups.
#' \item a character giving the name of a feature that is present in the input data 
#' \item a character giving the name of a column in the sample metadata slot
#' \item a vector of the same length as the number of samples specifying the value for each sample. 
#' }
#' @param add_dots logical indicating whether to add dots.
#' @param add_violin logical indicating whether to add violin plots
#' @param add_boxplot logical indicating whether to add box plots
#' @param dodge logical indicating whether to dodge the dots (default is FALSE).
#' @param show_missing logical indicating whether to remove samples for which \code{shape_by} or \code{color_by} is missing.
#' @param dot_size numeric indicating dot size.
#' @param dot_alpha numeric indicating dot transparency.
#' @param violin_alpha numeric indicating violin plot transparency.
#' @param color_violin logical indicating whether to color violin plots.
#' @param boxplot_alpha numeric indicating boxplot transparency.
#' @param color_boxplot logical indicating whether to color box plots.
#' @param color_name name for color legend (usually only used if color_by is not a character itself).
#' @param shape_name name for shape legend (usually only used if shape_by is not a character itself).
#' @param stroke numeric indicating the stroke size (the black border around the dots).
#' @param legend logical indicating whether to add a legend to the plot (default is TRUE).
#' @param rasterize logical indicating whether to rasterize the plot (default is FALSE).
#' @details One of the main steps for the annotation of factors is to visualise and color them using known covariates or phenotypic data. \cr
#' This function generates a Beeswarm plot of the sample values in a given latent factor. \cr
#' Similar functions are \code{\link{plot_factors}} for doing scatter plots.
#' @return Returns a \code{ggplot2} 
#' @import ggplot2 grDevices
#' @importFrom stats complete.cases
#' @importFrom forcats fct_explicit_na
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr summarise group_by
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Plot Factors 1 and 2 and colour by "group"
#' plot_factor(model, factors = c(1,2), color_by="group")
#' 
#' # Plot Factor 3 and colour by the value of a specific feature
#' plot_factor(model, factors = 3, color_by="feature_981_view_1")
#' 
#' # Add violin plots
#' plot_factor(model, factors = c(1,2), color_by="group", add_violin = TRUE)
#' 
#' # Scale factor values from -1 to 1
#' plot_factor(model, factors = c(1,2), scale = TRUE)
#' 
plot_factor <- function(object, factors = 1, groups = "all",
                        group_by = "group", color_by = "group", shape_by = NULL, 
                        add_dots = TRUE, dot_size = 2, dot_alpha = 1,
                        add_violin = FALSE, violin_alpha = 0.5, color_violin = TRUE,
                        add_boxplot = FALSE, boxplot_alpha = 0.5, color_boxplot = TRUE,
                        show_missing = TRUE, scale = FALSE, dodge = FALSE,
                        color_name = "", shape_name = "", stroke = NULL,
                        legend = TRUE, rasterize = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Define factors
  factors <- .check_and_get_factors(object, factors)

  # Get factor values
  Z <- get_factors(object, factors=factors, groups = groups, as.data.frame=TRUE)
  Z$factor <- factor(Z$factor, levels=factors)
  
  # Set group/color/shape
  if (length(color_by)==1 & is.character(color_by)) color_name <- color_by
  if (length(shape_by)==1 & is.character(shape_by)) shape_name <- shape_by
  object_group_by <- .set_groupby(object, group_by)
  object_color_by <- .set_colorby(object, color_by)
  object_shape_by <- .set_shapeby(object, shape_by)
  
  # Remove samples with missing values
  Z <- Z[complete.cases(Z),]
  
  # Merge factor values with group/color/shape information
  df <- merge(Z,  object_group_by, by="sample")
  df <- merge(df, object_color_by, by="sample")
  df <- merge(df, object_shape_by, by="sample")
  
  # QC
  if (length(unique(df$color_by))>10 & is.numeric(df$color_by)) {
    add_violin <- FALSE 
    add_boxplot <- FALSE
    dodge <- FALSE
  }
  
  if (length(unique(df$shape_by))>5) {
    warning("Maximum number of shapes is 5")
    df$shape_by <- "1"
  }
  
  # if (all(unique(df$color_by)==unique(df$group_by))) dodge <- TRUE
  
  # Remove samples with no sample metadata
  if (!show_missing) df <- filter(df, !is.na(color_by) & !is.na(shape_by))
  if (is.factor(df$color_by))
    df$color_by <- forcats::fct_explicit_na(df$color_by)
  if (is.factor(df$shape_by))
    df$shape_by <- forcats::fct_explicit_na(df$shape_by)
  
  # Scale values
  if (scale) df$value <- df$value/max(abs(df$value))
  
  # Generate plot
  p <- ggplot(df, aes_string(x="group_by", y="value", fill="color_by", shape="shape_by")) +
    theme_classic()
  
  # Defien facets as factors or groups
  if (length(factors) == 1) {
    p <- p + facet_wrap(~group_by, nrow=1, scales="free_x") +
      labs(x=group_by, y=as.character(factors))
    if (length(unique(df$group_by))==1) p <- p + theme(strip.text = element_blank()) # remove facet title
  } else {
    p <- p + facet_wrap(~factor, nrow=1, scales="free_x") +
      labs(x=group_by, y="Factor value")
    if (length(unique(df$factor))==1) p <- p + theme(strip.text = element_blank()) # remove facet title
  }

  # Add dots
  if (add_dots) {
    
    # Set stroke
    if (is.null(stroke)) stroke <- .select_stroke(N=length(unique(df$sample)))
    
    if (rasterize) {
      warning("geom_jitter is not available with rasterize==TRUE. We use instead ggrastr::geom_quasirandom_rast()")
      if (dodge) {
        p <- p + ggrastr::geom_quasirandom_rast(size = dot_size, position = "dodge", stroke = stroke,  alpha = dot_alpha, dodge.width = 1)
      } else {
        p <- p + ggrastr::geom_quasirandom_rast(size = dot_size, stroke = stroke,  alpha = dot_alpha)
      }
    } else {
      if (dodge) {
        p <- p + geom_jitter(colour = "black", size = dot_size, stroke = stroke, alpha = dot_alpha, 
                  position = position_jitterdodge(dodge.width=1, jitter.width=0.2))
      } else {
        p <- p + geom_jitter(colour = "black", size = dot_size, stroke = stroke, alpha = dot_alpha)
      }
    }
  }
  
  # Add violin plot
  if (add_violin) {
    if (color_violin && dodge) {
      tmp <- summarise(group_by(df, factor, color_by), n=n())
      if (min(tmp$n)==1) {
        warning("Warning: some 'color_by' groups have only one observation, violin plots cannot be added. Adding boxplots instead...")
        add_boxplot <- TRUE
        # p <- p + geom_violin(color="black", fill="grey", alpha=violin_alpha, trim=TRUE, scale="width", show.legend = FALSE)
        # p <- p + geom_violin(color="black", alpha=violin_alpha, trim=TRUE, scale="width", show.legend = FALSE)
      } else {
        p <- p + geom_violin(alpha=violin_alpha, trim=TRUE, scale="width", position=position_dodge(width=1))
      }
      # p <- p + geom_violin(color="black", alpha=violin_alpha, trim=TRUE, scale="width", position=position_dodge(width=1), show.legend = FALSE)
    } else {
      p <- p + geom_violin(color="black", fill="grey", alpha=violin_alpha, trim=TRUE, scale="width", show.legend = FALSE)
    }
  }
  
  # Add boxplot plot
  if (add_boxplot) {
    if (color_boxplot && dodge) {
      tmp <- summarise(group_by(df, factor, color_by), n=n())
      # if (min(tmp$n)==1) {
      #   warning("Warning: some 'color_by' groups have only one observation, boxplot plots cannot be coloured")
      #   p <- p + geom_boxplot(color="black", alpha=boxplot_alpha, show.legend = FALSE)
      # } else {
      #   p <- p + geom_boxplot(alpha=boxplot_alpha, position=position_dodge(width=1), show.legend = FALSE)
      # }
      p <- p + geom_boxplot(color="black", alpha=boxplot_alpha, position=position_dodge(width=1))
    } else {
      p <- p + geom_boxplot(color="black", fill="grey", alpha=boxplot_alpha, show.legend = FALSE)
    }
  }
  
  # Add theme
  p <- p +
    geom_hline(yintercept=0, linetype="dashed", size=0.2, alpha=0.5) +
    theme(
        panel.border = element_rect(color="black", size=0.1, fill=NA),
        strip.background = element_rect(colour = "black", size=0.25),
        panel.spacing = unit(0,"lines"),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text = element_text(size=rel(0.75), color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=rel(1.0), color="black"),
        axis.line = element_line(color="black", size=0.25),
        axis.ticks = element_line(color = "black")
    )
  
  if (length(unique(df$factor))>1) {
    # p <- p + scale_y_continuous(breaks=NULL)
  } else {
    # Remove strip labels for groups, they are laballed along X axis
    if (isFALSE(dodge)) {
      p <- p + theme(strip.text.x = element_blank())
    }
  }
  
  # If group_by has a single value, remove text
  if (length(unique(df$group_by))==1) {
    p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  # Add legend
  p <- .add_legend(p, df, legend, color_name, shape_name)
  
  return(p)
}


#' @title Scatterplots of two factor values
#' @name plot_factors
#' @description Scatterplot of the values of two latent factors.
#' @param object a trained \code{\link{MOFA}} object.
#' @param factors a vector of length two with the factors to plot. Factors can be specified either as a characters
#' @param groups character vector with the groups names, or numeric vector with the indices of the groups of samples to use, or "all" to use samples from all groups.
#' @param show_missing logical indicating whether to include samples for which \code{shape_by} or \code{color_by} is missing
#' @param scale logical indicating whether to scale factor values.
#' @param color_by specifies groups or values used to color the samples. This can be either:
#' (1) a character giving the name of a feature present in the training data.
#' (2) a character giving the name of a column present in the sample metadata.
#' (3) a vector of the name length as the number of samples specifying discrete groups or continuous numeric values.
#' @param shape_by specifies groups or values used to shape the samples. This can be either:
#' (1) a character giving the name of a feature present in the training data, 
#' (2) a character giving the name of a column present in the sample metadata.
#' (3) a vector of the same length as the number of samples specifying discrete groups.
#' @param color_name name for color legend.
#' @param shape_name name for shape legend.
#' @param dot_size numeric indicating dot size (default is 1.5).
#' @param stroke numeric indicating the stroke size (the black border around the dots, default is NULL, infered automatically).
#' @param alpha numeric indicating dot transparency (default is 1).
#' @param legend logical indicating whether to add legend.
#' @param return_data logical indicating whether to return the data frame to plot instead of plotting
#' @details One of the first steps for the annotation of factors is to visualise and group/color them using known covariates such as phenotypic or clinical data.
#' This method generates a single scatterplot for the combination of two latent factors.
#' TO-FINISH...
#' \code{\link{plot_factors}} for doing Beeswarm plots for factors.
#' @return Returns a \code{ggplot2} object
#' @import ggplot2 dplyr
#' @importFrom stats complete.cases
#' @importFrom tidyr spread
#' @importFrom magrittr %>% set_colnames
# #' @importFrom ggbeeswarm geom_quasirandom
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Scatterplot of factors 1 and 2
#' plot_factors(model, factors = c(1,2))
#' 
#' # Shape dots by a column in the metadata
#' plot_factors(model, factors = c(1,2), shape_by="group")
#' 
#' # Scale factor values from -1 to 1
#' plot_factors(model, factors = c(1,2), scale = TRUE)
#' 
plot_factors <- function(object, factors = c(1, 2), groups = "all",
                         show_missing = TRUE, scale = FALSE,
                         color_by = NULL, shape_by = NULL, color_name = NULL, shape_name = NULL,
                         dot_size = 1.5, alpha = 1, legend = TRUE, stroke = NULL, return_data = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # If plotting one or multiple factors, re-direct to other functions 
  if (length(unique(factors)) == 1) {
    .args <- as.list(match.call()[-1])
    .args <- .args[names(.args) != "factors"]
    return(do.call(plot_factor, c(.args, list(factors = unique(factors)))))
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
  
  # Define factors
  factors <- .check_and_get_factors(object, factors)
  
  # Get factors
  Z <- get_factors(object, factors=factors, groups = groups, as.data.frame=TRUE)
  
  # Set color and shape
  color_by <- .set_colorby(object, color_by)
  shape_by <- .set_shapeby(object, shape_by)
  
  # Remove samples with missing values
  Z <- Z[complete.cases(Z),]
  
  # Merge factor values with color and shape information
  df <- merge(Z, color_by, by="sample")
  df <- merge(df, shape_by, by="sample")
  df$shape_by <- as.character(df$shape_by)
  
  # Remove missing values
  if(isFALSE(show_missing)) df <- filter(df, !is.na(color_by) & !is.na(shape_by))
  
  # spread over factors
  df <- spread(df, key="factor", value="value")
  df <- df[,c(colnames(df)[seq_len(4)], factors)]
  df <- set_colnames(df, c(colnames(df)[seq_len(4)], "x", "y"))

  # Scale values from 0 to 1
  if (scale) {
    df$x <- df$x/max(abs(df$x))
    df$y <- df$y/max(abs(df$y))
  }
  
  # Return data if requested instead of plotting
  if (return_data) return(df)
  
  # Set stroke
  if (is.null(stroke)) {
    stroke <- .select_stroke(N=length(unique(df$sample)))
  }
  
  # Generate plot
  p <- ggplot(df, aes_string(x="x", y="y", fill="color_by", shape="shape_by")) + 
    geom_point(size=dot_size, alpha=alpha, stroke = stroke) +
    labs(x=factors[1], y=factors[2]) +
    theme_classic() +
    theme(
      axis.text = element_text(size = rel(0.8), color = "black"), 
      axis.title = element_text(size = rel(1.1), color = "black"), 
      axis.line = element_line(color = "black", size = 0.5), 
      axis.ticks = element_line(color = "black", size = 0.5)
    )
  
  p <- .add_legend(p, df, legend, color_name, shape_name)

  # Fix legend labels
  if (!is.null(color_name)) {
    p <- p + labs(fill = color_name)
  }

  if (!is.null(shape_name)) {
    p <- p + labs(shape = shape_name)
  }

  return(p)
}


  
# Plot multiple factors as pairwise scatterplots
#' @importFrom stats complete.cases
.plot_multiple_factors <- function(object, factors = "all", show_missing = TRUE, dot_size = 1,
                                   color_by = NULL, color_name = "", shape_by = NULL, shape_name = "") {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Define factors
  factors <- .check_and_get_factors(object, factors)
  
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
  if (length(unique(df$color))>1) { p <- p + labs(color=color_name) } else { p <- p + guides(color=FALSE) + scale_color_manual(values="black") }
  if (is.numeric(df$color)) p <- p + scale_color_gradientn(colors=colorRampPalette(rev(brewer.pal(n=5, name="RdYlBu")))(10)) 
  if (length(unique(df$shape))>1) { p <- p + labs(shape=shape_name) } else { p <- p + guides(shape = FALSE) }
  if (length(unique(df$color))>1 || length(unique(df$shape))>1) { legend <- GGally::grab_legend(p) } else { legend <- NULL }
  
  
  # Generate plot
  p <- GGally::ggpairs(df, 
    columns = factors,
    lower = list(continuous=GGally::wrap("points", size=dot_size)), 
    diag = list(continuous='densityDiag'), 
    upper = list(continuous=GGally::wrap("points", size=dot_size)), 
    mapping = aes_string(color="color_by", shape="shape_by"), 
    title = "", 
    legend = legend
    )
  p <- p + theme_bw() + theme(
    # axis.line = element_line(color="black", size=rel(1.0)),
    panel.grid.major = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank()
  )
  
  return(p)
}
  


#' @title Plot correlation matrix between latent factors
#' @name plot_factor_cor
#' @description Function to plot the correlation matrix between the latent factors.
#' @param object a trained \code{\link{MOFA}} object.
#' @param method a character indicating the type of correlation coefficient to be computed: pearson (default), kendall, or spearman.
#' @param ... arguments passed to \code{\link[corrplot]{corrplot}}
#' @details This method plots the correlation matrix between the latent factors. \cr 
#' The model encourages the factors to be uncorrelated, so this function usually yields a diagonal correlation matrix. \cr 
#' However, it is not a hard constraint such as in Principal Component Analysis and correlations between factors can occur, 
#' particularly with large number factors. \cr
#' Generally, correlated factors are redundant and should be avoided, as they make interpretation harder. Therefore, 
#' if you have too many correlated factors we suggest you try reducing the number of factors.
#' @return Returns a symmetric matrix with the correlation coefficient between every pair of factors.
# #' @importFrom corrplot corrplot
#' @importFrom corrplot corrplot
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Plot correlation between all factors
#' plot_factor_cor(model)
#' 
plot_factor_cor <- function(object, method = "pearson", ...) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Fetch factors
  Z <- get_factors(object)
  
  # Compute and plot correlation
  r <- abs(cor(x=do.call(rbind, Z), y=do.call(rbind, Z), method=method, use = "complete.obs"))
  corrplot(r, tl.col = "black", ...)
}

