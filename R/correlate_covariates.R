#' @title Plot correlation of factors with external covariates
#' @name correlate_factors_with_covariates
#' @description Function to correlate factor values with external covariates.
#' @param object a trained \code{\link{MOFA}} object.
#' @param covariates 
#' \itemize{
#'   \item{\strong{data.frame}:}{a data.frame where the samples are stored in the rows and the covariates are stored in the columns. 
#'   Use row names for sample names and column names for covariate names. Columns values must be numeric. }
#'   \item{\strong{character vector}:}{character vector with names of columns that are present in the sample metadata (\code{samples_metadata(model)}}
#' }
#' @param factors character vector with the factor name(s), or numeric vector with the index of the factor(s) to use. Default is 'all'.
#' @param groups character vector with the groups names, or numeric vector with the indices of the groups of samples to use, or "all" to use samples from all groups.
#' @param abs logical indicating whether to take the absolute value of the correlation coefficient (default is \code{TRUE}).
#' @param plot character indicating whether to plot Pearson correlation coefficiens (\code{plot="r"}) or log10 adjusted p-values (\code{plot="log_pval"}).
#' @param return_data logical indicating whether to return the correlation results instead of plotting
#' @param transpose logical indicating whether to transpose the plot
#' @param alpha p-value threshold
#' @param ... extra arguments passed to \code{\link[corrplot]{corrplot}} (if \code{plot=="r"}) or \code{\link[pheatmap]{pheatmap}} (if \code{plot=="log_pval"}).
#' @importFrom pheatmap pheatmap
#' @importFrom corrplot corrplot
#' @return A \code{\link[corrplot]{corrplot}} (if \code{plot=="r"}) or \code{\link[pheatmap]{pheatmap}} (if \code{plot=="log_pval"}) or the underlying data.frame if return_data is TRUE
#' @export
correlate_factors_with_covariates <- function(object, covariates, factors = "all", groups = "all", 
                                              abs = FALSE, plot = c("log_pval","r"), 
                                              alpha = 0.05, return_data = FALSE, transpose = FALSE, ...) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  groups <- .check_and_get_groups(object,groups)
  
  # Get covariates
  metadata <- samples_metadata(object)
  metadata <- metadata[metadata$group%in%groups,]
  if (is.character(covariates)) {
    stopifnot(all(covariates %in% colnames(metadata)))
    covariates <- metadata[,covariates,drop=FALSE]
  } else if (is.data.frame(covariates)) {
    samples <- metadata$sample
    if (is.null(rownames(covariates))) stop("The 'covariates' data.frame does not have samples names")
    stopifnot(all(rownames(covariates) %in% samples))
    covariates <- metadata[match(rownames(covariates), metadata$sample),]
  } else {
    stop("covariates argument not recognised. Please read the documentation: ?correlate_factors_with_covariates")
  }
  
  # convert character columns to factors
  cols <- which(sapply(covariates, is.character))
  if (length(cols>=1)) {
    covariates[cols] <- lapply(covariates[cols], as.factor)
  }
  
  # convert all columns to numeric
  cols <- which(!sapply(covariates,class)%in%c("numeric","integer"))
  if (length(cols>=1)) {
    cols.factor <- which(sapply(covariates,class)=="factor")
    covariates[cols] <- lapply(covariates[cols], as.numeric)
    warning("There are non-numeric values in the covariates data.frame, converting to numeric...")
    covariates[cols] <- lapply(covariates[cols], as.numeric)
  }
  stopifnot(all(sapply(covariates,class)%in%c("numeric","integer")))
  
  # Get factors
  factors <- .check_and_get_factors(object, factors)
  Z <- get_factors(object, factors = factors, groups = groups, as.data.frame=FALSE)
  Z <- do.call(rbind, Z)
  
  # correlation
  cor <- psych::corr.test(Z, covariates, method = "pearson", adjust = "BH")
  
  # plot  
  plot <- match.arg(plot)
  
  if (plot=="r") {
    stat <- cor$r
    if (abs) stat <- abs(stat)
    if (transpose) stat <- t(stat)
    if (return_data) return(stat)
    corrplot::corrplot(stat, tl.col = "black", title="Pearson correlation coefficient", ...)
    
  } else if (plot=="log_pval") {
    stat <- cor$p
    stat[stat>alpha] <- 1.0
    if (all(stat==1.0)) stop("All p-values are 1.0, cannot plot the histogram")
    stat <- -log10(stat)
    stat[is.infinite(stat)] <- 1000
    if (transpose) stat <- t(stat)
    if (return_data) return(stat)
    col <- colorRampPalette(c("lightgrey", "red"))(n=100)
    pheatmap::pheatmap(stat, main="log10 adjusted p-values", cluster_rows = FALSE, color=col, ...)
    
  } else {
    stop("'plot' argument not recognised. Please read the documentation: ?correlate_factors_with_covariates")
  }
  
}



#' @title Summarise factor values using external groups
#' @name summarise_factors
#' @description Function to summarise factor values using a discrete grouping of samples.
#' @param object a trained \code{\link{MOFA}} object.
#' @param df a data.frame with the columns "sample" and "level", where level is a factor with discrete group assigments for each sample.
#' @param factors character vector with the factor name(s), or numeric vector with the index of the factor(s) to use. Default is 'all'.
#' @param groups character vector with the groups names, or numeric vector with the indices of the groups of samples to use, or "all" to use samples from all groups.
#' @param abs logical indicating whether to take the absolute value of the factors (default is \code{FALSE}).
#' @param return_data logical indicating whether to return the fa instead of plotting
#' @import ggplot2
#' @importFrom dplyr group_by summarise mutate
#' @importFrom stats median
#' @importFrom magrittr %>%
#' @return A \code{\link{ggplot}} object or a \code{data.frame} if return_data is TRUE
#' @export
summarise_factors <- function(object, df, factors = "all", groups = "all", abs = FALSE, return_data = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(is.data.frame(df))
  stopifnot((c("sample","level")%in%colnames(df)))
  stopifnot(df$sample %in% unlist(samples_names(object)))
  stopifnot(length(df$level)>1)
  df$level <- as.factor(df$level)
  
  # Get factors
  factors <- .check_and_get_factors(object, factors)
  groups <- .check_and_get_groups(object, groups)
  factors_df <- get_factors(object, factors = factors, groups = groups, as.data.frame=TRUE) %>% 
    group_by(factor) %>% mutate(value=value/max(abs(value),na.rm=TRUE)) # Scale factor values
  
  # Merge data.frames
  to.plot <- merge(factors_df, df, by="sample") %>% 
    group_by(level,factor,group) %>%
    summarise(value=median(value,na.rm=TRUE))
  
  if (abs) {
    to.plot$value <- abs(to.plot$value)
  }
  
  
  # Plot
  if (length(unique(factors_df$group))>1) {
    to.plot$group <- factor(to.plot$group, levels=groups)
    p <- ggplot(to.plot, aes_string(x="group", y="level", fill="value")) +
      facet_wrap(~factor)
  } else {
    p <- ggplot(to.plot, aes_string(x="factor", y="level", fill="value"))
  }
  
  p <- p +
    geom_tile() +
    theme_classic() +
    labs(x="", y="", fill="Score") +
    theme(
      axis.text.x = element_text(color="black", angle=30, hjust=1),
      axis.text.y = element_text(color="black")
    )

  if (abs) {
    p <- p + scale_fill_gradient2(low = "white", high = "red")
  } else {
    # center the color scheme at 0
    p <- p + scale_fill_distiller(type = "div", limit = max(abs(to.plot$value),na.rm=TRUE)*c(-1,1))
  } 
  
  # Return data or plot
  if (return_data) {
    return(to.plot)
  } else {
    return(p)
  }
}


