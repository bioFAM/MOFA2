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
#' @param abs logical indicating whether to take the absolute value of the correlation coefficient (default is TRUE).
#' @param plot character indicating whether to plot Pearson correlation coefficiens (\code{plot="r"}) or log10 adjusted p-values (\code{plot="log_pval"}).
#' @param return_data logical indicating whether to return the correlation results instead of plotting
#' @param transpose logical indicating whether to transpose the plot
#' @param ... extra arguments passed to \code{\link[corrplot]{corrplot}} (if plot=="r") or \code{\link[pheatmap]{pheatmap}} (if plot=="log_pval").
#' @importFrom pheatmap pheatmap
#' @export
correlate_factors_with_covariates <- function(object, covariates, factors = "all", groups = "all", abs = TRUE, plot = c("log_pval","r"), return_data = FALSE, 
                                              transpose = FALSE, ...) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get covariates
  metadata <- samples_metadata(object)
  if (is.character(covariates)) {
    stopifnot(all(covariates %in% colnames(metadata)))
    covariates <- metadata[,covariates,drop=F]
  } else if (is.data.frame(covariates)) {
    samples <- metadata$sample
    if (is.null(rownames(covariates))) stop("The 'covariates' data.frame does not have samples names")
    stopifnot(all(rownames(covariates) %in% samples))
    metadata <- metadata[match(rownames(covariates), metadata$sample),]
  } else {
    stop("covariates argument not recognised. Please read the documentation: ?correlate_factors_with_covariates")
  }
    
  # convert character columns to factors
  cols <- which(sapply(covariates,class)=="character")
  if (length(cols>=1)) {
    covariates[cols] <- lapply(covariates[cols], as.factor)
  }
  
  # convert all columns to numeric
  cols <- which(sapply(covariates,class)!="numeric")
  if (length(cols>=1)) {
    cols.factor <- which(sapply(covariates,class)=="factor")
    covariates[cols] <- lapply(covariates[cols], as.numeric)
    warning("There are non-numeric values in the covariates data.frame, converting to numeric...")
    covariates[cols] <- lapply(covariates[cols], as.numeric)
  }
  stopifnot(all(sapply(covariates,class)=="numeric"))
  
  # Get factors
  factors <- .check_and_get_factors(object, factors)
  Z <- get_factors(object, factors = factors, groups = groups, as.data.frame=FALSE)
  Z <- do.call(rbind, Z)
  
  # correlate
  cor <- psych::corr.test(Z, covariates, method = "pearson", adjust = "BH")
  
  # plot  
  plot <- match.arg(plot)
  
  if (plot=="r") {
    stat <- cor$r
    if (isTRUE(transpose)) stat <- t(stat)
    if (isTRUE(return_data)) return(stat)
    corrplot::corrplot(stat, ...)
    
  } else if (plot=="log_pval") {
    stat <- -log10(cor$p)
    if (isTRUE(transpose)) stat <- t(stat)
    if (isTRUE(return_data)) return(stat)
    col <- colorRampPalette(c("lightgrey", "red"))(n=100)
    pheatmap::pheatmap(stat, color=col, ...)
    
  } else {
    stop("'plot' argument not recognised. Please read the documentation: ?correlate_factors_with_covariates")
  }
  
}

