
##################################################################
## Functions to do dimensionality reduction on the MOFA factors ##
##################################################################

#' @title Run t-SNE on the MOFA factors
#' @name run_tsne
#' @param object a trained \code{\link{MOFA}} object.
#' @param factors character vector with the factor names, or numeric vector with the indices of the factors to use, or "all" to plot all factors.
#' @param groups character vector with the groups names, or numeric vector with the indices of the groups of samples to use, or "all" to use samples from all groups.
#' @param ... arguments passed to \code{\link{Rtsne}}
#' @details use set.seed before the function call to get reproducible results.
#' @return Returns a \code{\link{MOFA}} object with the dim_red slot filled with the t-SNE output
#' @importFrom Rtsne Rtsne
#' @export
run_tsne <- function(object, factors = "all", groups = "all", ...) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get factor values
  Z <- get_factors(object, factors=factors, groups=groups)
  
  # Concatenate groups
  Z <- do.call(rbind, Z)
  
  # Replace missing values by zero
  Z[is.na(Z)] <- 0
  
  # Run t-SNE
  object@dim_red$TSNE <- Rtsne(Z, check_duplicates=FALSE, pca=FALSE, ...)
  
  return(object)
  
}



#' @title Run UMAP on the MOFA factors
#' @name run_umap
#' @param object a trained \code{\link{MOFA}} object.
#' @param factors character vector with the factor names, or numeric vector with the indices of the factors to use, or "all" to plot all factors.
#' @param groups character vector with the groups names, or numeric vector with the indices of the groups of samples to use, or "all" to use samples from all groups.
#' @param ... arguments passed to \code{\link{umap}}
#' @details use set.seed before the function call to get reproducible results.
#' @return Returns a \code{\link{MOFA}} object with the dim_red slot filled with the UMAP output
#' @importFrom uwot umap
#' @export
run_umap <- function(object, factors = "all", groups = "all", ...) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get factor values
  Z <- get_factors(object, factors=factors, groups=groups)
  
  # Concatenate groups
  Z <- do.call(rbind, Z)
  
  # Replace missing values by zero
  Z[is.na(Z)] <- 0
  
  # Run UMAP
  object@dim_red$UMAP <- umap(Z, ...)
  
  return(object)
  
}



#' @title Plot dimensionality reduction based on MOFA factors
#' @name plot_dimred
#' @param object a trained \code{\link{MOFA}} object.
#' @param method string indicating which method has been used for non-linear dimensionality reduction (either 'umap' or 'tsne')
#' @param groups character vector with the groups names, or numeric vector with the indices of the groups of samples to use, or "all" to use samples from all groups.
#' @param show_missing logical indicating whether to include samples for which \code{shape_by} or \code{color_by} is missing
#' @param color_by specifies groups or values used to color the samples. This can be either:
#' (1) a character giving the name of a feature present in the training data.
#' (2) a character giving the same of a column present in the sample metadata.
#' (3) a vector of the same length as the number of samples specifying discrete groups or continuous numeric values.
#' @param shape_by specifies groups or values used to shape the samples. This can be either:
#' (1) a character giving the name of a feature present in the training data, 
#' (2) a character giving the same of a column present in the sample metadata.
#' (3) a vector of the same length as the number of samples specifying discrete groups.
#' @param color_name name for color legend.
#' @param shape_name name for shape legend.
#' @param dot_size numeric indicating dot size.
#' @param alpha numeric indicating dot transparency.
#' @param legend logical indicating whether to add legend.
#' @param return_data logical indicating whether to return the long data frame to plot instead of plotting
#' @details TO-FINISH...
#' @return Returns a \code{ggplot2} object or a long data.frame (if return_data is TRUE)
#' @import ggplot2 dplyr
#' @importFrom stats complete.cases
#' @importFrom tidyr spread
#' @importFrom magrittr %>% set_colnames
#' @export
plot_dimred <- function(object, method = c("umap","tsne"), groups = "all", show_missing = TRUE,
                         color_by = NULL, shape_by = NULL, color_name = NULL, shape_name = NULL,
                         dot_size = 1.5, alpha = 1, legend = TRUE, return_data = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  method = match.arg(method)
  if (!method %in% names(object@dim_red)) stop(sprintf("%s not present in object@dim_red",method))
  
  # Remember color_name and shape_name if not provided
  if (!is.null(color_by) && (length(color_by) == 1) && is.null(color_name))
    color_name <- color_by
  if (!is.null(shape_by) && (length(shape_by) == 1) && is.null(shape_name))
    shape_name <- shape_by
  
  # Fetch latent manifold
  Z <- object@dim_red[[method]]
  
  # Set color and shape
  color_by <- .set_colorby(object, color_by)
  shape_by <- .set_shapeby(object, shape_by )
  
  # Merge factor values with color and shape information
  df <- merge(Z, color_by, by="sample")
  df <- merge(df, shape_by, by="sample")
  df$shape_by <- as.character(df$shape_by)
  
  # Remove missing values
  if(!show_missing) df <- filter(df, !is.na(color_by) & !is.na(shape_by))
  
  # spread over factors
  df <- spread(df, key="factor", value="value")
  df <- df[,c(colnames(df)[seq_len(4)], factors)]
  df <- set_colnames(df, c(colnames(df)[seq_len(4)], "x", "y"))
  
  # Return data if requested instead of plotting
  if (return_data) return(df)
  
  # Generate plot
  p <- ggplot(df, aes_string(x="x", y="y")) + 
    geom_point(aes_string(color = "color_by", shape = "shape_by"), size=dot_size, alpha=alpha) +
    # ggrastr::geom_point_rast(aes_string(color = "color_by", shape = "shape_by")) +
    labs(x=factors[1], y=factors[2]) +
    theme_classic() +
    theme(
      axis.text = element_text(size = rel(0.9), color = "black"), 
      axis.title = element_text(size = rel(1.2), color = "black"), 
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
    p <- p + guides(color=FALSE) + scale_color_manual(values="black")
  }
  
  # Add legend for shape
  if (length(unique(df$shape))>1) { 
    p <- p + labs(shape=shape_name)
  } else { 
    p <- p + guides(shape=FALSE) 
  }
  
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
