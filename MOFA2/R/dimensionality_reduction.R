
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
# #' @importFrom Rtsne Rtsne
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("exdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Run t-SNE
#' \dontrun{ model <- run_tsne(model) }
#' 
#' # Change hyperparameters passed to Rtsne
#' \dontrun{ model <- run_tsne(model, perplexity = 15) }
#' 
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
  tsne_embedding <- Rtsne::Rtsne(Z, check_duplicates = FALSE, pca = FALSE, ...)

  # Add sample names and enumerate latent dimensions (e.g. TSNE1 and TSNE2)
  object@dim_red$TSNE <- data.frame(rownames(Z), tsne_embedding$Y)
  colnames(object@dim_red$TSNE) <- c("sample", paste0("TSNE", 1:ncol(tsne_embedding$Y)))
  
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
# #' @importFrom uwot umap
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("exdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Run UMAP
#' \dontrun{ model <- run_umap(model) }
#' 
#' # Change hyperparameters passed to umap
#' \dontrun{ model <- run_umap(model, min_dist = 0.01, n_neighbors = 10) }
#' 
run_umap <- function(object, factors = "all", groups = "all", ...) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get factor values
  Z <- get_factors(object, factors = factors, groups = groups)
  
  # Concatenate groups
  Z <- do.call(rbind, Z)
  
  # Replace missing values by zero
  Z[is.na(Z)] <- 0
  
  # Run UMAP
  umap_embedding <- uwot::umap(Z, ...)

  # Add sample names and enumerate latent dimensions (e.g. UMAP1 and UMAP2)
  object@dim_red$UMAP <- data.frame(rownames(Z), umap_embedding)
  colnames(object@dim_red$UMAP) <- c("sample", paste0("UMAP", 1:ncol(umap_embedding)))
  
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
#' @param alpha_missing numeric indicating dot transparency of missing data.
#' @param legend logical indicating whether to add legend.
#' @param return_data logical indicating whether to return the long data frame to plot instead of plotting
#' @param ... extra arguments passed to \code{\link{run_umap}} or \code{\link{run_tsne}}.
#' @details This function plots dimensionality reduction projections that are stored in the \code{dim_red} slot.
#' Typically this contains UMAP or t-SNE projections computed using \code{\link{run_tsne}} or \code{\link{run_umap}}, respectively.
#' @return Returns a \code{ggplot2} object or a long data.frame (if return_data is TRUE)
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom stats complete.cases
#' @importFrom tidyr spread gather
#' @importFrom magrittr %>% set_colnames
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("exdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Run UMAP
#' model <- run_umap(model)
#' 
#' # Plot UMAP
#' plot_dimred(model, method = "UMAP")
#' 
#' # Plot UMAP, colour by Factor 1 values
#' plot_dimred(model, method = "UMAP", color_by = "Factor1")
#' 
#' # Plot UMAP, colour by the values of a specific feature
#' plot_dimred(model, method = "UMAP", color_by = "feature_0_view_0")
#' 
plot_dimred <- function(object, method = c("UMAP", "TSNE"), groups = "all", show_missing = TRUE,
                        color_by = NULL, shape_by = NULL, color_name = NULL, shape_name = NULL,
                        dot_size = 1.5, alpha_missing = 1, legend = TRUE, return_data = FALSE, ...) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")

  # If UMAP or TSNE is requested but were not computed, compute the requested embedding
  if ((method %in% c("UMAP", "TSNE")) && (!.hasSlot(object, "dim_red") || !(method %in% names(object@dim_red)))) {
    message(paste0(method, " embedding was not computed. Running run_", tolower(method), "()..."))
    if (method == "UMAP") {
      object <- run_umap(object, ...)
    } else {  # method == "TSNE"
      object <- run_tsne(object, ...)
    }
  }

  # make sure the slot for the requested method exists
  method <- match.arg(method, names(object@dim_red))  
  
  # Remember color_name and shape_name if not provided
  if (!is.null(color_by) && (length(color_by) == 1) && is.null(color_name))
    color_name <- color_by
  if (!is.null(shape_by) && (length(shape_by) == 1) && is.null(shape_name))
    shape_name <- shape_by
  
  # Fetch latent manifold
  Z <- object@dim_red[[method]]
  latent_dimensions_names <- colnames(Z)[-1]
  Z <- gather(Z, -sample, key="latent_dimension", value="value")
  
  # Subset groups
  groups <- .check_and_get_groups(object, groups)
  Z <- Z[Z$sample%in%unlist(samples_names(object)[groups]),]
  
  # Set color and shape
  color_by <- .set_colorby(object, color_by)
  shape_by <- .set_shapeby(object, shape_by)
  
  # Merge factor values with color and shape information
  df <- merge(Z, color_by, by="sample")
  df <- merge(df, shape_by, by="sample")
  df$shape_by <- as.character(df$shape_by)
  
  # Remove missing values
  if(!show_missing) df <- filter(df, !is.na(color_by) & !is.na(shape_by))
  df$observed <- as.factor(!is.na(df$color_by))
  
  # spread over latent dimensions
  df <- spread(df, key="latent_dimension", value="value")
  df <- set_colnames(df, c(colnames(df)[seq_len(4)], "x", "y"))
  
  # Return data if requested instead of plotting
  if (return_data) return(df)
  
  # Generate plot
  p <- ggplot(df, aes_string(x = "x", y = "y")) + 
    geom_point(aes_string(color = "color_by", shape = "shape_by", alpha = "observed"), size = dot_size) +
    labs(x = latent_dimensions_names[1], y = latent_dimensions_names[2]) +
    theme_classic() +
    theme(
      axis.text = element_text(size = rel(0.9), color = "black"), 
      axis.title = element_text(size = rel(1.2), color = "black"), 
      axis.line = element_line(color = "black", size = 0.5), 
      axis.ticks = element_line(color = "black", size = 0.5)
    )
  
  # Add legend for color
  if (length(unique(df$observed))>1) { 
    p <- p + scale_alpha_manual(values = c("TRUE"=1.0, "FALSE"=alpha_missing))
  } else { 
    p <- p + scale_alpha_manual(values = 1.0)
  }
  p <- p + guides(alpha=FALSE)
    
  # If color is numeric, define the default gradient
  if (is.numeric(df$color))
    # p <- p + scale_color_gradientn(colors=colorRampPalette(rev(brewer.pal(n = 5, name = "RdYlBu")))(10)) 
    p <- p + scale_color_gradientn(colours = c('lightgrey', 'blue'))
  
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
