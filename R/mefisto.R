##########################################################################
## Functions to use smooth covariates, as part of the MEFISTO framework ##
##########################################################################

#' @title Add covariates to a MOFA model
#' @name set_covariates
#' @description Function to add  continuous covariate to a \code{\link{MOFA}} object for smooth training (MEFISTO) 
#' @param object an untrained \code{\link{MOFA}}
#' @param covariates list of data_options (see \code{\link{get_default_data_options}} details). 
#' @return Returns an untrained \code{\link{MOFA}} with covariates filled in the corresponding slots
#' @details To activate the functional MEFISTO framework, specify smooth_options when preparing the training using \code{prepare_mofa} 
#' @export
#' @examples
#' #' # Simulate data
#' dd <- make_example_data(sample_cov = seq(0,1,length.out = 100), n_samples = 100, n_factors = 4)
#' 
#' # Create MOFA object
#' sm <- create_mofa(data = dd$data)
#' 
#' Add a covariate
#' sm <- set_covariates(sm, covariates = dd$sample_cov)
#' sm

set_covariates <- function(object, covariates) {

  # Sanity checks
  if (!is(object, "MOFA")) 
    stop("'object' has to be an instance of MOFA")
  if (object@status=="trained") 
    stop("The model is already trained! Covariates must be added before training")
  
  # get sample names
  samples_data <- lapply(object@data[[1]], colnames)
  # samples <- unlist(samples_data)
  samples_data_vec <- unlist(samples_names(object))
  
  # covariates passed as characters: extract from the metadata as dataframe
  if (is(covariates, "character")) {
    if (!all(covariates %in% colnames(samples_metadata(object)))) {
      stop("Columns specified in covariates do not exist in the MOFA object metadata slot.")
    }
    covariates <- samples_metadata(object)[,c("sample",covariates),drop=F]

    covariates <- gather(covariates, key = "covariate", value = "value", -sample)
    if(!is.numeric(covariates$value)){
      stop("Covariates need to be numeric")
    }
    # TO-DO: Check that they continuous
  
  # covariates passed in data.frame format
  }
  
  if (any(class(covariates) %in% c("data.frame", "tibble", "Data.Frame"))) { # TO-DO: USE is()
    if (!all(c("sample", "covariate", "value") %in% colnames(covariates)))
      stop("If covariates is provided as data.frame it needs to contain the columns: sample, covariate, value")
    if (!is.numeric(covariates$value)) {
      stop("Values in covariates need to be numeric")
    }
    samples <- covariates$sample
    # covariates <- covariates[!duplicated(covariates), ]
    covariates <- reshape2::acast(covariates, covariate ~ sample)
    
  # covariates passed in matrix format
  # TO-DO: CHECK THIS
  } else if (all(is.numeric(covariates)) || class(covariates) %in% c("dgTMatrix", "dgCMatrix")) {
    samples <- colnames(covariates)
    if (!is.null(samples)) {
      if(!(all(samples %in% samples_data_vec) && all(samples_data_vec %in% samples)))
        stop("Sample names of the data and the sample covariates do not match.")
      covariates <- covariates[ , samples_data_vec, drop = FALSE]
    } else {
      # warnings and checks if no matching sample names
      if(sum(object@dimensions[['N']]) != ncol(covariates))
        stop("Number of columns in sample covariates does not match the number of samples")
      if(!is.null(samples_data) && length(samples_data_vec) > 0) {
        warning("No sample names in covariates - we will use the sample names in data. Please ensure that the order matches.")
        colnames(covariates) <- samples
      } else {
        stop("No sample names found!")
      }
    }
    
  # covariates format not recognised
  } else {
    stop("covariates needs to be a character vector, a dataframe, a matrix or NULL.")
  }
    
  # Set covariate dimensionality
  object@dimensions[["C"]] <- nrow(covariates)
    
  # Set covariate names
  if (is.null(rownames(covariates))) {
    message("No covariates names provided - using generic: covariate1, covariate2, ...")
    rownames(covariates) <- paste0("covariate", seq_len(nrow(covariates)))
  }
  
  # split covariates by groups
  covariates <- lapply(samples_names(object), function(i)   covariates[, i, drop = FALSE])
  names(covariates) <- groups_names(object)
  
  # Sanity checks
  stopifnot(all(sapply(object@covariates, ncol) == object@dimensions[["N"]]))
  
  # add covariates to the MOFA object
  object@covariates <- covariates
  
  return(object)
}


#' @title Get sample covariates
#' @name get_covariates
#' @description Function to extract the covariates from a \code{\link{MOFA}} object.
#' @param object a \code{\link{MOFA}} object.
#' @param covariates character vector with the covariate name(s), or numeric vector with the covariate index(es). 
#' @param as.data.frame logical indicating whether to output the result as a long data frame, default is \code{FALSE}.
#' @param warped logical indicating whether to extract the aligned covariates
#' @return a matrix with dimensions (samples,covariates). If \code{as.data.frame} is \code{TRUE}, a long-formatted data frame with columns (sample,factor,value)
#' @export

get_covariates <- function(object, covariates = "all", as.data.frame = FALSE, warped = FALSE) {
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get and check covariate names
  covariates <- .check_and_get_covariates(object, covariates)
  
  # Get covariates
  if(warped){
    sample_cov <- lapply(object@covariates_warped, function(cmat) cmat[covariates,,drop=FALSE])
  } else {
    sample_cov <- lapply(object@covariates, function(cmat) cmat[covariates,,drop=FALSE])
    
  }
  
  if (as.data.frame) {
    if(!is.null(rownames(sample_cov[[1]]))){
      nms <- rownames(sample_cov[[1]]) 
    } else {
      nms <- paste0("covariate_", seq_along(covariates))
    }
    sample_cov <- Reduce(cbind, sample_cov) # remove group info
    sample_cov <- melt(sample_cov, varnames = c("covariate", "sample"))
  }
  
  return(sample_cov)
}


#' @title Get default options for smooth covariates
#' @name get_default_smooth_options
#' @description Function to obtain the default options for the usage of smooth covariates
#' @param object an untrained \code{\link{MOFA}} object
#' @details The options are the following: \cr
#' \itemize{
#'  \item{\strong{scale_cov}:}
#'  \item{\strong{start_opt}:}
#'  \item{\strong{n_grid}:}
#'  \item{\strong{opt_freq}:}
#'  \item{\strong{sparseGP}:}
#'  \item{\strong{n_inducing}:}
#'  \item{\strong{warping}:}
#'  \item{\strong{warping_freq}:}
#'  \item{\strong{warping_ref}:}
#'  \item{\strong{warping_open_begin}:}
#'  \item{\strong{warping_open_end}:}
#'  \item{\strong{model_groups}:}
#' }
#' @return Returns a list with default options for the smooth covariate(s) functionality.
#' @importFrom utils modifyList
#' @export
#' @examples
get_default_smooth_options <- function(object) {
  
  smooth_options <- list(
    
    # Standard options
    scale_cov = TRUE,            # (logical) Scale covariates?
    start_opt = 20,              # (integer) First iteration to start the optimisation of GP hyperparameters
    n_grid = 20,                 # (integer) Number of points for the grid search in the optimisation of GP hyperparameters
    opt_freq = 10,               # (integer) Frequency of optimisation of GP hyperparameters
    model_groups = TRUE,         # (logical) model covariance structure across groups
    
    # sparse GP options
    sparseGP = FALSE,            # (logical) Use sparse GPs to speed up the optimisation of the GP parameters?
    n_inducing = 20,             # (integer) Number of inducing points
    
    # warping
    warping = FALSE,             # (logical) Activate warping functionality to align covariates between groups (requires a multi-group design)
    warping_freq = 20,           # (numeric) Warping: frequency of the optimisation
    warping_ref = groups_names(object)[[1]],          # (character) Warping: reference group
    warping_open_begin = TRUE,   # (logical) Warping: Allow for open beginning?
    warping_open_end = TRUE      # (logical) Warping: Allow for open ending?
    
  )
  
  # if smooth_options already exist, replace the default values but keep the additional ones
  if (length(object@smooth_options)>0)
    smooth_options <- modifyList(smooth_options, object@smooth_options)
  
  return(smooth_options)
}



#' @title Heatmap plot showing the group-group correlations per factor
#' @name plot_group_kernel
#' @description Heatmap plot showing the group-group correlations inferred by the model per factor
#' @param object a trained \code{\link{MOFA}} object.
#' @param factors character vector with the factors names, or numeric vector indicating the indices of the factors to use
#' @param groups character vector with the groups names, or numeric vector with the indices of the groups of samples to use, or "all" to use samples from all groups.
#' @param ... additional parameters that can be passed to  \code{pheatmap} 
#' @details The heatmap gives insight into the clustering of the patterns that factors display along the covariate in each group. 
#' A correlation of 1 indicates that the module caputred by a factor shows identical patterns across groups, a correlation of zero that it shows distinct patterns,
#' a negative correlation that the patterns go in opposite directions.
#' @return Returns a \code{ggplot,gg} object containing the heatmaps
#' @import pheatmap 
#' @import cowplot
#' @export

plot_group_kernel <- function(object, factors = "all", groups = "all", ...) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Define factors
  factors <- .check_and_get_factors(object, factors)

  # Define groups
  groups <- .check_and_get_groups(object, groups)
  
  # Get group kernels
  Kg <- get_group_kernel(object)
  
  hmlist <- lapply(factors, function(f){
    tmp <- Kg[[f]][groups,groups]
    # set breaks for heatmaps
    ncols <- 100
    seq_breaks <- c(seq(-1, 0, 1/ncols * 2), seq(0, 1, 1/ncols * 2)[-1])
    
    p <- pheatmap::pheatmap(tmp, color = rev(colorRampPalette((RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(ncols)), breaks = seq_breaks, silent = TRUE,...)
    
    p$gtable
  })
  # subset to groups
  
  p <- cowplot::plot_grid(plotlist = hmlist)

  return(p)
}



#' @title Barplot showing the smoothness per factor
#' @name plot_smoothness
#' @description Barplot indicating a smoothness score (between 0 (non-smooth) and 1 (smooth)) per factor
#' @param object a trained \code{\link{MOFA}} object.
#' @param factors character vector with the factors names, or numeric vector indicating the indices of the factors to use
#' @param color for the smooth part of the bar
#' @details The smoothness score is given by the scale parameter for the underlying Gaussian process of each factor.
#' @return Returns a \code{ggplot2} object
#' @import ggplot2
#' @export

plot_smoothness <- function(object, factors = "all", color = "cadetblue") {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Define factors
  factors <- .check_and_get_factors(object, factors)
  
  # Get scale parameters
  ss <- get_scales(object)[factors]
  df <- data.frame(factor = names(ss), smooth = ss, non_smooth = 1- ss)
  df <- gather(df, -factor, key = "smoothness", value = "value")
  gg_bar <- ggplot(df, aes(x= 1, y = value, fill = smoothness)) +
    geom_bar(stat="identity") +
    facet_wrap(~factor, nrow = 1) +
    theme_void() + coord_flip() +
    guides(fill=FALSE) + scale_fill_manual(values = c("non_smooth" = "gray", "smooth" = color)) +
    geom_text(x=1, y = 0.5, label = "smoothness", size = 3)

  return(gg_bar)
}


#' @title Barplot showing the sharedness per factor
#' @name plot_sharedness
#' @description Barplot indicating a sharedness score (between 0 (non-shared) and 1 (shared)) per factor
#' @param object a trained \code{\link{MOFA}} object.
#' @param factors character vector with the factors names, or numeric vector indicating the indices of the factors to use
#' @param color for the shared part of the bar
#' @details The sharedness score is calculated as the distance of the learnt group correlation matrix to the identity matrix
#'  in terms of the mean absolute distance on the off-diagonal elements.
#' @return Returns a \code{ggplot2} object
#' @import ggplot2
#' @export

plot_sharedness <- function(object, factors = "all", color = "#B8CF87") {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Define factors
  factors <- .check_and_get_factors(object, factors)
  
  # Get group kernels
  Kgs <- get_group_kernel(object)[factors]
  
  # Calculate distance
  idmat <- diag(1, ncol(Kgs[[1]]))
  gr <- sapply(Kgs, function(k) mean(abs(k - idmat)[lower.tri(idmat)]))
  
  # make plot
  df <- data.frame(factor = names(gr), group = gr, non_group = 1-gr)
  df <- gather(df, -factor, key = "sharedness", value = "value")
  df <- mutate(df, sharedness = factor(sharedness, levels = rev(c("group", "non_group"))))
  gg_bar <- ggplot(df, aes(x= 1, y=value, fill = sharedness)) + geom_bar(stat="identity") +
    facet_wrap(~factor, nrow = 1) +
    theme_void() + coord_flip() +
    guides(fill=FALSE) + scale_fill_manual(values = c("non_group" = "gray", "group" = color)) +
    geom_text(x=1, y = 0.5, label = "sharedness", size = 3)
  
  return(gg_bar)
}
#' 
#' #' @title Plot interpolated factors versus covariate (1-dimensional)
#' #' @name plot_sharedness
#' #' @description Barplot indicating a sharedness score (between 0 (non-shared) and 1 (shared)) per factor
#' #' @param object a trained \code{\link{MOFA}} object.
#' #' @param covariate covariate to use for plotting
#' #' @param factors character vector with the factors names, or numeric vector indicating the indices of the factors to use
#' #' @param only_mean show onyl mean or include uncertainties?
#' #' @param show_observed include observed factor values as dots on the plot
#' #' @details to be filled
#' #' @return Returns a \code{ggplot2} object
#' #' @import ggplot2
#' #' @export
#' 
#' plot_interpolation_vs_covariate <- function(object, covariate = 1, factors = "all", only_mean = FALSE, show_observed = TRUE){
#'   
#'   # Sanity checks
#'   if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
#'   
#'   # get and check covariate
#'   covariate <- .check_and_get_covariates(object, covariate)
#'   
#'   # get interpolated factors
#'   df <- get_interpolated_factors(object, as.data.frame = TRUE)
#'   
#'   # calculate ribbon borders
#'   if(!only_mean) {
#'     df %<>% mutate(sd = sqrt(variance), ymin = mean -1.96 * sd, ymax = mean + 1.96 * sd)
#'   }
#' 
#'   if(show_observed) {
#'     # add the factor values of the observed time point  to the plot
#'     df_observed <- plot_factors_vs_cov(object, covariate = covariate, return_data = TRUE)
#'   }
#'   
#'   gg_interpol <- ggplot(df, aes_string(x="covariate", y = "mean", col = "group")) +
#'     geom_line(aes(y=mean,  col = group)) +
#'     facet_wrap(~ factor) + theme_classic()
#'   
#'   if(show_observed) {
#'     gg_interpol <- gg_interpol + geom_point(data = df_observed, aes(x= covariate_value,
#'                                                                   y = value, col = group), size = 1) 
#'   }
#'   if(!only_mean) {
#'     gg_interpol <- gg_interpol + geom_ribbon(aes(ymin=ymin, ymax = ymax, fill = group),
#'                                              alpha = .2, col = "gray", size = 0.1)   
#'   }
#'   
#'   gg_interpol
#' }




#' @title Scatterplots of feature values against sample covariates
#' @name plot_data_scatter_vs_cov
#' @description Function to do a scatterplot of features against sample covariate values.
#' @param object a \code{\link{MOFA}} object.
#' @param covariate string with the covariate name or a samples_metadata column, or an integer with the index of the covariate
#' @param factor string with the factor name, or an integer with the index of the factor to take top features from
#' @param view string with the view name, or an integer with the index of the view. Default is the first view.
#' @param groups groups to plot. Default is "all".
#' @param features if an integer (default), the total number of features to plot (given by highest weights). If a character vector, a set of manually-defined features.
#' @param sign can be 'positive', 'negative' or 'all' (default) to show only features with highest positive, negative or all weights, respectively.
#' @param color_by specifies groups or values (either discrete or continuous) used to color the dots (samples). This can be either: 
#' \itemize{
#' \item the string "group": dots are coloured with respect to their predefined groups.
#' \item a character giving the name of a feature that is present in the input data 
#' \item a character giving the same of a column in the sample metadata slot
#' \item a vector of the same length as the number of samples specifying the value for each sample. 
#' \item a dataframe with two columns: "sample" and "color"
#' }
#' @param shape_by specifies groups or values (only discrete) used to shape the dots (samples). This can be either: 
#' \itemize{
#' \item the string "group": dots are shaped with respect to their predefined groups.
#' \item a character giving the name of a feature that is present in the input data 
#' \item a character giving the same of a column in the sample metadata slot
#' \item a vector of the same length as the number of samples specifying the value for each sample. 
#' \item a dataframe with two columns: "sample" and "shape"
#' }
#' @param legend logical indicating whether to add a legend
#' @param dot_size numeric indicating dot size (default is 5).
#' @param text_size numeric indicating text size (default is 5).
#' @param stroke numeric indicating the stroke size (the black border around the dots, default is NULL, infered automatically).
#' @param alpha numeric indicating dot transparency (default is 1).
#' @param add_lm logical indicating whether to add a linear regression line for each plot
#' @param lm_per_group logical indicating whether to add a linear regression line separately for each group
#' @param imputed logical indicating whether to include imputed measurements
#' @param return_data logical indicating whether to return a data frame instead of a plot
#' @details One of the first steps for the annotation of factors is to visualise the weights using \code{\link{plot_weights}} or \code{\link{plot_top_weights}}
#' and inspect the relationshio of the factor to the covariate(s) using  \code{\link{plot_factors_vs_cov}}.
#' However, one might also be interested in visualising the direct relationship between features and covariate(s), rather than looking at "abstract" weights and
#' possibly look at the interpolated and extrapolated values by setting imputed to True.
#' @import ggplot2
# #' @importFrom ggpubr stat_cor
#' @importFrom dplyr left_join
#' @importFrom utils tail
#' @importFrom stats quantile
#' @export
plot_data_scatter_vs_cov <- function(object, covariate = 1, factor = 1, view = 1, groups = "all", features = 10, sign = "all",
                              color_by = "group", legend = TRUE, alpha = 1, shape_by = NULL, stroke = NULL,
                              dot_size = 2.5, text_size = NULL, add_lm = FALSE, lm_per_group = FALSE, imputed = FALSE, return_data = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(length(factor)==1)
  stopifnot(length(covariate)==1)
  stopifnot(length(view)==1)
  if (lm_per_group) add_lm = TRUE
  
  # Define views, factors and groups
  groups <- .check_and_get_groups(object, groups)
  factor <- .check_and_get_factors(object, factor)
  view <- .check_and_get_views(object, view)
  
  # Check and fetch covariates
  if(covariate %in% colnames(samples_metadata(object))){
    covari <- samples_metadata(object)[,c("sample", covariate)]
  } else {
    covariate <- .check_and_get_covariates(object, covariates = covariate)
    covari <- get_covariates(object, covariates = covariate, as.data.frame = TRUE)
    covari <- covari[,c("sample","value")]
  }
  colnames(covari) <- c("sample","x")
  
  # Collect relevant data
  N <- get_dimensions(object)[["N"]]
  W <- get_weights(object)[[view]][,factor]
  
  # Get features
  if (sign=="all") {
    W <- abs(W)
  } else if (sign=="positive") {
    W <- W[W>0]
  } else if (sign=="negative") {
    W <- W[W<0]
  }
  
  if (is(features, "numeric")) {
    if (length(features) == 1) {
      features <- names(tail(sort(abs(W)), n=features))
    } else {
      features <- names(sort(-abs(W))[features])
    }
    stopifnot(all(features %in% features_names(object)[[view]]))  
  } else if (is(features, "character")) {
    stopifnot(all(features %in% features_names(object)[[view]]))
  } else {
    stop("Features need to be either a numeric or character vector")
  }

  # Set group/color/shape
  if (length(color_by)==1 && is.character(color_by)) color_name <- color_by
  if (length(shape_by)==1 && is.character(shape_by)) shape_name <- shape_by
  color_by <- .set_colorby(object, color_by)
  shape_by <- .set_shapeby(object, shape_by)
  
  # Merge factor values with color and shape information
  df1 <- merge(covari, color_by, by="sample")
  df1 <- merge(df1, shape_by, by="sample")
  
  # Create data frame 
  foo <- list(features); names(foo) <- view
  if (imputed) {
    df2 <- get_imputed_data(object, groups = groups, views = view, features = foo, as.data.frame = TRUE)
  } else {
    df2 <- get_data(object, groups = groups, features = foo, as.data.frame = TRUE)
  }
  
  df2$sample <- as.character(df2$sample)
  df <- left_join(df1, df2, by = "sample")
  
  # (Q) Remove samples with missing values in Factor values
  df <- df[!is.na(df$value),]
  
  if(return_data){
    return(df)
  }
  
  # Set stroke
  if (is.null(stroke)) {
    stroke <- .select_stroke(N=length(unique(df$sample)))
  }
  
  # Set Pearson text size
  if (add_lm && is.null(text_size)) {
    text_size <- .select_pearson_text_size(N=length(unique(df$feature)))
  }
  
  # Set axis text size
  axis.text.size <- .select_axis.text.size(N=length(unique(df$feature)))
  
  # Generate plot
  p <- ggplot(df, aes_string(x = "x", y = "value")) + 
    geom_point(aes_string(fill = "color_by", shape = "shape_by"), colour = "black", size = dot_size, stroke = stroke, alpha = alpha) +
    labs(x=covariate, y="") +
    facet_wrap(~feature, scales="free_y") +
    theme_classic() + 
    theme(
      axis.text = element_text(size = rel(axis.text.size), color = "black"), 
      axis.title = element_text(size = rel(1.0), color="black")
    )
  
  # Add linear regression line
  if (add_lm) {
    if (lm_per_group && length(groups)>1) {
      p <- p +
        stat_smooth(formula=y~x, aes_string(color="group"), method="lm", alpha=0.4) +
        ggpubr::stat_cor(aes_string(color="group", label = "..r.label.."), method = "pearson", label.sep="\n", output.type = "latex", size = text_size)# +
      # guides(color = FALSE)
    } else {
      p <- p +
        stat_smooth(formula=y~x, method="lm", color="grey", fill="grey", alpha=0.4) +
        ggpubr::stat_cor(method = "pearson", label.sep="\n", output.type = "latex", size = text_size, color = "black")
    }
  }
  
  # Add legend
  p <- .add_legend(p, df, legend, color_name, shape_name)
  
  return(p)
}


#' @title Scatterplots of a factor's values againt the sample covariates
#' @name plot_factors_vs_cov
#' @description  Scatterplots of a factor's values againt the sample covariates
#' @param object a trained \code{\link{MOFA}} object.
#' @param factors character or numeric specifying the factor(s) to plot, default is "all"
#' @param covariate specifies sample covariate to plot against:
#' (1) a character giving the name of a column present in the sample covariates or sample metadata.
#' (2) a character giving the name of a feature present in the training data.
#' (3) a vector of the same length as the number of samples specifying continuous numeric values per sample.
#' Default is the first sample covariates in covariates slot
#' @param show_missing logical indicating whether to include samples for which \code{shape_by} or \code{color_by} is missing
#' @param scale logical indicating whether to scale factor values.
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
#' @param original logical indicating whether to use unscaled covariates (only if used with GP prior)
#' @param return_data logical indicating whether to return the data frame to plot instead of plotting
#' @param show_variance logical indicating whether to show the marginal variance of inferred factor values 
#' (only relevant for 1-dimensional covariates)
#' @details To investigate the factors pattern along the covariates (such as time or a spatial coordinate) 
#' this function an be used to plot a scatterplot of the factor againt the values of each covariate
#' @return Returns a \code{ggplot2} object
#' @import ggplot2 dplyr
#' @importFrom stats complete.cases
#' @importFrom tidyr spread
#' @importFrom magrittr %>% set_colnames
#' @importFrom ggbeeswarm geom_quasirandom
#' @export
plot_factors_vs_cov <- function(object, factors = "all", covariates = NULL, show_missing = TRUE, scale = FALSE,
                                color_by = NULL, shape_by = NULL, color_name = NULL, shape_name = NULL,
                                dot_size = 1.5, alpha = 1, stroke = NULL, legend = TRUE, return_data = FALSE, show_variance = FALSE) {
  
  # TO-DO: ORIGINAL = FALSE
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Define covariates
  if (is.null(covariates)) {
    if (any(object@dimensions[["C"]] < 1, is.null(object@covariates)))  
      stop("No covariates found in object. Please specify one.")
    covariates <- covariates_names(object)
  }
  
  # Get factors
  factors <- .check_and_get_factors(object, factors)
  Z <- get_factors(object, factors=factors, as.data.frame=TRUE)
  
  # Remove samples with missing values
  Z <- Z[complete.cases(Z),]
  
  # Get covariates (TO-DO: AD DWARP OPTIONS)
  df <- get_covariates(object, covariates, as.data.frame = TRUE) %>%
    merge(Z, by="sample", suffixes = c(".covariate",".factor"))
  
  # Remember color_name and shape_name if not provided
  if (!is.null(color_by) && (length(color_by) == 1) && is.null(color_name))
    color_name <- color_by
  if (!is.null(shape_by) && (length(shape_by) == 1) && is.null(shape_name))
    shape_name <- shape_by
  
  # Set color and shape
  color_by <- .set_colorby(object, color_by)
  shape_by <- .set_shapeby(object, shape_by )
  
  # Merge factor values with color and shape information
  df <- df %>%
    merge(color_by, by="sample") %>%
    merge(shape_by, by="sample") %>%
    mutate(shape_by = as.character(shape_by))
  
  # Remove missing values
  if (!show_missing) df <- filter(df, !is.na(color_by) && !is.na(shape_by))
  
  # Return data if requested instead of plotting
  if (return_data) return(df)
  
  # Set stroke
  if (is.null(stroke)) stroke <- .select_stroke(N=length(unique(df$sample)))
  
  # Select 1D or 2D plots
  if (length(covariates) == 1) {
    
    # Include marginal variance
    if (show_variance) {
      if("E2" %in% names(object@expectations$Z)){
        ZZ = object@expectations$Z$E2
        ZZ <- reshape2::melt(ZZ, na.rm=T)
        colnames(ZZ) <- c("sample", "factor", "E2")
        df <- left_join(df, ZZ, by = c("sample", "factor"))
        df <- mutate(df, var = E2 - value^2)
      } else {
        show_variance <- FALSE
        warning("No second moments saved in the trained model - variance can not be shown.")
      }
    }
    p <- .plot_factors_vs_cov_1d(df,
            color_name = color_name,
            shape_name = shape_name,
            scale = scale, 
            dot_size = dot_size, 
            alpha = alpha, 
            stroke = stroke,
            show_variance = show_variance,
            legend = legend
          ) 
  } else if (length(covariates) == 2) {
    p <- .plot_factors_vs_cov_2d(df,
           color_name = color_name,
           shape_name = shape_name,
           scale = scale, 
           dot_size = dot_size, 
           alpha = alpha, 
           stroke = stroke,
           legend = legend
          )
  } else {
    stop("too many covariates provided")
  }
  
  return(p)
}


.plot_factors_vs_cov_1d <- function(df, color_name = "", shape_name = "", scale = FALSE, dot_size = 1.5, alpha = 1, stroke = 1, show_variance = FALSE, legend = TRUE) {
  
  # Sanity checks
  stopifnot(length(unique(df$covariate))==1)
  
  
  # Scale values from 0 to 1
  if (scale) {
    df <- df %>% 
      group_by(factor) %>%
      mutate(value_scaled = value.factor/max(abs(value.factor)))
    if(show_variance) df <- mutate(df, var = var/(max(abs(value.factor))^2))
    df <- df %>% 
      mutate(value.factor = value_scaled) %>%
      select(-value_scaled) %>%
      ungroup
  }
  
  # Generate plot
  p <- ggplot(df, aes(x=value.covariate, y=value.factor)) + 
    # geom_point(aes_string(color = "color_by", shape = "shape_by"), size=dot_size, alpha=alpha) +
    geom_point(aes_string(fill = "color_by", shape = "shape_by"), colour="black", stroke = stroke, size=dot_size, alpha=alpha) +
    facet_grid(~ factor) +
    theme_classic() +
    theme(
      axis.text = element_text(size = rel(0.9), color = "black"), 
      axis.title = element_text(size = rel(1.2), color = "black"), 
      axis.line = element_line(color = "black", size = 0.5), 
      axis.ticks = element_line(color = "black", size = 0.5)
    )
  
  if (show_variance){
    p <- p + geom_errorbar(aes(ymin = value - sqrt(var)*1.96, ymax =value + sqrt(var)*1.96), col = "red", alpha = 0.7)
  }
  
  p <- .add_legend(p, df, legend, color_name, shape_name)
  
  return(p)
}

.plot_factors_vs_cov_2d <- function(df, color_name = "", shape_name = "", scale = FALSE, dot_size = 1.5, alpha = 1, stroke = 1, legend = TRUE) {
  
  # Sanity checks
  stopifnot(length(unique(df$covariate))==2)
  
  # pivot covariate values
  covariates_dt <- df %>%
    tidyr::pivot_wider(names_from="covariate", values_from="value.covariate") 
  
  covariates.names <- c(colnames(covariates_dt)[ncol(covariates_dt)-1], colnames(covariates_dt)[ncol(covariates_dt)])
  
  # Scale factor values from 0 to 1
  if (scale) {
    covariates_dt <- covariates_dt %>%
      group_by(factor) %>%
      mutate(value.factor = value.factor/max(abs(value.factor))) %>%
      ungroup
  }
  
  covariates_dt <- mutate(covariates_dt, color_by = value.factor) # for compatibility with .add_legend
  # Generate plot
  p <- ggplot(covariates_dt, aes_string(x=covariates.names[1], y=covariates.names[2], fill = "color_by")) + 
    geom_tile() + 
    facet_grid( ~ factor) + 
    theme_bw() +
    theme(
      axis.text = element_text(size = rel(0.9), color = "black"), 
      axis.title = element_text(size = rel(1.2), color = "black"), 
      axis.line = element_line(color = "black", size = 0.5), 
      axis.ticks = element_line(color = "black", size = 0.5)
    ) + coord_fixed()
  
  # p <- ggplot(covariates_dt, aes_string(x=covariates.names[1], y=covariates.names[2])) + 
  #   geom_point(aes_string(fill = "value.factor"), shape=21, colour = "black", size = dot_size, stroke = stroke, alpha = alpha) + 
  #   # geom_point( col = "grey", alpha =0.05)+
  #   facet_grid( ~ factor) + 
  #   theme_bw() +
  #   theme(
  #     axis.text = element_text(size = rel(0.9), color = "black"), 
  #     axis.title = element_text(size = rel(1.0), color = "black"), 
  #     axis.line = element_line(color = "black", size = 0.5), 
  #     axis.ticks = element_line(color = "black", size = 0.5)
  #   )
  
  p <- .add_legend(p, covariates_dt, legend, color_name, shape_name)
  
  return(p)
}

