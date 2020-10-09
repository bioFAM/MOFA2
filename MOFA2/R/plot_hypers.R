
########################################################################
## Functions to visualise the hyperparameters of the latent processes ##
########################################################################

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
