
################################################
## Functions to compare different MOFAs ##
################################################


#' @title Plot the robustness of the latent factors across diferent trials
#' @name compare_factors
#' @description Different \code{\link{MOFA}} objects are compared in terms of correlation between their factors.
#' @param models a list with \code{\link{MOFA}} objects.
#' @param ... extra arguments passed to pheatmap
#' @details If assessing model robustness across trials, the output should look like a block diagonal matrix, 
#' suggesting that all factors are robustly detected in all model instances.
#' @return Plots a heatmap of the Pearson correlation between latent factors across all input models.
#' @importFrom stats cor
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @export
compare_factors <- function(models, ...) {
  
  # Sanity checks
  if(!is.list(models))
    stop("'models' has to be a list")
  if (!all(sapply(models, function (l) class(l)=="MOFA")))
    stop("Each element of the the list 'models' has to be an instance of MOFA")

  # Give generic names if no names present
  if(is.null(names(models))) names(models) <- paste("model", 1: length(models), sep="")

  # Get latent factors
  LFs <- lapply(seq_along(models), function(i){
    do.call(rbind, get_factors(models[[i]]))
    })
  
  # Rename factors
  for (i in seq_along(LFs))
    colnames(LFs[[i]]) <- paste(names(models)[i], colnames(LFs[[i]]), sep="_")
  
  # Sanity checks
  if (is.null(Reduce(intersect,lapply(LFs, rownames))))
    stop("No common samples in all models for comparison")

  # calculate correlation between factors across models
  corLFs <- cor(Reduce(cbind, LFs), use="complete.obs")
  corLFs[is.na(corLFs)] <- 0
  corLFs <- abs(corLFs)

  # Plot heatmap
  breaksList <- seq(0,1, by=0.01)
  colors <- colorRampPalette(c("white",RColorBrewer::brewer.pal(9,name="YlOrRd")))(length(breaksList))
  pheatmap(corLFs, color = colors, breaks = breaksList, ...)
}



#' @title Compare different trained \code{\link{MOFA}} objects in terms of the final value of the ELBO statistics and number of inferred factors
#' @name compare_elbo
#' @description Different objects of \code{\link{MOFA}} are compared in terms of the final value of the ELBO statistics.
#' For model selection the model with the highest ELBO value is selected.
#' @param models a list containing \code{\link{MOFA}} objects.
#' @export
compare_elbo <- function(models) {
  
  # Sanity checks
  if(!is.list(models))
    stop("'models' has to be a list")
  if (!all(sapply(models, function (l) class(l)=="MOFA")))
    stop("Each element of the the list 'models' has to be an instance of MOFA")
  
  # Give generic names if no names present
  if (is.null(names(models))) names(models) <- paste0("model_", seq_along(models))
  
  # Get ELBO values
  elbo_vals <- sapply(models, get_elbo)
  
  # Generate plot
  df <- data.frame(
    ELBO = elbo_vals, 
    model = names(models)
  )
  df$logELBO <- -log2(-df$ELBO)
  
  gg <- ggplot(df, aes_string(x="model", y="logELBO")) + 
    geom_bar(stat="identity", color="black", fill="grey70") +
    labs(x="", y="log Evidence Lower Bound (ELBO)") +
    theme_classic()
  
  gg <- gg + theme(
    axis.text.x=element_text(angle=60, vjust=1, hjust=1)
  )
  
  return(gg)
}



#' @title Select a model from a list of trained \code{\link{MOFA}} objects based on the best ELBO value
#' @name select_model
#' @description Different objects of \code{\link{MOFA}} are compared in terms of the final value of the ELBO statistics
#' and the model with the highest ELBO value is selected.
#' @param models a list containing \code{\link{MOFA}} objects.
#' @export
select_model <- function(models) {
  # Sanity checks
  if(!is.list(models))
    stop("'models' has to be a list")
  if (!all(sapply(models, function (l) class(l)=="MOFA")))
    stop("Each element of the the list 'models' has to be an instance of MOFA")

  elbo_vals <- sapply(models, get_elbo)
  if(plotit) compare_models(models)
  models[[which.max(elbo_vals)]]
}
