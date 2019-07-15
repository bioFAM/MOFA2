
################################################
## Functions to compare different BioFAModels ##
################################################


#' @title Plot the robustness of the latent factors across diferent trials
#' @name compare_factors
#' @description Different objects of \code{\link{BioFAModel}} are compared in terms of correlation between
#' their latent factors. The correlation is calculated only on those samples which are present in all models.
#' Ideally, the output should look like a block diagonal matrix, suggesting that all detected factors are robust under different initialisations.
#' If not, it suggests that some factors are weak and not captured by all models.
#' @param models a list containing \code{\link{BioFAModel}} objects.
#' @param comparison tye of comparison, either 'pairwise' or 'all'
#' @param ... extra arguments passed to pheatmap
#' @details TO-FILL
#' @return Plots a heatmap of correlation of Latent Factors in all models when 'comparison' is 'all'.
#' Otherwise, for each pair of models, a seperate heatmap is produced comparing one model againt the other.
#' The corresponding correlation matrix or list or pairwise correlation matrices is returned
#' @references fill this
#' @importFrom stats cor
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @export
compare_factors <- function(models, ...) {
  
  # Sanity checks
  if(!is.list(models))
    stop("'models' has to be a list")
  if (!all(sapply(models, function (l) class(l)=="BioFAModel")))
    stop("Each element of the the list 'models' has to be an instance of BioFAModel")

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



#' @title Compare different trained \code{\link{BioFAModel}} objects in terms of the final value of the ELBO statistics and number of inferred factors
#' @name compare_elbo
#' @description Different objects of \code{\link{BioFAModel}} are compared in terms of the final value of the ELBO statistics.
#' For model selection the model with the highest ELBO value is selected.
#' @param models a list containing \code{\link{BioFAModel}} objects.
#' @export
compare_elbo <- function(models) {
  
  # Sanity checks
  if(!is.list(models))
    stop("'models' has to be a list")
  if (!all(sapply(models, function (l) class(l)=="BioFAModel")))
    stop("Each element of the the list 'models' has to be an instance of BioFAModel")
  
  # Give generic names if no names present
  if (is.null(names(models))) names(models) <- paste0("model_", seq_along(models))
  
  # Get ELBO values
  elbo_vals <- sapply(models, get_elbo)
  
  # Generate plot
  df <- data.frame(
    ELBO = elbo_vals, 
    model = names(models)
  )
  
  gg <- ggplot(df, aes(x=model, y=-log2(-ELBO))) + 
    geom_bar(stat="identity", color="black", fill="grey70") +
    labs(x="", y="log Evidence Lower Bound (ELBO)") +
    theme_classic()
  
  gg <- gg + theme(axis.text.x=element_text(angle=60, vjust=1, hjust=1))
  
  return(gg)
}



#' @title Select a model from a list of trained \code{\link{BioFAModel}} objects based on the best ELBO value
#' @name select_model
#' @description Different objects of \code{\link{BioFAModel}} are compared in terms of the final value of the ELBO statistics
#' and the model with the highest ELBO value is selected.
#' @param models a list containing \code{\link{BioFAModel}} objects.
#' @export
select_model <- function(models, plotit = TRUE) {
  # Sanity checks
  if(!is.list(models))
    stop("'models' has to be a list")
  if (!all(sapply(models, function (l) class(l)=="BioFAModel")))
    stop("Each element of the the list 'models' has to be an instance of BioFAModel")

  elbo_vals <- sapply(models, get_elbo)
  if(plotit) compare_models(models)
  models[[which.max(elbo_vals)]]
}
