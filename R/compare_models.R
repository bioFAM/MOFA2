
################################################
## Functions to compare different MOFA models ##
################################################


#' @title Plot the correlation of factors between different models
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
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model1 <- load_model(file)
#' model2 <- load_model(file)
#' 
#' # Compare factors between models
#' compare_factors(list(model1,model2))
compare_factors <- function(models, ...) {
  
  # Sanity checks
  if(!is.list(models))
    stop("'models' has to be a list")
  if (!all(sapply(models, function (l) is(l, "MOFA"))))
    stop("Each element of the the list 'models' has to be an instance of MOFA")

  # Give generic names if no names present
  if(is.null(names(models))) names(models) <- paste("model", seq_len(length(models)), sep="")

  # Get latent factors
  LFs <- lapply(seq_along(models), function(i){
    do.call(rbind, get_factors(models[[i]]))
  })
  
  # Sanity checks
  if (is.null(Reduce(intersect,lapply(LFs, rownames))))
    stop("No common samples in all models for comparison")

  # Align samples between models
  samples_names <- Reduce(intersect, lapply(LFs, rownames))
  LFs <- lapply(LFs, function(z) {
    z[samples_names,,drop=FALSE]
  })
  
  # Rename factors
  for (i in seq_along(LFs))
    colnames(LFs[[i]]) <- paste(names(models)[i], colnames(LFs[[i]]), sep="_")

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
#' @param log logical indicating whether to plot the log of the ELBO.
#' @param return_data logical indicating whether to return a data.frame with the ELBO values per model
#' @return A \code{\link{ggplot}} object or the underlying data.frame if return_data is TRUE
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model1 <- load_model(file)
#' model2 <- load_model(file)
#' 
#' # Compare ELBO between models
#' \dontrun{compare_elbo(list(model1,model2))}
compare_elbo <- function(models, log = FALSE, return_data = FALSE) {
  
  # Sanity checks
  if(!is.list(models))
    stop("'models' has to be a list")
  if (!all(sapply(models, function (l) is(l, "MOFA"))))
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
  
  
  # take the log
  if (log) {
    message("Plotting the log2 of the negative of the ELBO (the higher the better)")
    df$ELBO <- log2(-df$ELBO)
  }
  
  if (all(df$ELBO<0)) {
    df$ELBO <- abs(df$ELBO)
    message("Plotting the absolute value of the ELBO for every model (the smaller the better)")
} else {
    message("Plotting the ELBO for every model (the higher the better)")
  }
  
  # return data
  if (return_data) return(df)
  
  gg <- ggplot(df, aes_string(x="model", y="ELBO")) + 
    geom_bar(stat="identity", color="black", fill="grey70") +
    labs(x="", y="Evidence Lower Bound (ELBO)") +
    theme_classic()
  
  return(gg)
}



#' @title Select a model from a list of trained \code{\link{MOFA}} objects based on the best ELBO value
#' @name select_model
#' @description Different objects of \code{\link{MOFA}} are compared in terms of the final value of the ELBO statistics
#' and the model with the highest ELBO value is selected.
#' @param models a list containing \code{\link{MOFA}} objects.
#' @param plot boolean indicating whether to show a plot of the ELBO for each model instance
#' @return A \code{\link{MOFA}} object
#' @export
select_model <- function(models, plot = FALSE) {
  # Sanity checks
  if(!is.list(models))
    stop("'models' has to be a list")
  if (!all(sapply(models, function (l) is(l, "MOFA"))))
    stop("Each element of the the list 'models' has to be an instance of MOFA")

  elbo_vals <- sapply(models, get_elbo)
  if(plot) compare_elbo(models)
  models[[which.max(elbo_vals)]]
}
