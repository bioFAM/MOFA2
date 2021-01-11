#' @title Calculate contribution scores for each data modality in each sample
#' @description  This function takes a trained MOFA model as input and calculates, **for each sample** a contribution score
#' that quantifies how much each data modality is by the latent space. In other words, how much each data modality "contributes to the latent space.
#' @name calculate_contribution_scores
#' @param object a \code{\link{MOFA}} object.
#' @param views character vector with the view names, or numeric vector with view indexes. Default is 'all'
#' @param groups character vector with the group names, or numeric vector with group indexes. Default is 'all'
#' @param factors character vector with the factor names, or numeric vector with the factor indexes. Default is 'all'
#' @param scale logical indicating whether to scale the sample-wise variance explained values by the total amount of variance explained per view. 
#' This is useful when the collection of all factors explain different amounts of variance for each data modality (check with `plot_variance_explained(..., plot_total=T)`)
#' @details TO-FILL
#' @return adds the contribution scores to the metadata slot (see `samples_metadata(MOFAobject)`) and it also returns a data.frame with the contribution score for each sample and data modality
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' r2 <- calculate_contribution_scores(model, scale = FALSE)
#' r2 <- calculate_contribution_scores(model, scale = TRUE)
#'
calculate_contribution_scores <- function(object, views = "all", groups = "all", factors = "all", scale = TRUE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (any(object@model_options$likelihoods!="gaussian"))
    stop("Not possible to compute contribution scores when using non-gaussian likelihoods.")

  # Define factors, views and groups
  views  <- .check_and_get_views(object, views)
  if (length(views)<2) stop("contribution scores only make sense when having at least 2 views")
  groups <- .check_and_get_groups(object, groups)
  factors <- .check_and_get_factors(object, factors)
  if (length(factors)<2) stop("contribution scores only make sense when having at least 2 factors")
  
  # fetch variance explained values
  r2.per.sample <- calculate_variance_explained_per_sample(object, factors=factors, views = views, groups = groups)
  
  # scale the variance explained values to the total amount of variance explained per view
  if (scale) {
    r2.per.view <- get_variance_explained(object, factors=factors, views = views, groups = groups)[["r2_total"]]
    r2.per.sample <- lapply(1:length(groups), function(g) sweep(r2.per.sample[[g]], 2, r2.per.view[[g]],"/"))
  }
  
  # concatenate groups
  r2.per.sample <- do.call("rbind",r2.per.sample)
  
  # Calculate the fraction of (relative) variance explained for each data modality in each cell -> the contribution score
  contribution_scores <- r2.per.sample / rowSums(r2.per.sample)
  
  # Add contribution scores to the sample metadata
  for (i in colnames(contribution_scores)) {
    object <- .add_column_to_metadata(object, contribution_scores[,i], paste0(i,"_contribution"))
  }
  # Add contribution scores to the cache
  object@cache[["contribution_scores"]] <- contribution_scores
  
  
  return(object)
  
}


get_contribution_scores <- function(object, groups = "all", views = "all", factors = "all", 
                                   as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get factors and groups
  groups <- .check_and_get_groups(object, groups)
  views <- .check_and_get_views(object, views)
  factors <- .check_and_get_factors(object, factors)
  
  # Fetch
  if (.hasSlot(object, "cache") && ("contribution_scores" %in% names(object@cache))) {
    scores_list <- object@cache[["contribution_scores"]]
  } else {
    scores_list <- calculate_contribution_scores(object, factors = factors, views = views, groups = groups)
  }
  
  # Convert to data.frame format
  if (as.data.frame) {
    scores <- reshape2::melt( do.call("rbind",scores_list) )
    colnames(scores) <- c("sample", "view", "value")
  } else {
    scores <- scores_list
  }
  
  return(scores)
  
}

plot_contribution_scores <- function(object, samples = "all", group_by = NULL, return_data = FALSE, ...) {

  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # (TO-DO) get samples
  
  # get contribution scores
  scores <- get_contribution_scores(object, as.data.frame = TRUE, ...)
  
  # TO-DO: CHECK THAT GROUP IS A CHARACTER/FACTOR
  # individual samples
  if (is.null(group_by)) {
    
    to.plot <- scores
    if (return_data) return(to.plot)
    p <- ggplot(to.plot, aes_string(x="view", y="value")) +
      geom_bar(aes(fill=view), stat="identity", color="black") +
      facet_wrap(~sample) +
      labs(x="", y="Contribution score") +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank()
      )
    return(p)
    
  # group samples
  } else {
    
    to.plot <- merge(scores, object@samples_metadata[,c("sample",group_by)], by="sample")
    if (return_data) return(to.plot)
    p <- ggplot(to.plot, aes_string(x="view", y="value")) +
      geom_boxplot(aes(fill=view)) +
      facet_wrap(as.formula(paste("~", group_by))) +
      labs(x="", y="Contribution score") +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "top",
        legend.title = element_blank()
      )
    return(p)
  }
}