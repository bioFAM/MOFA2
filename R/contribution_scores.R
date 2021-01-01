calculate_contribution_scores <- function(object, views = "all", groups = "all", factors = "all", scale = TRUE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (any(object@model_options$likelihoods!="gaussian"))
    stop("Not possible to compute contribution scores when using non-gaussian likelihoods.")
  
  # Define factors, views and groups
  views  <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)
  factors <- .check_and_get_factors(object, factors)
  
  # fetch variance explained values
  r2.per.factor <- get_variance_explained(object, factors=factors, views = views, groups = groups)[["r2_per_factor"]]
  r2.per.view <- get_variance_explained(object, factors=factors, views = views, groups = groups)[["r2_total"]]
  
  # calculate relative variance explained values (each view goes from 0 to 1)
  r2.per.factor <- lapply(1:length(groups), function(g) sweep(r2_per_factor[[1]], 2, r2.per.view,"/"))
  
  # fetch factors
  Z <- get_factors(object, groups=groups, factors=factors)
  
  # scale factors
  Z <- lapply(Z, function(z) apply(z, 2, function(x) abs(x)/max(abs(x)) ))
  
  # compute contribution scores per sample
  contribution_scores <- lapply(1:length(groups), function(g) Z[[g]] %*% r2.per.factor[[g]])
  
  # Add contributions to the sample metadata
  contribution_scores <- do.call("rbind",contribution_scores)
  for (i in colnames(contribution_scores)) {
    object <- .add_column_to_metadata(object, contribution_scores[,i], paste0(i,"_contribution"))
  }
  
  # Store in cache
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