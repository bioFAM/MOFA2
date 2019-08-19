
#######################################################
## Functions to perform imputation of missing values ##
#######################################################

#' @title Impute missing values from a fitted MOFA
#' @name impute
#' @description This function uses the latent factors and the loadings to impute missing values.
#' @param object a \code{\link{MOFA}} object.
#' @param views character vector with the view name(s), or numeric vector with view index(es).
#' @param groups character vector with the group name(s), or numeric vector with group index(es).
#' @param factors character vector with the factor names, or numeric vector with the factor index(es).
#' @param type type of imputation.
#' \itemize{
#' \item \strong{response}:{gives mean for gaussian and poisson and probabilities for bernoulli.}
#' \item \strong{link}: {gives the linear predictions.}
#' \item \strong{inRange}: {rounds the fitted values from "terms" for integer-valued distributions to the next integer (default).}
#' }
#' @details matrix factorization models generate a denoised and condensed low-dimensional representation of the data which capture the main sources of heterogeneity of the data.
#' These representation can be used to do predictions using the equation \code{Y = WX}. For more details read the supplementary methods of the manuscript. \cr
#' This method fills the \code{ImputedData} slot by replacing the missing values in the input data with the model predictions.
#' @export
impute <- function(object, views = "all", groups = "all", factors = "all", type = c("inRange", "response", "link")) {

  # Get views and groups
  views  <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)

  # Select imputation type
  type <- match.arg(type)

  # Do predictions
  pred <- predict(object, views=views, factors=factors, type=type)

  # Replace NAs with predicted values
  # TO-DO: HOW TO REPLACE  VALUES WITH DELAYEDARRAYS??
  imputed <- get_data(object, views=views, groups=groups, add_intercept = FALSE)
  for (m in views) {
    for (g in groups) {
      non_observed <- which(is.na(imputed[[m]][[g]]), arr.ind = T)
      imputed[[m]][[g]][non_observed] <- pred[[m]][[g]][non_observed]
    }
  }

  # Save imputed data in the corresponding slot
  object@imputed_data <- imputed

  return(object)
}

# 
# impute2 <- function(object, views = "all", groups = "all", factors = "all") {
#
#   # Get views and groups
#   views  <- .check_and_get_views(object, views)
#   groups <- .check_and_get_groups(object, groups)
#
#
#   # Do predictions by sampling from the variational distributions
#   # TODO fix that
#   file = "/Users/ricard/data/mofaplus/hdf5/test.hdf5"
#   par =  h5read(file, "parameters")
#
#   tmp <- predict2(object, par, views=views, factors=factors)
#   pred_mean = tmp[[1]]
#   pred_var = tmp[[2]]
#
#   # imputed should be a list containing mean and variance for each
#   imputed_mean <- get_data(object, views=views, groups=groups, add_intercept = FALSE)
#   imputed_var <- get_data(object, views=views, groups=groups, add_intercept = FALSE)
#   for (m in views) {
#     for (g in groups) {
#       non_observed <- which(is.na(imputed_mean[[m]][[g]]), arr.ind = T)
#
#       imputed_mean[[m]][[g]][non_observed] <- pred_mean[[m]][[g]][non_observed]
#       imputed_var[[m]][[g]][non_observed] <- pred_var[[m]][[g]][non_observed]
#       imputed_var[[m]][[g]][-non_observed] <- 0 # TODO check that
#     }
#   }
#
#   # Save imputed data in the corresponding slot
#   object@imputed_data <- imputed
#   object@imputed_uncertainty <- imputed
#
#   return(object)
# }


impute.scatterplot <- function(object, view, features, factor, groups = "all") {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get view
  if (is.numeric(view)) view <- views_names(object)[view]
  stopifnot(view %in% views_names(object))
  
  # Get groups
  groups <- .check_and_get_groups(object, groups)
  
  # Get factor
  if (is.numeric(factor)) factor <- factors_names(object)[factor]
  stopifnot(factor %in% factors_names(object)) 
  Z.df <- get_factors(object, groups=groups, factors=factor, as.data.frame = TRUE)
  colnames(Z.df) <- c("sample", "factor", "factor.value", "group")
  
  # Fetch data
  imputed.df <- get_imputed_data2(object, views = view, groups = groups, features = "all", as.data.frame = TRUE)
  df <- merge(Z.df, imputed.df, by=c("sample","group"))
  df1 <- df[df$feature=="acc_H3K27ac_distal_E7.5_Mes_intersect12_1329",]
  # df1 <- df[df$feature=="Mesp1",]
  
  df2 <- reshape2::dcast(df1, formula=sample+view+group+feature+factor+factor.value~estimate, value.var="value")
  df2$stdev <- sqrt(df2$variance)
  
  df2$stdev[is.na(df2$stdev)] <- 0
  ggplot(df2, aes_string(x="factor.value", y="mean")) +
    facet_wrap(~feature) +
    geom_pointrange(aes(ymin=mean-stdev, ymax=mean+stdev), data=df2[df2$stdev==0,], size=0.1, alpha=0.5, color="grey70") +
    geom_pointrange(aes(ymin=mean-stdev, ymax=mean+stdev), data=df2[df2$stdev>0,], size=0.5) +
    theme_classic()
}





impute.heatmap <- function(object, view, features, factor, groups = "all") {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get view
  if (is.numeric(view)) view <- views_names(object)[view]
  stopifnot(view %in% views_names(object))
  
  # Get groups
  groups <- .check_and_get_groups(object, groups)
  
  # Get factor
  if (is.numeric(factor)) factor <- factors_names(object)[factor]
  stopifnot(factor %in% factors_names(object)) 
  Z.df <- get_factors(object, groups=groups, factors=factor, as.data.frame = TRUE)
  colnames(Z.df) <- c("sample", "factor", "factor.value", "group")
  
  # Fetch data
  df <- get_imputed_data2(object, views = view, groups = groups, features = "all", as.data.frame = TRUE)
  df <- merge(Z.df, df, by=c("sample","group"))
  # df <- df[df$feature=="feature0_view0" | df$feature=="feature0_view1",]
  
  # df2 <- reshape2::dcast(df, formula=sample+view+group+feature+factor+factor.value~estimate, value.var="value")
  # df2$stdev <- sqrt(df2$variance)
  
  p1 <- ggplot(df[df$estimate=="variance",], aes_string(x="sample", y="feature")) + 
    geom_tile(aes_string(fill = "value"), colour = "white") +
    facet_wrap(~estimate) +
    scale_fill_gradient(low = "white", high = "steelblue", na.value = 'grey') +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
  p1
  
  
}


filter.imputation <- function(object, cutoff = c("stringent", "medium", "lenient")) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  cutoff = match.arg(cutoff)
  
  # Define numeric cutoff
  predefined.cutoffs <- c(
    "stringent" = 2,
    "medium" = 5,
    "lenient" = 10
  )
  cutoff = predefined.cutoffs[cutoff]
  
  # calculate signal-to-noise ratio as sdev/mean and filter imputation results by signal-to-noise ratio
  sig2noise <- list()
  for (m in views_names(object)) {
    sig2noise[[m]] <- list()
    for (g in groups_names(object)) {
      sig2noise[[m]][[g]] <- sqrt(object@imputed_data[[m]][[g]]$var) / abs(object@imputed_data[[m]][[g]]$mean)
      object@imputed_data[[m]][[g]]$mean[ sig2noise[[m]][[g]]>cutoff ] <- NaN
      object@imputed_data[[m]][[g]]$variance[ sig2noise[[m]][[g]]>cutoff ] <- NaN
    }
  }
  
  return(object)
}