
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


impute.plot <- function(object, view, features, factor, groups = "all") {
  
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
  df <- df[df$feature=="feature0_view0",]
  

  df2 <- reshape2::dcast(df, formula=sample+view+group+feature+factor+factor.value~estimate, value.var="value")
  df2$stdev <- sqrt(df2$variance)
  
  ggplot(df2, aes_string(x="factor.value", y="mean")) +
    geom_pointrange(aes(ymin=mean-stdev, ymax=mean+stdev)) +
    theme_classic()
