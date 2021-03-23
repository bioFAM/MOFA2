#######################################################
## Functions to prepare data for input to MOFA ##
#######################################################

#' @title Process count data for use with MOFA
#' @name preprocess_data
#' @description Function to preprocess data for input to MOFA. 
#' Note that this step is optional and users can provide their already preprocessed
#' data directly to `create_mofa_object` (recommended). 
#' @param data A list of matrices, where each entry corresponds to one view. 
#' Samples are stored in columns and features in rows. 
#' Missing values must be filled in prior to creating the MOFA object 
#' (see for example the CLL tutorial).
#' @param groups group information, only relevant when using the multi-group framework. 
#' @param n_top number of features to keep for each view based on highest variance,
#' needs to be a vector of same length as data, by default all features are used
#' @param type Can be on of 'gaussian', 'count' or 'binary'. 
#' If not provided, is guessed from the data.
#' @details Selects the most informative features based on variance across samples (within groups) and 
#' transforms count data to Pearson residuals.
#' @return Returns an data in the same format as provided.
#' @export
#' @examples
#' data <- make_example_data(likelihood = "poisson")$data
#' pp_data <- preprocess_data(data)

preprocess_data <- function(data, groups = NULL, n_top = NULL, type = NULL) {
  
  if(!is.null(groups)){
    message("Note that using the groups option is a advanced option.
            Make sure you use the same group specification in create_mofa.
            This will select variables based on variation left after taking the group effect into account
            and will remove mean differences between groups and center per group.
            Set to NULL if you are interested in these.")
  }
  if(is.null(n_top)){
    n_top <- sapply(data, nrow)
  }
  if(length(n_top) == 1){
    n_top = rep(n_top, length(data))
  }
  if(is.null(type)){
    type <- .guess_type(data)
  }
  # Filtering the according to top N and preprocess for count data
  data_new <- lapply(seq_len(length(data)), function(m){
    
    dd <- data[[m]]
    
    if (type[m] %in% c("gaussian", "binary")) {
      if( nrow(dd) < n_top[m]){
        return(dd) # do nothing
      } else { 
        # filter top n_top features by variance
        if(is.null(groups) | length(unique(groups)) == 1){
          vars <- apply(dd,1,var)
          dd_new <- dd[order(vars, decreasing = TRUE)[seq_len(n_top[m])],]
        } else {
          # filter top n_top features by variance *within group*
          fit <- lm(t(dd) ~ groups)
          vars <- apply(fit$residuals,1,var)
          dd_new <- dd[order(vars, decreasing = TRUE)[seq_len(n_top[m])],]
        }
        return(dd_new)
      }
        
      } else if(type[m] == "counts") {
        # Pearson residuals as suggested by Townes et al, Gen Biol 2019
        sz <- Matrix::colSums(dd)
        if (!is.null(groups)) {
          dd_new <- matrix(NA, ncol = ncol(dd), nrow = nrow(dd))
          for (g in unique(groups)){
          phat <- Matrix::rowSums(dd[,groups == g])/sum(sz[groups == g])
          mhat <- outer(phat, sz[groups == g])
          dd_new[,groups == g] <- (dd[,groups == g] - mhat) / sqrt(mhat * (1-phat)) # binomial
          # dd_new[,groups == g] <- (dd[,groups == g] - mhat) / sqrt(mhat + 100 * mhat**2)  # neg. binomial with fixed theta
          }
        } else {
            phat <- Matrix::rowSums(dd)/sum(sz)
            mhat <- outer(phat, sz)
            dd_new <- (dd - mhat) / sqrt(mhat * (1-phat)) # binomial
            # dd_new <- (dd - mhat) / sqrt(mhat + 100 * mhat**2) # neg. binomial with fixed theta: (9) in Lause et al
        }
        if( nrow(dd) < n_top[m]){
          return(dd_new)
        } else{
          return(dd_new[order(apply(dd_new,1,var), decreasing = TRUE)[seq_len(n_top[m])],])
        }
      }
  })

  return(data_new)
}



.guess_type <- function(data){
  {
    type <- rep(x = "gaussian", times = length(data))
    
    for (m in seq_len(length(data))) {
      dd <- Reduce(cbind, data[[m]])
      if (length(unique(dd[!is.na(dd)])) == 2) {
        type[m] <- "binary"
      }
      else if (all(dd[!is.na(dd)]%%1 == 0)) {
        type[m] <- "counts"
      }
    }
    return(type)
  }
}
