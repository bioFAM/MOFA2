
#' @title Calculate variance explained by the model
#' @name calculate_variance_explained
#' @description Method to calculate variance explained by the BioFAModel for each view and latent factor. \cr
#' As a measure of variance explained for gaussian data we adopt the coefficient of determination (R2). \cr
#' For non-gaussian views the calculations are based on the normally-distributed pseudo-data 
#' (for more information on the non-gaussian model see Supplementary Methods of the MOFA paper or Seeger & Bouchard, 2012).
#' @param object a \code{\link{BioFAModel}} object.
#' @param views character vector with the view names, or numeric vector with view indexes. Default is 'all'
#' @param factors character vector with the factor names, or numeric vector with the factor indexes. Default is 'all'
#' @param include_intercept include the intercept factor for calculation of variance explained (only used when an intercept was learned)
#' @param groups character vector with the group names, or numeric vector with group indexes. Default is 'all'
#' @details This function takes a trained BioFAModel as input and calculates for each view the coefficient of determination (R2),
#' i.e. the proportion of variance in the data explained by the BioFAM factor(s) (both jointly and for each individual factor). 
#' In case of non-Gaussian data the variance explained on the Gaussian pseudo-data is calculated. 
#' @return a list with matrices with the amount of variation explained per factor and view, and optionally total variance explained per view and variance explained by each feature alone
#' @export
calculate_variance_explained <- function(object, views = "all", groups = "all", factors = "all", include_intercept = TRUE, ...) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  
  # check whether the intercept was learned
  if(!object@model_options$learn_intercept & include_intercept) {
    include_intercept <- FALSE
    # warning("No intercept was learned in BioFAM.\n Intercept is not included in the model prediction.")
  }
  
  # Define views
  if (paste0(views, sep="", collapse="") == "all") { 
    views <- views_names(object) 
  } else {
    stopifnot(all(views %in% views_names(object)))  
  }
  M <- length(views)

  # Define groups
  if (paste0(groups, sep="", collapse="") == "all") { 
    groups <- groups_names(object) 
  } else {
    stopifnot(all(groups %in% groups_names(object)))  
  }
  H <- length(groups)

  # Define factors
  if (paste0(factors, collapse="") == "all") { 
    factors <- factors_names(object) 
  } else if (is.numeric(factors)) {
    if (include_intercept == T) {
      factors <- factors_names(object)[factors+1] 
    } else {
      factors <- factors_names(object)[factors]
    }
  } else { 
    stopifnot(all(factors %in% factors_names(object))) 
  }
  factors <- factors[factors!="intercept"]
  K <- length(factors)

  # Collect relevant expectations
  W <- get_weights(object, views, factors)
  Z <- get_factors(object, groups, factors)
  Y <- get_expectations(object, "Y")  # for non-Gaussian likelihoods the pseudodata is considered
  Y <- lapply(Y, function(m) lapply(m, t))
  
  # Calulcate feature-wise means as null model
  feature_mean <- lapply(views, function(m) {
    lapply(groups, function(h) {
      apply(Y[[m]][[h]], 2, mean, na.rm=T) 
    })
  })
  feature_mean <- .name_views_and_groups(feature_mean, views, groups)

  # Sweep out the feature-wise mean to calculate null model residuals
  res_null_model <- lapply(views, function(m) {
    lapply(groups, function(p) {
      sweep(Y[[m]][[p]], 2, feature_mean[[m]][[p]], "-")
    })
  })
  res_null_model <- .name_views_and_groups(res_null_model, views, groups)
  
  # replace masked values on Z by 0 (so that they do not contribute to predictions)
  for (group in groups) { Z[[group]][is.na(Z[[group]])] <- 0 }
    
  # Calculate predictions under the MOFA model using all (non-intercept) factors
  Ypred_m <- lapply(views, function(m) {
    lapply(groups, function(h) {
      Z[[h]] %*% t(W[[m]])
    })
  })
  Ypred_m <- .name_views_and_groups(Ypred_m, views, groups)

  for (view in views) { names(res_null_model[[view]]) <- groups }

  # Calculate predictions under the MOFA model using each (non-intercept) factors on its own
  Ypred_mk <- lapply(views, function(m) {
    lapply(groups, function(h) {
      ltmp <- lapply(factors, function(k) Z[[h]][,k] %*% t(W[[m]][,k]) )
      names(ltmp) <- factors
      ltmp
    })
  })
  Ypred_mk <- .name_views_and_groups(Ypred_mk, views, groups)
  
  # If an intercept is included, regress out the intercept from the data
  if (include_intercept) {
      intercept <- get_weights(object,views,"intercept")
      Y <- lapply(views, function(m) lapply(groups, function(h) sweep(Y[[m]][[h]], 2, intercept[[m]], "-")))
      Y <- .name_views_and_groups(Y, views, groups)
  }

  # Calculate coefficient of determination per view
  fvar_m <- lapply(groups, function(h) lapply(views, function(m) 1 - sum((Y[[m]][[h]]-Ypred_m[[m]][[h]])**2, na.rm=T) / sum(res_null_model[[m]][[h]]**2, na.rm=T)))
  fvar_m <- .name_views_and_groups(fvar_m, groups, views)

  # Calculate coefficient of determination per factor and view
  tmp <- lapply(groups, function(h) {
    sapply(views, function(m) {
      sapply(factors, function(k) {
        1 - sum((Y[[m]][[h]]-Ypred_mk[[m]][[h]][[k]])**2, na.rm=T) / sum(res_null_model[[m]][[h]]**2, na.rm=T) 
      })
    })
  })
  fvar_mk <- lapply(tmp, function(e) matrix(e, ncol=length(views), nrow=length(factors)))
  names(fvar_mk) <- groups
  for (h in groups) { colnames(fvar_mk[[h]]) <- views; rownames(fvar_mk[[h]]) <- factors }

  # Store results
  r2_list <- list(r2_total = fvar_m, r2_per_factor = fvar_mk)
  
  return(r2_list)
 
}


#' @title Plot variance explained by the model
#' @name plot_variance_explained
#' @description Method to plot variance explained (R-squared) by the MOFA model for each view, each group, and latent factor. \cr
#' As a measure of variance explained for gaussian data we adopt the coefficient of determination (R2). \cr
#' For details on the computation see the help of the \code{\link{calculate_variance_explained}} function
#' @param object a \code{\link{MOFAmodel}} object.
#' @param cluster logical indicating whether to do hierarchical clustering on the plot
#' @param ... extra arguments to be passed to \code{\link{calculate_variance_explained}}
#' @return ggplot object
#' @import pheatmap ggplot2 reshape2
#' @importFrom cowplot plot_grid
#' @export
plot_variance_explained <- function(object, views = "all", groups = "all", cluster = T, ...) {

  views  <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)
  
  # Calculate Variance Explained
  R2_list <- calculate_variance_explained(object, ...)
  fvar_m  <- lapply(R2_list$R2Total[groups], function(e) e[views])
  fvar_mk <- lapply(R2_list$R2PerFactor[groups], function(e) e[,views])

  ## Plot variance explained by factor ##
  
  # convert matrix to data frame for ggplot2
  fvar_mk_df <- reshape2::melt(
    lapply(fvar_mk, function(gr) 
      reshape2::melt(as.matrix(gr), varnames = c("factor", "view"))
    ), id.vars=c("factor", "view", "value")
  )
  colnames(fvar_mk_df)[ncol(fvar_mk_df)] <- "group"
  fvar_mk_df$factor <- factor(fvar_mk_df$factor)
  ## Plot variance explained per view ##
  
  # Create data.frame for ggplot
  fvar_m_df <- melt(lapply(fvar_m, function(e) lapply(e, function(x) x)), 
                    varnames=c("view", "group"),
                    value.name="R2")
  colnames(fvar_m_df)[(ncol(fvar_m_df)-1):ncol(fvar_m_df)] <- c("view", "group")

  hms   <- list()
  bplts <- list()

  for (gr in unique(fvar_mk_df$group)) {

    # Grid plot with the variance explained per factor and view
    hm <- ggplot(fvar_mk_df[fvar_mk_df$group == gr,], aes(view, factor)) + 
      geom_tile(aes(fill=value), color="black") +
      guides(fill=guide_colorbar("R2")) +
      scale_fill_gradientn(colors=c("gray97","darkblue"), guide="colorbar") +
      ylab("Latent factor") +
      theme(
        # plot.margin = margin(5,5,5,5),
        plot.title = element_text(size=17, hjust=0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=11, angle=60, hjust=1, vjust=1, color="black"),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.y = element_text(size=15),
        axis.line = element_blank(),
        axis.ticks =  element_blank(),
        panel.background = element_blank()
      )
    hm <- hm + ggtitle(paste0("Variance explained per factor\nin group ", gr))  + 
      guides(fill=guide_colorbar("R2"))
    hms[[gr]] <- hm
    
    # Barplot with variance explained per view
    bplt <- ggplot(fvar_m_df[fvar_m_df$group == gr,], aes(x=view, y=R2)) + 
      ggtitle(paste0("Total variance explained per view\nin group ", gr)) +
      geom_bar(stat="identity", fill="deepskyblue4", width=0.9) +
      xlab("") + ylab("R2") +
      scale_y_continuous(expand=c(0.01,0.01)) +
      theme(
        plot.margin = unit(c(1,2.4,0,0), "cm"),
        panel.background = element_blank(),
        plot.title = element_text(size=17, hjust=0.5),
        axis.ticks.x = element_blank(),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12, color="black"),
        axis.title.y = element_text(size=13, color="black"),
        axis.line = element_line(size=rel(1.0), color="black")
      )
    bplts[[gr]] <- bplt
  }
  
  # Join the two plots
  p <- plot_grid(plotlist = c(bplts, hms), align="v", nrow=2, ncol=length(unique(fvar_mk_df$group)), rel_heights=c(1/3,2/3), axis="l")
  
  return(p)
}


#' @title Plot variance explained by the model
#' @name plot_variance_explained_in_views
#' @description Method to plot variance explained (R-squared) by the MOFA model for each view and latent factor. \cr
#' As a measure of variance explained for gaussian data we adopt the coefficient of determination (R2). \cr
#' For details on the computation see the help of the \code{\link{calculate_variance_explained}} function
#' @param object a \code{\link{MOFAmodel}} object.
#' @param cluster logical indicating whether to do hierarchical clustering on the plot
#' @param ... extra arguments to be passed to \code{\link{calculate_variance_explained}}
#' @return ggplot object
#' @import pheatmap ggplot2 reshape2
#' @importFrom cowplot plot_grid
#' @export
plot_variance_explained_in_views <- function(object, groups = "all", cluster = T, ...) {

  groups <- .check_and_get_groups(object, groups)
  
  # Calculate Variance Explained
  R2_list <- calculate_variance_explained(object, ...)
  fvar_m  <- Reduce(function(a, b) Map('+', a, b), R2_list$R2Total[groups])
  fvar_mk <- Reduce('+', R2_list$R2PerFactor[groups])

  ## Plot variance explained by factor ##
  
  # convert matrix to data frame for ggplot2  
  fvar_mk_df <- reshape2::melt(fvar_mk, varnames = c("factor", "view"))
  fvar_mk_df$factor <- factor(fvar_mk_df$factor)
  
  # If multiple views, sort factors according to hierarchical clustering
  if (cluster==TRUE & ncol(fvar_mk)>1) {
    hc <- hclust(dist(t(fvar_mk)))
    fvar_mk_df$view <- factor(fvar_mk_df$view, levels = colnames(fvar_mk)[hc$order])
  }
  
  # Grid plot with the variance explained per factor and view
  hm <- ggplot(fvar_mk_df, aes(view,factor)) + 
    geom_tile(aes(fill=value), color="black") +
    guides(fill=guide_colorbar("R2")) +
    scale_fill_gradientn(colors=c("gray97","darkblue"), guide="colorbar") +
    ylab("Latent factor") +
    theme(
      # plot.margin = margin(5,5,5,5),
      plot.title = element_text(size=17, hjust=0.5),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size=11, angle=60, hjust=1, vjust=1, color="black"),
      axis.text.y = element_text(size=12, color="black"),
      axis.title.y = element_text(size=15),
      axis.line = element_blank(),
      axis.ticks =  element_blank(),
      panel.background = element_blank()
    )
  hm <- hm + ggtitle("Variance explained per factor")  + 
    guides(fill=guide_colorbar("R2"))
  
  ## Plot variance explained per view ##
  
  # Create data.frame for ggplot
  fvar_m_df <- data.frame(view=names(fvar_m), R2=unlist(fvar_m))
  
  # If multiple views, sort factors according to hierarchical clustering
  if (cluster==TRUE & ncol(fvar_mk)>1) {
    fvar_m_df$view <- factor(fvar_m_df$view, levels = colnames(fvar_mk)[hc$order])
  }
  
  # Barplot with variance explained per view
  bplt <- ggplot( fvar_m_df, aes(x=view, y=R2)) + 
    ggtitle("Total variance explained per view") +
    geom_bar(stat="identity", fill="deepskyblue4", width=0.9) +
    xlab("") + ylab("R2") +
    scale_y_continuous(expand=c(0.01,0.01)) +
    theme(
      plot.margin = unit(c(1,2.4,0,0), "cm"),
      panel.background = element_blank(),
      plot.title = element_text(size=17, hjust=0.5),
      axis.ticks.x = element_blank(),
      # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)
      axis.text.x = element_blank(),
      axis.text.y = element_text(size=12, color="black"),
      axis.title.y = element_text(size=13, color="black"),
      axis.line = element_line(size=rel(1.0), color="black")
    )
  
  # Join the two plots
  p <- plot_grid(bplt, hm, align="v", nrow=2, rel_heights=c(1/3,2/3), axis="l")
  
  return(p)
}


#' @title Plot variance explained by the model in groups
#' @name plot_variance_explained_in_groups
#' @description Method to plot variance explained (R-squared) by the MOFA model for each group and latent factor. \cr
#' As a measure of variance explained for gaussian data we adopt the coefficient of determination (R2). \cr
#' For details on the computation see the help of the \code{\link{calculate_variance_explained}} function
#' @param object a \code{\link{MOFAmodel}} object.
#' @param cluster logical indicating whether to do hierarchical clustering on the plot
#' @param ... extra arguments to be passed to \code{\link{calculate_variance_explained}}
#' @return ggplot object
#' @import pheatmap ggplot2 reshape2
#' @importFrom cowplot plot_grid
#' @export
plot_variance_explained_in_groups <- function(object, views = "all", cluster = T, ...) {

  views <- .check_and_get_views(object, views)
  
  # Calculate Variance Explained
  R2_list <- calculate_variance_explained(object, ...)
  fvar_m  <- lapply(R2_list$R2Total, function(p) Reduce('+', p[views]))
  fvar_mk <- sapply(R2_list$R2PerFactor, function(p) rowSums(as.matrix(p[,views])))

  ## Plot variance explained by factor ##
  
  # convert matrix to data frame for ggplot2  
  fvar_mk_df <- reshape2::melt(fvar_mk, varnames = c("factor", "group"))
  fvar_mk_df$factor <- factor(fvar_mk_df$factor)
  
  # If multiple views, sort factors according to hierarchical clustering
  if (cluster==TRUE & ncol(fvar_mk)>1) {
    hc <- hclust(dist(t(fvar_mk)))
    fvar_mk_df$group <- factor(fvar_mk_df$group, levels = colnames(fvar_mk)[hc$order])
  }
  
  # Grid plot with the variance explained per factor and group
  hm <- ggplot(fvar_mk_df, aes(group, factor)) + 
    geom_tile(aes(fill=value), color="black") +
    guides(fill=guide_colorbar("R2")) +
    scale_fill_gradientn(colors=c("gray97","darkblue"), guide="colorbar") +
    ylab("Latent factor") +
    theme(
      # plot.margin = margin(5,5,5,5),
      plot.title = element_text(size=17, hjust=0.5),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size=11, angle=60, hjust=1, vjust=1, color="black"),
      axis.text.y = element_text(size=12, color="black"),
      axis.title.y = element_text(size=15),
      axis.line = element_blank(),
      axis.ticks =  element_blank(),
      panel.background = element_blank()
    )
  hm <- hm + ggtitle("Variance explained per factor")  + 
    guides(fill=guide_colorbar("R2"))
  
  ## Plot variance explained per group ##
  
  # Create data.frame for ggplot
  fvar_m_df <- data.frame(group=names(fvar_m), R2=unlist(fvar_m))
  
  # If multiple groups, sort factors according to hierarchical clustering
  if (cluster==TRUE & ncol(fvar_mk)>1) {
    fvar_m_df$group <- factor(fvar_m_df$group, levels = colnames(fvar_mk)[hc$order])
  }
  
  # Barplot with variance explained per group
  bplt <- ggplot( fvar_m_df, aes(x=group, y=R2)) + 
    ggtitle("Total variance explained per group") +
    geom_bar(stat="identity", fill="deepskyblue4", width=0.9) +
    xlab("") + ylab("R2") +
    scale_y_continuous(expand=c(0.01,0.01)) +
    theme(
      plot.margin = unit(c(1,2.4,0,0), "cm"),
      panel.background = element_blank(),
      plot.title = element_text(size=17, hjust=0.5),
      axis.ticks.x = element_blank(),
      # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)
      axis.text.x = element_blank(),
      axis.text.y = element_text(size=12, color="black"),
      axis.title.y = element_text(size=13, color="black"),
      axis.line = element_line(size=rel(1.0), color="black")
    )
  
  # Join the two plots
  p <- plot_grid(bplt, hm, align="v", nrow=2, rel_heights=c(1/3,2/3), axis="l")
  
  return(p)
}

