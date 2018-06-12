
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
#' @param groupwise boolean indicating whether the variance explained is calculated using a groupwise intercept for the features or an overall
#' @details This function takes a trained BioFAModel as input and calculates for each view the coefficient of determination (R2),
#' i.e. the proportion of variance in the data explained by the BioFAM factor(s) (both jointly and for each individual factor). 
#' In case of non-Gaussian data the variance explained on the Gaussian pseudo-data is calculated. 
#' @return a list with matrices with the amount of variation explained per factor and view, and optionally total variance explained per view and variance explained by each feature alone
#' @export
calculate_variance_explained <- function(object, views = "all", groups = "all", factors = "all", 
                                         include_intercept = TRUE, only = NULL, flatten = FALSE, groupwise = FALSE, ...) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  
  # check whether the intercept was learned
  if(!object@model_options$learn_intercept & include_intercept) {
    include_intercept <- FALSE
    # warning("No intercept was learned in BioFAM.\n Intercept is not included in the model prediction.")
  }
  
  # Define views and groups
  views  <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)
  
  if (is.null(only) & flatten) {
    if (length(groups) == 1) {
      only <- "views"
    } else if (length(views) == 1) {
      only <- "groups"
    }
  }
  
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
  
  if (is.null(only)) {
    
    # Calulcate feature-wise means as null model
    feature_mean <- lapply(views, function(m) {
      lapply(groups, function(h) {
       if (groupwise) {
         #apply(Y[[m]][[h]], 2, mean, na.rm=T) 
         apply(Reduce(rbind, Y[[m]]), 2, mean, na.rm=T)
       } 
       else {
         #apply(Reduce(rbind,Y[[m]]), 2, mean, na.rm=T) 
         apply(Y[[m]][[h]], 2, mean, na.rm=T)
       }
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
    
    # Calculate predictions under the BioFAModel using all (non-intercept) factors
    Ypred_m <- lapply(views, function(m) {
      lapply(groups, function(h) {
        Z[[h]] %*% t(W[[m]])
      })
    })
    Ypred_m <- .name_views_and_groups(Ypred_m, views, groups)
    
    for (view in views) { names(res_null_model[[view]]) <- groups }
    
    # Calculate predictions under the BioFAModel using each (non-intercept) factor on its own
    Ypred_mk <- lapply(views, function(m) {
      lapply(groups, function(p) {
        ltmp <- lapply(factors, function(k) Z[[p]][,k] %*% t(W[[m]][,k]) )
        names(ltmp) <- factors
        ltmp
      })
    })
    Ypred_mk <- .name_views_and_groups(Ypred_mk, views, groups)
    
    # If an intercept is included, regress out the intercept from the data
    if (include_intercept) {
      intercept <- get_weights(object, views, "intercept")
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
    
  } else if ((only == "views") | (only == "groups")) {
    
    if (only == "groups") {
      
      W <- Reduce(rbind, W)
      Y <- lapply(groups, function(p) {
        views_list <- lapply(views, function(m) {
          Y[[m]][[p]]
        })
        Reduce(cbind, views_list)
      })
      lnames <- groups
      
      # Replace masked values on Z by 0 (so that they do not contribute to predictions)
      for (group in groups) { Z[[group]][is.na(Z[[group]])] <- 0 }
      
      # Calculate predictions under the BioFAModel using all (non-intercept) factors
      Ypred_s <- lapply(groups, function(p) {
        Z[[p]] %*% t(W)
      })
      
      # Calculate predictions under the BioFAModel using each (non-intercept) factors on its own
      # while combining different views
      Ypred_sk <- lapply(lnames, function(p) {
        tmp <- lapply(factors, function(k) Z[[p]][,k] %*% t(W[,k]) )
        names(tmp) <- factors
        tmp
      })
      
    } else {
      
      Z <- Reduce(rbind, Z)
      Y <- lapply(views, function(m) {
        Reduce(rbind, Y[[m]])
      })
      lnames <- views
      
      # Replace masked values on Z by 0 (so that they do not contribute to predictions)
      Z[is.na(Z)] <- 0
      
      # Calculate predictions under the BioFAModel using all (non-intercept) factors
      Ypred_s <- lapply(views, function(m) {
        Z %*% t(W[[m]])
      })
      
      
      # Calculate predictions under the BioFAModel using each (non-intercept) factors on its own
      # while combining different views
      Ypred_sk <- lapply(views, function(m) {
        tmp <- lapply(factors, function(k) Z[,k] %*% t(W[[m]][,k]) )
        names(tmp) <- factors
        tmp
      })
      
      
    }
    
    names(Y) <- lnames
    names(Ypred_s) <- lnames
    names(Ypred_sk) <- lnames
    
    # Calulcate feature-wise means as null model
    feature_mean <- lapply(lnames, function(s) {
      apply(Y[[s]], 2, mean, na.rm=T) 
    })
    names(feature_mean) <- lnames
    
    # Sweep out the feature-wise mean to calculate null model residuals
    res_null_model <- lapply(lnames, function(s) {
      sweep(Y[[s]], 2, feature_mean[[s]], "-")
    })
    names(res_null_model) <- lnames
    
    # If an intercept is included, regress out the intercept from the data
    if (include_intercept) {
      if (only == "groups") {
        intercept <- get_weights(object, views, "intercept")
        intercept <- Reduce(rbind, intercept)
        Y <- lapply(groups, function(p) sweep(Y[[p]], 2, intercept, "-"))
        names(Y) <- groups
      } else {
        intercept <- get_weights(object, views, "intercept")
        Y <- lapply(views, function(m) sweep(Y[[m]], 2, intercept[[m]], "-"))
        names(Y) <- views
      }
    }
    
    # Calculate coefficient of determination per group
    fvar_s <- lapply(lnames, function(s)
      1 - sum((Y[[s]] - Ypred_s[[s]])**2, na.rm=T) / sum(res_null_model[[s]]**2, na.rm=T)
    )
    names(fvar_s) <- lnames
    
    # Calculate coefficient of determination per factor and group
    fvar_sk <- lapply(lnames, function(s) {
      sapply(factors, function(k) {
        1 - sum((Y[[s]] - Ypred_sk[[s]][[k]])**2, na.rm=T) / sum(res_null_model[[s]]**2, na.rm=T) 
      })
    })
    names(fvar_sk) <- lnames
    
    # Store results
    r2_list <- list(r2_total = sapply(fvar_s, function(e) e), r2_per_factor = sapply(fvar_sk, function(e) e))
    
    return(r2_list)
    
  } else {
    stop("Please provide `only` as `views` or `groups` or leave it empty")
  }
  
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
plot_variance_explained <- function(object, views = "all", groups = "all", only = NULL, split_by = NULL, cluster = TRUE, ...) {

  views  <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)

  if (is.null(only)) {
    if (length(groups) == 1) {
      only <- "views"
    } else if (length(views) == 1) {
      only <- "groups"
    }
  }

  # Sanity checks
  if (is.null(split_by)) {
    if (is.null(only)) {
      split_by <- "group"  # default
    }
  } else {
    if (!is.null(only)) {
      print(paste0("Warning: `split_by` can't be set together with `only` or when there's only 1 group or 1 view. Not using it for plotting."))
    }
    if (split_by != "view" & split_by != "group") {
      print(paste0("Warning: `split_by` can be either `view` or `group`. Setting it to `group` for plotting..."))
      split_by <- "group"  # default
    }
  }

  if (!is.null(only)) {
    if (only == "views") { sub <- "view"; scnd <- "group" } else { sub <- "group"; scnd <- "view" }
  } else {
    # There are >1 views and >1 groups
    if (split_by == "group") { sub <- "view"; scnd <- "group" } else { sub <- "group"; scnd <- "view" }
  }

  # Calculate Variance Explained
  r2_list <- calculate_variance_explained(object, views = views, groups = groups, only = only, flatten = TRUE, ...)
  
  if (is.null(only)) {

    fvar_m  <- lapply(r2_list$r2_total[groups], function(e) e[views])
    fvar_mk <- lapply(r2_list$r2_per_factor[groups], function(e) e[,views])

    ## Plot variance explained by factor ##
    
    # Convert matrix to data frame for ggplot2
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


  } else if (only == "views" | only == "groups") {

    fvar_m  <- r2_list$r2_total
    fvar_mk <- r2_list$r2_per_factor

    # Convert matrix to data frame for ggplot2  
    fvar_mk_df <- reshape2::melt(fvar_mk, varnames = c("factor", sub))
    fvar_mk_df$factor <- factor(fvar_mk_df$factor)
    
    # If multiple views/groups, sort factors according to hierarchical clustering
    if (cluster==TRUE & ncol(fvar_mk)>1) {
      hc <- hclust(dist(t(fvar_mk)))
      fvar_mk_df[[sub]] <- factor(fvar_mk_df[[sub]], levels = colnames(fvar_mk)[hc$order])
    }

    ## Plot variance explained per view ##
    
    # Create data.frame for ggplot
    fvar_m_df <- data.frame(names(fvar_m), unlist(fvar_m))
    colnames(fvar_m_df) <- c(sub, "R2")
    
    # If multiple views, sort factors according to hierarchical clustering
    if (cluster==TRUE & ncol(fvar_mk)>1) {
      fvar_mk_df[[sub]] <- factor(fvar_mk_df[[sub]], levels = colnames(fvar_mk)[hc$order])
    }

    fvar_mk_df[[scnd]] <- "all"
    fvar_m_df[[scnd]]  <- "all"


  } else {
    stop("Please provide `only` as views or groups or leave it empty")
  }
  
  hms   <- list()
  bplts <- list()
  min_lim_hm <- min(fvar_mk_df$value)
  max_lim_hm <- max(fvar_mk_df$value)
  min_lim_bplt <- min(0, fvar_m_df$R2)
  max_lim_bplt <- max(fvar_m_df$R2)

  for (s in unique(fvar_mk_df[[scnd]])) {

    mk_title <- paste0("Variance explained per factor")
    if (is.null(only)) mk_title <- paste0(mk_title, "\nin ", scnd ," ", s)

    # Grid plot with the variance explained per factor and view
    hm <- ggplot(fvar_mk_df[fvar_mk_df[[scnd]] == s,], aes_string(sub, "factor")) + 
      geom_tile(aes(fill=value), color="black") +
      guides(fill=guide_colorbar("R2")) +
      ylab("Latent factor") + 
      scale_fill_gradientn(colors=c("gray97","darkblue"), guide="colorbar", limits=c(min_lim_hm,max_lim_hm)) +
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
    hm <- hm + ggtitle(mk_title)  + 
      guides(fill=guide_colorbar("R2"))
    hms[[s]] <- hm

    m_title <- paste0("Total variance explained per ", sub)
    if (is.null(only)) m_title <- paste0(m_title, "\nin ", scnd, " ", s)
    
    # Barplot with variance explained per view
    bplt <- ggplot(fvar_m_df[fvar_m_df[[scnd]] == s,], aes_string(x=sub, y="R2")) + 
      ggtitle(m_title) +
      geom_bar(stat="identity", fill="deepskyblue4", width=0.9) +
      xlab("") + ylab("R2") +
      scale_y_continuous(limits=c(min_lim_bplt, max_lim_bplt), expand=c(0.01, 0.01)) +
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
    bplts[[s]] <- bplt
  }
  
  # Join the two plots
  p <- plot_grid(plotlist = c(bplts, hms), align="v", nrow=2, ncol=length(unique(fvar_mk_df[[scnd]])), rel_heights=c(1/3,2/3), axis="l")
  
  return(p)
}


