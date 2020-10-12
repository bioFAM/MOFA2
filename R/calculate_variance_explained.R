#' @title Calculate variance explained by the model
#' @description  This function takes a trained MOFA model as input and calculates the proportion of variance explained 
#' (i.e. the coefficient of determinations (R^2)) by the MOFA factors across the different views.
#' @name calculate_variance_explained
#' @param object a \code{\link{MOFA}} object.
#' @param views character vector with the view names, or numeric vector with view indexes. Default is 'all'
#' @param groups character vector with the group names, or numeric vector with group indexes. Default is 'all'
#' @param factors character vector with the factor names, or numeric vector with the factor indexes. Default is 'all'
#' @return a list with matrices with the amount of variation explained per factor and view.
#' @importFrom utils relist as.relistable
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Calculate variance explained (R2)
#' r2 <- calculate_variance_explained(model)
#' 
#' # Plot variance explained values (view as x-axis, and factor as y-axis)
#' plot_variance_explained(model, x="view", y="factor")
#' 
#' # Plot variance explained values (view as x-axis, and group as y-axis)
#' plot_variance_explained(model, x="view", y="group")
#' 
#' # Plot variance explained values for factors 1 to 3
#' plot_variance_explained(model, x="view", y="group", factors=1:3)
#' 
#' # Scale R2 values
#' plot_variance_explained(model, max_r2 = 0.25)
calculate_variance_explained <- function(object, views = "all", groups = "all", factors = "all") {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (any(object@model_options$likelihoods!="gaussian"))
    stop("Not possible to recompute the variance explained estimates when using non-gaussian likelihoods.")
  if (any(object@model_options$likelihoods!="gaussian"))
    if (isFALSE(object@data_options$loaded)) stop("Data is not loaded, cannot compute variance explained.")
  
  # Define factors, views and groups
  views  <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)
  factors <- .check_and_get_factors(object, factors)
  K <- length(factors)
  
  # Collect relevant expectations
  W <- get_weights(object, views=views, factors=factors)
  Z <- get_factors(object, groups=groups, factors=factors)
  Y <- lapply(get_data(object, add_intercept = FALSE)[views], function(view) view[groups])
  Y <- lapply(Y, function(x) lapply(x,t))
  
  # Replace masked values on Z by 0 (so that they do not contribute to predictions)
  for (g in groups) {
    Z[[g]][is.na(Z[[g]])] <- 0
  }
  
  # Calculate coefficient of determination per group and view
  r2_m <- tryCatch({
    lapply(groups, function(g) sapply(views, function(m) {
      a <- sum((as.matrix(Y[[m]][[g]]) - tcrossprod(Z[[g]], W[[m]]))**2, na.rm = TRUE)
      b <- sum(Y[[m]][[g]]**2, na.rm = TRUE)
      return(1 - a/b)
    })
    )}, error = function(err) {
      stop(paste0("Calculating explained variance doesn't work with the current version of DelayedArray.\n",
                  "  Do not sort factors if you're trying to load the model (sort_factors = FALSE),\n",
                  "  or load the full dataset into memory (on_disk = FALSE)."))
      return(err)
    })
  r2_m <- .name_views_and_groups(r2_m, groups, views)
  
  # Lower bound is zero
  r2_m = lapply(r2_m, function(x){
    x[x < 0] = 0
    return(x)
  })
  
  # Calculate coefficient of determination per group, factor and view
  r2_mk <- lapply(groups, function(g) {
    tmp <- sapply(views, function(m) { sapply(factors, function(k) {
      a <- sum((as.matrix(Y[[m]][[g]]) - tcrossprod(Z[[g]][,k], W[[m]][,k]))**2, na.rm = TRUE)
      b <- sum(Y[[m]][[g]]**2, na.rm = TRUE)
      return(1 - a/b)
    })
    })
    tmp <- matrix(tmp, ncol = length(views), nrow = length(factors))
    colnames(tmp) <- views
    rownames(tmp) <- factors
    return(tmp)
  })
  names(r2_mk) <- groups
  
  # Lower bound is 0
  r2_mk = lapply(r2_mk, function(x){
    x[x < 0] = 0
    return(x)
  })
  
  # Transform from fraction to percentage
  r2_mk = utils::relist(unlist(utils::as.relistable(r2_mk)) * 100 ) 
  r2_m = utils::relist(unlist(utils::as.relistable(r2_m)) * 100 )
  
  # Store results
  r2_list <- list(r2_total = r2_m, r2_per_factor = r2_mk)
  
  return(r2_list)
}









#' @title Plot variance explained by the model
#' @description plots the variance explained by the MOFA factors across different views and groups, as specified by the user.
#' Consider using cowplot::plot_grid(plotlist = ...) to combine the multiple plots that this function generates.
#' @name plot_variance_explained
#' @param object a \code{\link{MOFA}} object
#' @param x character specifying the dimension for the x-axis ("view", "factor", or "group").
#' @param y character specifying the dimension for the y-axis ("view", "factor", or "group").
#' @param split_by character specifying the dimension to be faceted ("view", "factor", or "group").
#' @param factors character vector with a factor name(s), or numeric vector with the index(es) of the factor(s). Default is "all".
#' @param plot_total logical value to indicate if to plot the total variance explained (for the variable in the x-axis)
#' @param min_r2 minimum variance explained for the color scheme (default is 0).
#' @param max_r2 maximum variance explained for the color scheme.
#' @param legend logical indicating whether to add a legend to the plot  (default is TRUE).
#' @param use_cache logical indicating whether to use cache (default is TRUE)
#' @param ... extra arguments to be passed to \code{\link{calculate_variance_explained}}
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom stats as.formula
#' @importFrom reshape2 melt
#' @return A list of \code{\link{ggplot}} objects (if \code{plot_total} is TRUE) or a single \code{\link{ggplot}} object
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Calculate variance explained (R2)
#' r2 <- calculate_variance_explained(model)
#' 
#' # Plot variance explained values (view as x-axis, and factor as y-axis)
#' plot_variance_explained(model, x="view", y="factor")
#' 
#' # Plot variance explained values (view as x-axis, and group as y-axis)
#' plot_variance_explained(model, x="view", y="group")
#' 
#' # Plot variance explained values for factors 1 to 3
#' plot_variance_explained(model, x="view", y="group", factors=1:3)
#' 
#' # Scale R2 values
#' plot_variance_explained(model, max_r2=0.25)
plot_variance_explained <- function(object, x = "view", y = "factor", split_by = NA, plot_total = FALSE, 
                                    factors = "all", min_r2 = 0, max_r2 = NULL, legend = TRUE, use_cache = TRUE, ...) {
  
  # Sanity checks 
  if (length(unique(c(x, y, split_by))) != 3) { 
    stop(paste0("Please ensure x, y, and split_by arguments are different.\n",
                "  Possible values are `view`, `group`, and `factor`."))
  }
  
  # Automatically fill split_by in
  if (is.na(split_by)) split_by <- setdiff(c("view", "factor", "group"), c(x, y, split_by))
  
  # Calculate variance explained
  if ((use_cache) & .hasSlot(object, "cache") && ("variance_explained" %in% names(object@cache))) {
    r2_list <- object@cache$variance_explained
  } else {
    r2_list <- calculate_variance_explained(object, factors = factors, ...)
  }
  
  r2_mk <- r2_list$r2_per_factor
  
  # convert matrix to long data frame for ggplot2
  r2_mk_df <- melt(
    lapply(r2_mk, function(x)
      melt(as.matrix(x), varnames = c("factor", "view"))
    ), id.vars=c("factor", "view", "value")
  )
  colnames(r2_mk_df)[ncol(r2_mk_df)] <- "group"
  
  # Subset factors for plotting
  if ((length(factors) == 1) && (factors[1] == "all")) {
    factors <- factors_names(object)
  } else {
    if (is.numeric(factors)) {
      factors <- factors_names(object)[factors]
    } else { 
      stopifnot(all(factors %in% factors_names(object)))
    }
    r2_mk_df <- r2_mk_df[r2_mk_df$factor %in% factors,]
  }
  
  r2_mk_df$factor <- factor(r2_mk_df$factor, levels = factors)
  r2_mk_df$group <- factor(r2_mk_df$group, levels = groups_names(object))
  r2_mk_df$view <- factor(r2_mk_df$view, levels = views_names(object))
  
  # Detect whether to split by group or by view
  groups <- names(r2_list$r2_total)
  views <- colnames(r2_list$r2_per_factor[[1]])
  
  # Set R2 limits
  if (!is.null(min_r2)) r2_mk_df$value[r2_mk_df$value<min_r2] <- 0.001
  min_r2 = 0
  
  if (!is.null(max_r2)) {
    r2_mk_df$value[r2_mk_df$value>max_r2] <- max_r2
  } else {
    max_r2 = max(r2_mk_df$value)
  }
  
  
  # Grid plot with the variance explained per factor and view/group
  p1 <- ggplot(r2_mk_df, aes_string(x=x, y=y)) + 
    geom_tile(aes_string(fill="value"), color="black") +
    facet_wrap(as.formula(sprintf('~%s',split_by)), nrow=1) +
    labs(x="", y="", title="") +
    scale_fill_gradientn(colors=c("gray97","darkblue"), guide="colorbar", limits=c(min_r2,max_r2)) +
    guides(fill=guide_colorbar("Var. (%)")) +
    theme(
      axis.text.x = element_text(size=rel(1.0), color="black"),
      axis.text.y = element_text(size=rel(1.1), color="black"),
      axis.line = element_blank(),
      axis.ticks =  element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size=rel(1.0))
    )
  
  if (isFALSE(legend)) p1 <- p1 + theme(legend.position = "none")
  
  # remove facet title
  if (length(unique(r2_mk_df[,split_by]))==1) p1 <- p1 + theme(strip.text = element_blank())
  
  # Add total variance explained bar plots
  if (plot_total) {
    
    r2_m_df <- melt(lapply(r2_list$r2_total, function(x) lapply(x, function(z) z)),
                    varnames=c("view", "group"), value.name="R2")
    colnames(r2_m_df)[(ncol(r2_m_df)-1):ncol(r2_m_df)] <- c("view", "group")
    
    r2_m_df$group <- factor(r2_m_df$group, levels = MOFA2::groups_names(object))
    r2_m_df$view <- factor(r2_m_df$view, levels = views_names(object))
    
    # Barplots for total variance explained
    min_lim_bplt <- min(0, r2_m_df$R2)
    max_lim_bplt <- max(r2_m_df$R2)
    
    # Barplot with variance explained per view/group (across all factors)
    p2 <- ggplot(r2_m_df, aes_string(x=x, y="R2")) + 
      # ggtitle(sprintf("%s\nTotal variance explained per %s", i, x)) +
      geom_bar(stat="identity", fill="deepskyblue4", color="black", width=0.9) +
      facet_wrap(as.formula(sprintf('~%s',split_by)), nrow=1) +
      xlab("") + ylab("Variance explained (%)") +
      scale_y_continuous(limits=c(min_lim_bplt, max_lim_bplt), expand=c(0.005, 0.005)) +
      theme(
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size=rel(1.1), color="black"),
        axis.text.y = element_text(size=rel(1.0), color="black"),
        axis.title.y = element_text(size=rel(1.0), color="black"),
        axis.line = element_line(size=rel(1.0), color="black"),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=rel(1.0))
      )
    
    # remove facet title
    if (length(unique(r2_m_df[,split_by]))==1) p2 <- p2 + theme(strip.text = element_blank())
    
    # Bind plots      
    plot_list <- list(p1,p2)
    
  } else {
    plot_list <- p1
  }
  
  return(plot_list)
}


#' @title Plot variance explained by the model for a set of features
#' 
#' Returns a tile plot with a group on the X axis and a feature along the Y axis
#' 
#' @name plot_variance_explained_per_feature
#' @param object a \code{\link{MOFA}} object.
#' @param view a view name or index.
#' @param features a vector with indices or names for features from the respective view, 
#' or number of top features to be fetched by their loadings across specified factors. 
#' "all" to plot all features.
#' @param split_by_factor logical indicating whether to split R2 per factor or plot R2 jointly
#' @param group_features_by column name of features metadata to group features by
#' @param groups a vector with indices or names for sample groups (default is all)
#' @param factors a vector with indices or names for factors (default is all)
#' @param min_r2 minimum variance explained for the color scheme (default is 0).
#' @param max_r2 maximum variance explained for the color scheme.
#' @param legend logical indicating whether to add a legend to the plot  (default is TRUE).
#' @param return_data logical indicating whether to return the data frame to plot instead of plotting
#' @param ... extra arguments to be passed to \code{\link{calculate_variance_explained}}
#' @return ggplot object
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom stats as.formula
#' @importFrom reshape2 melt
#' @export
#' @examples 
#' # Using an existing trained model
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' plot_variance_explained_per_feature(model, view = 1)

plot_variance_explained_per_feature <- function(object, view, features = 10,
                                                split_by_factor = FALSE, group_features_by = NULL,
                                                groups = "all", factors = "all",
                                                min_r2 = 0, max_r2 = NULL, legend = TRUE,
                                                return_data = FALSE, ...) {
  
  # Check that one view is requested
  view  <- .check_and_get_views(object, view)
  if (length(view) != 1)
    stop("Please choose a single view to plot features from")
  
  # Fetch loadings, factors, and data  
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Fetch relevant features)
  if (is.numeric(features) && (length(features) == 1)) {
    features <- as.integer(features)
    features <- .get_top_features_by_loading(object, view = view, factors = factors, nfeatures = features)
  } else if (is.character(features)) {
    if (features[1]=="all") features <- 1:object@dimensions$D[[view]]
  }
  features <- .check_and_get_features_from_view(object, view = view, features)
  
  # Collect relevant expectations
  groups <- .check_and_get_groups(object, groups)
  factors <- .check_and_get_factors(object, factors)
  # 1. Loadings: choose a view, one or multiple factors, and subset chosen features
  W <- get_weights(object, views = view, factors = factors)
  W <- lapply(W, function(W_m) W_m[rownames(W_m) %in% features,,drop=FALSE])
  # 2. Factor values: choose one or multiple groups and factors
  Z <- get_factors(object, groups = groups, factors = factors)
  # 3. Data: Choose a view, one or multiple groups, and subset chosen features
  # Y <- lapply(get_expectations(object, "Y")[view], function(Y_m) lapply(Y_m[groups], t))
  Y <- lapply(get_data(object, add_intercept = FALSE)[view], function(Y_m) lapply(Y_m[groups], t))
  Y <- lapply(Y, function(Y_m) lapply(Y_m, function(Y_mg) Y_mg[,colnames(Y_mg) %in% features,drop=FALSE]))
  
  # Replace masked values on Z by 0 (so that they do not contribute to predictions)
  for (g in groups) {
    Z[[g]][is.na(Z[[g]])] <- 0
  }
  
  m <- view  # Use shorter notation when calculating R2
  
  if (split_by_factor) {
    
    # Calculate coefficient of determination per group, factor and feature
    r2_gdk <- lapply(groups, function(g) {
      r2_g <- sapply(features, function(d) { 
        sapply(factors, function(k) {
        a <- sum((as.matrix(Y[[m]][[g]][,d,drop=FALSE]) - tcrossprod(Z[[g]][,k,drop=FALSE], W[[m]][d,k,drop=FALSE]))**2, na.rm = TRUE)
        b <- sum(Y[[m]][[g]][,d,drop=FALSE]**2, na.rm = TRUE)
        return(1 - a/b)
      })
      })
      r2_g <- matrix(r2_g, ncol = length(features), nrow = length(factors))
      colnames(r2_g) <- features
      rownames(r2_g) <- factors
      # Lower bound is zero
      r2_g[r2_g < 0] <- 0
      r2_g
    })
    names(r2_gdk) <- groups
    
    # Convert matrix to long data frame for ggplot2
    r2_gdk_df <- do.call(rbind, r2_gdk)
    r2_gdk_df <- data.frame(r2_gdk_df, 
                            "group" = rep(groups, lapply(r2_gdk, nrow)),
                            "factor" = rownames(r2_gdk_df))
    r2_gdk_df <- melt(r2_gdk_df, id.vars = c("group", "factor"))
    colnames(r2_gdk_df) <- c("group", "factor", "feature", "value")
    
    r2_gdk_df$group <- factor(r2_gdk_df$group, levels = unique(r2_gdk_df$group))
    
    r2_df <- r2_gdk_df
    
  } else {
    
    # Calculate coefficient of determination per group and feature
    r2_gd <- lapply(groups, function(g) {
      r2_g <- lapply(features, function(d) {
        a <- sum((as.matrix(Y[[m]][[g]][,d,drop=FALSE]) - tcrossprod(Z[[g]], W[[m]][d,,drop=FALSE]))**2, na.rm = TRUE)
        b <- sum(Y[[m]][[g]][,d,drop=FALSE]**2, na.rm = TRUE)
        return(1 - a/b)
      })
      names(r2_g) <- features
      # Lower bound is zero
      r2_g[r2_g < 0] <- 0
      r2_g
    })
    names(r2_gd) <- groups
    
    # Convert matrix to long data frame for ggplot2
    tmp <- as.matrix(data.frame(lapply(r2_gd, unlist))) 
    colnames(tmp) <- groups
    r2_gd_df <- melt(tmp)
    colnames(r2_gd_df) <- c("feature", "group", "value")
    
    r2_gd_df$group <- factor(r2_gd_df$group, levels = unique(r2_gd_df$group))
    
    r2_df <- r2_gd_df
    
  }
  
  # Transform from fraction to percentage
  r2_df$value <- 100*r2_df$value
  
  # Calculate minimum R2 to display
  if (!is.null(min_r2)) {
    r2_df$value[r2_df$value<min_r2] <- 0.001
  }
  min_r2 <- 0
  
  # Calculate maximum R2 to display
  if (!is.null(max_r2)) {
    r2_df$value[r2_df$value>max_r2] <- max_r2
  } else {
    max_r2 <- max(r2_df$value)
  }
  
  # Group features
  if (!is.null(group_features_by)) {
    features_indices <- match(r2_df$feature, features_metadata(object)$feature)
    features_grouped <- features_metadata(object)[,group_features_by,drop=FALSE][features_indices,,drop=FALSE]
    # If features grouped using multiple variables, concatenate them
    if (length(group_features_by) > 1) {
      features_grouped <- apply(features_grouped, 1, function(row) paste0(row, collapse="_"))
    } else {
      features_grouped <- features_grouped[,group_features_by,drop=TRUE]
    }
    r2_df["feature_group"] <- features_grouped
  }
  
  if (return_data)
    return(r2_df)
  
  if (split_by_factor) {
    r2_df$factor <- factor(r2_df$factor, levels = factors_names(object))
  }
  
  # Grid plot with the variance explained per feature in every group
  p <- ggplot(r2_df, aes_string(x = "group", y = "feature")) + 
    geom_tile(aes_string(fill = "value"), color = "black") +
    guides(fill = guide_colorbar("R2 (%)")) +
    labs(x = "", y = "", title = "") +
    scale_fill_gradientn(colors=c("gray97","darkblue"), guide="colorbar", limits=c(min_r2, max_r2)) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 12),
      axis.line = element_blank(),
      axis.ticks =  element_blank(),
      strip.text = element_text(size = 12),
    )
  
  if (!is.null(group_features_by) && split_by_factor) {
    p <- p + facet_grid(feature_group ~ factor, scales = "free_y")
  } else if (split_by_factor) {
    p <- p + facet_wrap(~factor, nrow = 1)
  } else if (!is.null(group_features_by)) {
    p <- p + facet_wrap(~feature_group, ncol = 1, scales = "free")
  }
  
  if (!legend)
    p <- p + theme(legend.position = "none")
  
  return(p)
}
