#' @title Calculate variance explained by the model
#' @name calculate_variance_explained
#' @param object a \code{\link{BioFAModel}} object.
#' @param views character vector with the view names, or numeric vector with view indexes. Default is 'all'
#' @param factors character vector with the factor names, or numeric vector with the factor indexes. Default is 'all'
#' @param groups character vector with the group names, or numeric vector with group indexes. Default is 'all'
#' @details This function takes a trained BioFAModel as input and calculates for each view the coefficient of determination (R^2),
#' i.e. the proportion of variance in the data explained by the BioFAM factor(s) (both jointly and for each individual factor).
#' In case of non-Gaussian data the variance explained on the Gaussian pseudo-data is calculated.
#' @return a list with matrices with the amount of variation explained per factor and view, and optionally total variance explained per view and variance explained by each feature alone
#' @import DelayedArray
#' @export
calculate_variance_explained <- function(object, views = "all", groups = "all", factors = "all", ...) {

  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")

  # Define views and groups
  views  <- .check_and_get_views(object, views)
  groups <- .check_and_get_groups(object, groups)

  # Define factors
  if (paste0(factors, collapse="") == "all") {
    factors <- factors_names(object)
  } else if (is.numeric(factors)) {
      factors <- factors_names(object)[factors]
  } else {
    stopifnot(all(factors %in% factors_names(object)))
  }
  K <- length(factors)

  # Collect relevant expectations
  W <- get_weights(object, views=views, factors=factors)
  Z <- get_factors(object, groups=groups, factors=factors)
  Y <- get_expectations(object, "Y")[views]
  Y <- lapply(Y, function(x) lapply(x,t))

  # Replace masked values on Z by 0 (so that they do not contribute to predictions)
  for (g in groups) {
    Z[[g]][is.na(Z[[g]])] <- 0
  }

  # Check that data observations are centered
  # for (m in views) { 
  #   for (g in groups) {
  #     if (object@model_options$likelihoods[m]=="gaussian" & 
  #         !all(DelayedArray::colMeans(Y[[m]][[g]], na.rm = TRUE) < 1e-2, na.rm = TRUE))
  #       cat(sprintf("Warning: data for view %s and group %s is not centered\n", m, g))
  #   }
  # }

  Y <- .name_views_and_groups(Y, views, groups)

  # Calculate coefficient of determination per group and view
  r2_m <- tryCatch({
    lapply(groups, function(g) lapply(views, function(m) {
        # a <- sum((as.matrix(Y[[m]][[g]]) - DelayedArray::tcrossprod(Z[[g]], W[[m]]))**2, na.rm = TRUE)
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

  # Calculate coefficient of determination per group, factor and view
  r2_mk <- lapply(groups, function(g) {
    tmp <- sapply(views, function(m) { sapply(factors, function(k) {
        # a <- sum((as.matrix(Y[[m]][[g]]) - DelayedArray::tcrossprod(Z[[g]][,k], W[[m]][,k]))**2, na.rm = TRUE)
        a <- sum((as.matrix(Y[[m]][[g]]) - tcrossprod(Z[[g]][,k], W[[m]][,k]))**2, na.rm = TRUE)
        b <- sum(Y[[m]][[g]]**2, na.rm = TRUE)
        return(1 - a/b)
      })
    })
    tmp <- matrix(tmp, ncol = length(views), nrow = length(factors))
    colnames(tmp) <- views; rownames(tmp) <- factors
    return(tmp)
  }); names(r2_mk) <- groups

  # Store results
  # r2_mk = lapply(r2_mk, function(x){x[x < 0] = 0; return(x)})
  r2_list <- list(r2_total = r2_m, r2_per_factor = r2_mk)

  return(r2_list)
}









#' @title Plot variance explained by the model
#' 
#' Returns a list of plots with specifies axes.
#' 
#' Consider using cowplot::plot_grid(plotlist = ...) in order to combine plots.
#' 
#' @name plot_variance_explained
#' @param object a \code{\link{MOFAmodel}} object
#' @param x string specifying the x-axis ("view", "factor", or "group")
#' @param y string specifying the y-axis ("view", "factor", or "group")
#' @param split_by string specifying the dimension to be faceted ("view", "factor", or "group")
#' @param cluster logical value indicating whether to do hierarchical clustering on the plot
#' @param factors TO-FILL. Default is "all"
#' @param plot_total logical value to indicate if to plot the total variance explained (for the variable in the x-axis)
#' @param legend logical indicating whether to add a legend
#' @param ... extra arguments to be passed to \code{\link{calculate_variance_explained}}
#' @return ggplot object
#' @import pheatmap ggplot2 reshape2
#' @importFrom cowplot plot_grid
#' @export
plot_variance_explained <- function(object, x = "view", y = "factor", split_by = NA, plot_total = FALSE, 
                                    factors = "all", min_r2=0, legend = TRUE, ...) {
  
  # Sanity checks 
  if (length(unique(c(x, y, split_by))) != 3) { 
    stop(paste0("Please ensure x, y, and split_by arguments are different.\n",
                "  Possible values are `view`, `group`, and `factor`."))
  }

  # Automatically fill split_by in
  if (is.na(split_by)) split_by <- setdiff(c("view", "factor", "group"), c(x, y, split_by))
  
  # Calculate variance explained
  if (.hasSlot(object, "cache") && ("variance_explained" %in% names(object@cache))) {
    r2_list <- object@cache$variance_explained
  } else {
    r2_list <- calculate_variance_explained(object, ...)
  }

  r2_m <- r2_list$r2_total
  r2_mk <- r2_list$r2_per_factor
  # r2_m  <- lapply(r2_list$r2_total[groups], function(e) e[views])
  # r2_mk <- lapply(r2_list$r2_per_factor[groups], function(e) e[,views])

  # convert matrix to long data frame for ggplot2
  r2_mk_df <- reshape2::melt(
    lapply(r2_mk, function(x)
      reshape2::melt(as.matrix(x), varnames = c("factor", "view"))
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

  r2_m_df <- reshape2::melt(lapply(r2_m, function(x) lapply(x, function(z) z)),
                              varnames=c("view", "group"), value.name="R2")
  colnames(r2_m_df)[(ncol(r2_m_df)-1):ncol(r2_m_df)] <- c("view", "group")

  # sort views according to hierarchical clustering on the variance explained pattern
  # g <- which.max(sapply(r2_m, function(x) sum(unlist(x)))) # use group with the highest variance explained
  # if (cluster & ncol(r2_mk[[g]])>1) {
  #   hc <- hclust(dist(t(r2_mk[[g]])))
  #   r2_mk_df$view <- factor(r2_mk_df$view, levels = colnames(r2_mk[[g]])[hc$order])
  #   r2_m_df$view <- factor(r2_m_df$view, levels = colnames(r2_mk[[g]])[hc$order])
  # }

  # Heatmaps (grid plots) for variance explained per factor
  # min_lim_p1 <- min(r2_mk_df$value)
  min_lim_p1 <- 0
  max_lim_p1 <- max(r2_mk_df$value)

  # Barplots for total variance explained
  min_lim_bplt <- min(0, r2_m_df$R2)
  max_lim_bplt <- max(r2_m_df$R2)

  # Detect whether to split by group or by view
  groups <- names(r2_list$r2_total)
  views <- colnames(r2_list$r2_per_factor[[1]])

  r2_mk_df$value[r2_mk_df$value<min_r2] <- 0.00001
  
  # Grid plot with the variance explained per factor and view/group
  p1 <- ggplot(r2_mk_df, aes_string(x=x, y=y)) + 
    geom_tile(aes(fill=value), color="black") +
    facet_wrap(as.formula(sprintf('~%s',split_by)), nrow=1) +
    guides(fill=guide_colorbar("R2")) +
    labs(x="", y="", title="") +
    scale_fill_gradientn(colors=c("gray97","darkblue"), guide="colorbar", limits=c(min_lim_p1, max_lim_p1)) +
    guides(fill=guide_colorbar("R2")) +
    theme(
      axis.title.x = element_blank(),
      # axis.text.x = element_text(size=12, angle=30, hjust=0.5, vjust=1, color="black"),
      axis.text.x = element_text(size=12, color="black"),
      axis.text.y = element_text(size=12, color="black"),
      axis.title.y = element_text(size=14),
      axis.line = element_blank(),
      axis.ticks =  element_blank(),
      panel.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size=12)
    )
  
    if (!legend)
      p1 <- p1 + theme(legend.position = "none")
    
    # Add total variance explained bar plots
    if (plot_total) {

      # Barplot with variance explained per view/group (across all factors)
      p2 <- ggplot(r2_m_df, aes_string(x=x, y="R2")) + 
        # ggtitle(sprintf("%s\nTotal variance explained per %s", i, x)) +
        geom_bar(stat="identity", fill="deepskyblue4", width=0.9) +
        facet_wrap(as.formula(sprintf('~%s',split_by)), nrow=1) +
        xlab("") + ylab("Variance explained (R^2)") +
        scale_y_continuous(limits=c(min_lim_bplt, max_lim_bplt), expand=c(0.005, 0.005)) +
        theme(
          plot.title = element_text(size=17, hjust=0.5),
          axis.ticks.x = element_blank(),
          axis.text.x = element_text(size=12, color="black"),
          axis.text.y = element_text(size=rel(1.0), color="black"),
          axis.title.y = element_text(size=13, color="black"),
          axis.line = element_line(size=rel(1.0), color="black"),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size=13)
        )
      
      plot_list <- list(p1,p2)
      
      # plot_list <- plot_grid(
      #   plotlist = list(p2, p1), 
      #   align = "v", 
      #   nrow = 2, ncol = 1, 
      #   rel_heights = c(1/3,2/3), 
      #   axis="l"
      # )
      
    } else {
      plot_list <- p1
    }
  
  return(plot_list)
}
