#' @title Calculate variance explained by the model
#' @name calculate_variance_explained
#' @param object a \code{\link{BioFAModel}} object.
#' @param views character vector with the view names, or numeric vector with view indexes. Default is 'all'
#' @param factors character vector with the factor names, or numeric vector with the factor indexes. Default is 'all'
#' @param groups character vector with the group names, or numeric vector with group indexes. Default is 'all'
#' @details This function takes a trained BioFAModel as input and calculates for each view the coefficient of determination (R2),
#' i.e. the proportion of variance in the data explained by the BioFAM factor(s) (both jointly and for each individual factor).
#' In case of non-Gaussian data the variance explained on the Gaussian pseudo-data is calculated.
#' @return a list with matrices with the amount of variation explained per factor and view, and optionally total variance explained per view and variance explained by each feature alone
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
  W <- get_weights(object, views, factors)
  Z <- get_factors(object, groups, factors)
  Y <- get_expectations(object, "Y")  # for non-Gaussian likelihoods the pseudodata is considered
  Y <- lapply(Y, function(x) lapply(x,t))

  # Replace masked values on Z by 0 (so that they do not contribute to predictions)
  for (p in groups) {
    Z[[p]][is.na(Z[[p]])] <- 0
  }

  for (m in views) { for (p in groups) {
    if (!all(colMeans(Y[[m]][[p]],na.rm=T)<1e-2,na.rm=T))
      cat(sprintf("Warning: data for view %s is not centered\n",m))
  }}

  Y <- .name_views_and_groups(Y, views, groups)

  # Calculate coefficient of determination per group and view
  fvar_m <- lapply(groups, function(p)
    lapply(views, function(m) {
      a <- sum((Y[[m]][[p]]-tcrossprod(Z[[p]],W[[m]]))**2, na.rm=T)
      b <- sum(Y[[m]][[p]]**2, na.rm=T)
      return(1 - a/b)
    }
    ))
  fvar_m <- .name_views_and_groups(fvar_m, groups, views)

  # Calculate coefficient of determination per group, factor and view
  fvar_mk <- lapply(groups, function(p) {
    tmp <- sapply(views, function(m) {
      sapply(factors, function(k) {
        a <- sum((Y[[m]][[p]]-tcrossprod(Z[[p]][,k],W[[m]][,k]))**2, na.rm=T)
        b <- sum(Y[[m]][[p]]**2, na.rm=T)
        return(1 - a/b)
      })
    })
    tmp <- matrix(tmp, ncol=length(views), nrow=length(factors))
    colnames(tmp) <- views; rownames(tmp) <- factors
    return(tmp)
  }); names(fvar_mk) <- groups

  # Store results
  # fvar_mk = lapply(fvar_mk, function(x){x[x < 0] = 0; return(x)})
  r2_list <- list(r2_total = fvar_m, r2_per_factor = fvar_mk)

  return(r2_list)
}









#' @title Plot variance explained by the model
#' @name plot_variance_explained
#' @param object a \code{\link{MOFAmodel}} object.
#' @param cluster logical indicating whether to do hierarchical clustering on the plot
#' @param ... extra arguments to be passed to \code{\link{calculate_variance_explained}}
#' @return ggplot object
#' @import pheatmap ggplot2 reshape2
#' @importFrom cowplot plot_grid
#' @export
plot_variance_explained <- function(object, cluster = TRUE, ...) {

  # Calculate variance explained
  if (.hasSlot(object, "cache") && ("variance_explained" %in% names(object@cache))) {
    r2_list <- object@cache[["variance_explained"]]
  } else {
    r2_list <- calculate_variance_explained(object, ...)
  }

  fvar_m <- r2_list$r2_total
  fvar_mk <- r2_list$r2_per_factor
  # fvar_m  <- lapply(r2_list$r2_total[groups], function(e) e[views])
  # fvar_mk <- lapply(r2_list$r2_per_factor[groups], function(e) e[,views])

  # convert matrix to long data frame for ggplot2
  fvar_mk_df <- reshape2::melt(
    lapply(fvar_mk, function(x)
      reshape2::melt(as.matrix(x), varnames = c("factor", "view"))
    ), id.vars=c("factor", "view", "value")
  )
  colnames(fvar_mk_df)[ncol(fvar_mk_df)] <- "group"
  fvar_mk_df$factor <- factor(fvar_mk_df$factor)
  fvar_mk_df$group <- factor(fvar_mk_df$group)

  fvar_m_df <- reshape2::melt(lapply(fvar_m, function(x) lapply(x, function(z) z)),
                              varnames=c("view", "group"), value.name="R2")
  colnames(fvar_m_df)[(ncol(fvar_m_df)-1):ncol(fvar_m_df)] <- c("view", "group")

  # sort views according to hierarchical clustering on the variance explained pattern
  g <- which.max(sapply(fvar_m, function(x) sum(unlist(x)))) # use group with the highest variance explained
  if (cluster & ncol(fvar_mk[[g]])>1) {
    hc <- hclust(dist(t(fvar_mk[[g]])))
    fvar_mk_df$view <- factor(fvar_mk_df$view, levels = colnames(fvar_mk[[g]])[hc$order])
    fvar_m_df$view <- factor(fvar_m_df$view, levels = colnames(fvar_mk[[g]])[hc$order])
  }

  # Heatmaps (grid plots) for fvar_mk
  hms   <- list()
  min_lim_hm <- min(fvar_mk_df$value)
  max_lim_hm <- max(fvar_mk_df$value)

  # Barplots for fvar_m
  bplts <- list()
  min_lim_bplt <- min(0, fvar_m_df$R2)
  max_lim_bplt <- max(fvar_m_df$R2)

  # Detect whether to split by group or by view
  groups <- names(r2_list$r2_total)
  views <- colnames(r2_list$r2_per_factor[[1]])
  x="view";  split_by="group"
  if ( length(groups)>1 ) { x="group"; split_by="view" }

  plot_list <- list()
  for (i in levels(fvar_mk_df[[split_by]])) {

    # Barplot with variance explained per view/group (across all factors)
    m_title <- sprintf("%s\nTotal variance explained per %s", i, x)
    bplt <- ggplot(fvar_m_df[fvar_m_df[[split_by]] == i,], aes_string(x=x, y="R2")) + 
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
    # bplts[[i]] <- bplt
    
    # Grid plot with the variance explained per factor and view/group
    mk_title <- paste0("Variance explained per factor")
    hm <- ggplot(fvar_mk_df[fvar_mk_df[[split_by]] == i,], aes_string(x=x, y="factor")) + 
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
    hm <- hm + ggtitle(mk_title) + guides(fill=guide_colorbar("R2"))
    # hms[[i]] <- hm
    
    # Join the two plots
    plot_list[[i]] <- plot_grid(
      plotlist = list(bplt, hm), 
      align = "v", 
      nrow = 2, ncol = 1, 
      rel_heights = c(1/3,2/3), 
      axis="l"
    )
    # print(plot_list[[i]])
  }
  return(plot_list)
}
