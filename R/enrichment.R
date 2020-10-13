##########################################################
## Functions to perform Feature Set Enrichment Analysis ##
##########################################################

#' @title Run feature set Enrichment Analysis
#' @name run_enrichment 
#' @description Method to perform feature set enrichment analysis. Here we use a slightly modified version of the \link[PCGSE]{pcgse} function.
#' @param object a \code{\link{MOFA}} object.
#' @param view a character with the view name, or a numeric vector with the index of the view to use.
#' @param feature.sets data structure that holds feature set membership information. 
#' Must be a binary membership matrix (rows are feature sets and columns are features). See details below for some pre-built gene set matrices.
#' @param factors character vector with the factor names, or numeric vector with the index of the factors for which to perform the enrichment.
#' @param set.statistic the set statisic computed from the feature statistics. Must be one of the following: "mean.diff" (default) or "rank.sum".
#' @param statistical.test the statistical test used to compute the significance of the feature set statistics under a competitive null hypothesis.
#' Must be one of the following: "parametric" (default), "cor.adj.parametric", "permutation".
#' @param sign use only "positive" or "negative" weights. Default is "all".
#' @param min.size Minimum size of a feature set (default is 10).
#' @param nperm number of permutations. Only relevant if statistical.test is set to "permutation". Default is 1000
#' @param p.adj.method Method to adjust p-values factor-wise for multiple testing. Can be any method in p.adjust.methods(). Default uses Benjamini-Hochberg procedure.
#' @param alpha FDR threshold to generate lists of significant pathways. Default is 0.1
#' @param verbose boolean indicating whether to print messages on progress 
#' @details 
#'  The aim of this function is to relate each factor to pre-defined biological pathways by performing a gene set enrichment analysis on the feature weights. \cr
#'  This function is particularly useful when a factor is difficult to characterise based only on the genes with the highest weight. \cr
#'  We provide a few pre-built gene set matrices in the MOFAdata package. See \code{https://github.com/bioFAM/MOFAdata} for details. \cr
#'  The function we implemented is based on the \code{\link[PCGSE]{pcgse}} function with some modifications. 
#'  Please read this paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4543476 for details on the math.
#' @return a list with five elements: 
#' \item{\strong{pval}:}{ matrices with nominal p-values. }
#' \item{\strong{pval.adj}:}{ matrices with FDR-adjusted p-values. }
#' \item{\strong{feature.statistics}:}{ matrices with the local (feature-wise) statistics.  }
#' \item{\strong{set.statistics}:}{ matrices with the global (gene set-wise) statistics.  }
#' \item{\strong{sigPathways}}{ list with significant pathways per factor. }
#' @importFrom stats p.adjust var p.adjust.methods
#' @export

run_enrichment <- function(object, view, feature.sets, factors = "all",
                           set.statistic = c("mean.diff", "rank.sum"),
                           statistical.test = c("parametric", "cor.adj.parametric", "permutation"), sign = c("all","positive","negative"),
                           min.size = 10, nperm = 1000, p.adj.method = "BH", alpha = 0.1, verbose = TRUE) {
  
  # Quality control
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (!(is(feature.sets, "matrix") & all(feature.sets %in% c(0,1)))) stop("feature.sets has to be a list or a binary matrix.")
  
  # Define views
  view  <- .check_and_get_views(object, view)
  
  # Define factors
  factors  <- .check_and_get_factors(object, factors)
  
  # Parse inputs
  sign <- match.arg(sign)
  set.statistic <- match.arg(set.statistic)
  statistical.test <- match.arg(statistical.test)
  
  # Collect observed data
  data <- get_data(object, views = view, as.data.frame = FALSE)[[1]]
  if(is(data, "list")) data <- Reduce(cbind, data) # concatenate groups
  data <- t(data)
  
  # Collect relevant expectations
  W <- get_weights(object, views=view, factors=factors, scale = TRUE)[[1]]
  Z <- get_factors(object, factors=factors)
  if(is(Z, "list")) Z <- Reduce(rbind, Z)
  stopifnot(rownames(data) == rownames(Z))
  
  # Remove features with no variance
  # if (statistical.test %in% c("cor.adj.parametric")) {
  idx <- apply(data,2, function(x) var(x,na.rm=TRUE))==0
  if (sum(idx)>=1) {
    warning(sprintf("%d fetures were removed because they had no variance in the data.\n",sum(idx)))
    data <- data[,!idx]
    W <- W[!idx,]
  }
  
  # Check if some features do not intersect between the feature sets and the observed data and remove them
  features <- intersect(colnames(data),colnames(feature.sets))
  if(length(features)== 0) stop("Feature names in feature.sets do not match feature names in model.")
  if(verbose) {
    message(sprintf("Intersecting features names in the model and the gene set annotation results in a total of %d features.",length(features)))
  }
  data <- data[,features]
  W <- W[features,,drop=FALSE]
  feature.sets <- feature.sets[,features]
  
  # Filter feature sets with small number of features
  feature.sets <- feature.sets[rowSums(feature.sets)>=min.size,]
  
  # Subset weights by sign
  if (sign=="positive") {
    W[W<0] <- NA
  } else if (sign=="negative") {
    W[W>0] <- NA
    W <- abs(W)
  }
  
  # Print options
  if(verbose) {
    message("\nRunning feature set Enrichment Analysis with the following options...\n",
            sprintf("View: %s \n", view),
            sprintf("Number of feature sets: %d \n", nrow(feature.sets)),
            sprintf("Set statistic: %s \n", set.statistic),
            sprintf("Statistical test: %s \n", statistical.test)
    )
    if (sign%in%c("positive","negative"))
      message(sprintf("Subsetting weights with %s sign",sign))
    if (statistical.test=="permutation") {
      message(sprintf("Number of permutations: %d", nperm))
    }
    message("\n")
  }
  
  if (nperm<100) 
    warning("A large number of permutations (at least 1000) is required for the permutation approach!\n")
  
  # Non-parametric permutation test
  if (statistical.test == "permutation") {

        null_dist_tmp <- lapply(seq_len(nperm), function(i) {
      print(sprintf("Running permutation %d/%d...",i,nperm))
      perm <- sample(ncol(data))
      
      # Permute rows of the weight matrix to obtain a null distribution
      W_null <- W[perm,]
      rownames(W_null) <- rownames(W)
      colnames(W_null) <- colnames(W)
      
      # Permute columns of the data matrix correspondingly (only matters for cor.adjusted test)
      data_null <- data[,perm]
      rownames(data_null) <- rownames(data)
      
      # Compute null (or background) statistic
      s.background <- .pcgse(
        data = data_null, 
        prcomp.output = list(rotation=W_null, x=Z),
        pc.indexes = seq_along(factors), 
        feature.sets = feature.sets,
        set.statistic = set.statistic,
        set.test = "parametric")$statistic
      return(abs(s.background))
    })
    null_dist <- do.call("rbind", null_dist_tmp)
    colnames(null_dist) <- factors
    
    # Compute foreground statistics
    results <- .pcgse(
      data = data, 
      prcomp.output = list(rotation=W, x=Z),
      pc.indexes = seq_along(factors), 
      feature.sets = feature.sets,
      set.statistic = set.statistic,
      set.test = "parametric")
    s.foreground <- results$statistic
    
    # Calculate p-values based on fraction true statistic per factor and feature set is larger than permuted
    xx <- array(unlist(null_dist_tmp), dim = c(nrow(null_dist_tmp[[1]]), ncol(null_dist_tmp[[1]]), length(null_dist_tmp)))
    ll <- lapply(seq_len(nperm), function(i) xx[,,i] > abs(s.foreground))
    results$p.values <- Reduce("+",ll)/nperm
    
    # Parametric test
  } else {
    results <- .pcgse(
      data = data,
      prcomp.output = list(rotation=W, x=Z),
      pc.indexes = seq_along(factors),
      feature.sets = feature.sets,
      set.statistic = set.statistic,
      set.test = statistical.test
    )
  }
  
  # Parse results
  pathways <- rownames(feature.sets)
  colnames(results$p.values) <- colnames(results$statistics) <- colnames(results$feature.statistics) <- factors
  rownames(results$p.values) <- rownames(results$statistics) <- pathways
  rownames(results$feature.statistics) <- colnames(data)
  
  # adjust for multiple testing
  if(!p.adj.method %in%  p.adjust.methods) 
    stop("p.adj.method needs to be an element of p.adjust.methods")
  adj.p.values <- apply(results$p.values, 2,function(lfw) p.adjust(lfw, method = p.adj.method))
  
  # If we specify a direction, we are only interested in overrepresented pathawys in the selected direction
  if (sign%in%c("positive","negative")) {
    results$p.values[results$statistics<0] <- 1.0
    adj.p.values[results$statistics<0] <- 1.0
    results$statistics[results$statistics<0] <- 0
  }
  
  
  # obtain list of significant pathways
  sigPathways <- lapply(factors, function(j) rownames(adj.p.values)[adj.p.values[,j] <= alpha])
  
  # prepare output
  output <- list(
    feature.sets = feature.sets, 
    pval = results$p.values, 
    pval.adj = adj.p.values, 
    feature.statistics = results$feature.statistics,
    set.statistics = results$statistics,
    sigPathways = sigPathways
  )
  return(output)
}


########################
## Plotting functions ##
########################


#' @title Plot output of gene set Enrichment Analysis
#' @name plot_enrichment
#' @description Method to plot the results of the gene set Enrichment Analyisis
#' @param enrichment.results output of \link{run_enrichment} function
#' @param factor a string with the factor name or an integer with the factor index
#' @param alpha p.value threshold to filter out gene sets
#' @param max.pathways maximum number of enriched pathways to display
#' @param text_size text size
#' @param dot_size dot size
#' @details it requires \code{\link{run_enrichment}} to be run beforehand.
#' @return a \code{ggplot2} object
#' @import ggplot2
#' @importFrom utils head
#' @export
plot_enrichment <- function(enrichment.results, factor, alpha = 0.1, max.pathways = 25,
                            text_size = 1.0, dot_size = 5.0) {
  
  # Sanity checks
  stopifnot(is.numeric(alpha)) 
  stopifnot(length(factor)==1) 
  if (is.numeric(factor)) factor <- colnames(enrichment.results$pval.adj)[factor]
  if(!factor %in% colnames(enrichment.results$pval)) 
    stop(paste0("No gene set enrichment calculated for factor ", factor))
  
  # get p-values
  p.values <- enrichment.results$pval.adj
  
  # Get data  
  tmp <- data.frame(
    pvalues = p.values[,factor, drop=TRUE], 
    pathway = rownames(p.values)
  )
  
  # Filter out pathways
  tmp <- tmp[tmp$pvalue<=alpha,,drop=FALSE]
  if (nrow(tmp)==0) stop("No siginificant pathways at the specified alpha threshold")
  
  # If there are too many pathways enriched, just keep the 'max_pathways' more significant
  if (nrow(tmp)>max.pathways) tmp <- head(tmp[order(tmp$pvalue),],n=max.pathways)
  
  # Convert pvalues to log scale
  tmp$logp <- -log10(tmp$pvalue+1e-100)
  
  #order according to significance
  tmp$pathway <- factor(tmp$pathway <- rownames(tmp), levels = tmp$pathway[order(tmp$pvalue, decreasing = TRUE)])
  tmp$start <- 0
  
  p <- ggplot(tmp, aes_string(x="pathway", y="logp")) +
    geom_point(size=dot_size) +
    geom_hline(yintercept=-log10(alpha), linetype="longdash") +
    scale_color_manual(values=c("black","red")) +
    geom_segment(aes_string(xend="pathway", yend="start")) +
    ylab("-log pvalue") +
    coord_flip() +
    theme(
      axis.text.y = element_text(size=rel(text_size), hjust=1, color='black'),
      axis.text.x = element_text(size=rel(1.2), vjust=0.5, color='black'),
      axis.title.y=element_blank(),
      legend.position='none',
      panel.background = element_blank()
    )
  
  return(p)
}

#' @title Heatmap of Feature Set Enrichment Analysis results
#' @name plot_enrichment_heatmap
#' @description This method generates a heatmap with the adjusted p.values that
#'  result from the the feature set enrichment analysis. Rows are feature sets and columns are factors.
#' @param enrichment.results output of \link{run_enrichment} function
#' @param alpha FDR threshold to filter out unsignificant feature sets which are
#'  not represented in the heatmap. Default is 0.10.
#' @param cap cap p-values below this threshold
#' @param log_scale logical indicating whether to plot the -log of the p.values.
#' @param ... extra arguments to be passed to the \link{pheatmap} function
#' @return produces a heatmap
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @export
plot_enrichment_heatmap <- function(enrichment.results, alpha = 0.1, cap = 1e-50, log_scale = TRUE, ...) {
  
  # get p-values
  p.values <- enrichment.results$pval.adj
  
  # remove factors that are full of NAs
  p.values <- p.values[,colMeans(is.na(p.values))<1]
  
  # cap p-values 
  p.values[p.values<cap] <- cap
  
  # Apply Log transform
  if (log_scale) {
    p.values <- -log10(p.values+1e-50)
    alpha <- -log10(alpha)
    col <- colorRampPalette(c("lightgrey","red"))(n=100)
  } else {
    col <- colorRampPalette(c("red","lightgrey"))(n=100)
  }
  
  # Generate heatmap
  pheatmap(p.values, color = col, cluster_cols = FALSE, show_rownames = FALSE, ...)
}


#' @title Plot detailed output of the Feature Set Enrichment Analysis
#' @name plot_enrichment_detailed
#' @description Method to plot a detailed output of the Feature Set Enrichment Analyisis (FSEA). \cr
#' Each row corresponds to a significant pathway, sorted by statistical significance, and each dot corresponds to a gene. \cr
#' For each pathway, we display the top genes of the pathway sorted by the corresponding feature statistic (by default, the absolute value of the weight)
#' The top genes with the highest statistic (max.genes argument) are displayed and labeled in black. The remaining genes are colored in grey.
#' @param enrichment.results output of \link{run_enrichment} function
#' @param factor string with factor name or numeric with factor index
#' @param alpha p.value threshold to filter out feature sets
#' @param max.pathways maximum number of enriched pathways to display
#' @param max.genes maximum number of genes to display, per pathway
#' @param text_size size of the text to label the top genes
#' @return a \code{ggplot2} object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom dplyr top_n
#' @importFrom ggrepel geom_text_repel
#' @export
plot_enrichment_detailed <- function(enrichment.results, factor, 
                                     alpha = 0.1, max.genes = 5, max.pathways = 10, text_size = 3) {
  
  # Sanity checks
  stopifnot(is.list(enrichment.results))
  stopifnot(length(factor)==1) 
  if (!is.numeric(factor)) {
    if(!factor %in% colnames(enrichment.results$pval)) 
      stop(paste0("No feature set enrichment calculated for ", factor))
  }
  
  # Fetch and prepare data  
  
  # foo
  foo <- reshape2::melt(enrichment.results$feature.statistics[,factor], na.rm=TRUE, value.name="feature.statistic")
  foo$feature <- rownames(foo)
  
  # bar
  feature.sets <- enrichment.results$feature.sets
  feature.sets[feature.sets==0] <- NA
  bar <- reshape2::melt(feature.sets, na.rm=TRUE)[,c(1,2)]
  colnames(bar) <- c("pathway","feature")
  bar$pathway <- as.character(bar$pathway)
  bar$feature <- as.character(bar$feature)
  
  # baz
  baz <- reshape2::melt(enrichment.results$pval.adj[,factor], value.name="pvalue", na.rm=TRUE)
  baz$pathway <- rownames(baz)
  
  # Filter out pathways by p-values
  baz <- baz[baz$pvalue<=alpha,,drop=FALSE]
  if(nrow(baz)==0) {
    stop("No siginificant pathways at the specified alpha threshold. \n
         For an overview use plot_enrichment_heatmap().")
  } else {
    if (nrow(baz)>max.pathways)
      baz <- head(baz[order(baz$pvalue),],n=max.pathways)
  }
  
  # order pathways according to significance
  baz$pathway <- factor(baz$pathway, levels = baz$pathway[order(baz$pvalue, decreasing = TRUE)])
  
  # Merge
  foobar <- merge(foo, bar, by="feature")
  tmp <- merge(foobar, baz, by="pathway")
  
  # Select the top N features with the largest feature.statistic (per pathway)
  tmp_filt <- top_n(group_by(tmp, pathway), n=max.genes, abs(feature.statistic))
  
  # Add number of features and p-value per pathway
  pathways <- unique(tmp_filt$pathway)
  
  # Add Ngenes and p-values to the pathway name
  df <- data.frame(pathway=pathways, nfeatures=rowSums(feature.sets,na.rm=TRUE)[pathways])
  df <- merge(df, baz, by="pathway")
  df$pathway_long_name <- sprintf("%s\n (Ngenes = %d) \n (p-val = %0.2g)",df$pathway, df$nfeatures, df$pvalue)
  tmp <- merge(tmp, df[,c("pathway","pathway_long_name")], by="pathway")
  tmp_filt <- merge(tmp_filt, df[,c("pathway","pathway_long_name")], by="pathway")
  
  # sort pathways by p-value
  order_pathways <- df$pathway_long_name[order(df$pvalue,decreasing=TRUE) ]
  tmp$pathway_long_name <- factor(tmp$pathway_long_name, levels=order_pathways)
  tmp_filt$pathway_long_name <- factor(tmp_filt$pathway_long_name, levels=order_pathways)
  
  p <- ggplot(tmp, aes_string(x="pathway_long_name", y="feature.statistic")) +
    geom_text_repel(aes_string(x="pathway_long_name", y="feature.statistic", label="feature"), size=text_size, color="black", force=1, data=tmp_filt) +
    geom_point(size=0.5, color="lightgrey") +
    geom_point(aes_string(x="pathway_long_name", y="feature.statistic"), size=1, color="black", data=tmp_filt) +
    labs(x="", y="Weight (scaled)", title="") +
    coord_flip() +
    theme(
      axis.line = element_line(color="black"),
      axis.text.y = element_text(size=rel(0.75), hjust=1, color='black'),
      axis.text.x = element_text(size=rel(1.0), vjust=0.5, color='black'),
      axis.title.y=element_blank(),
      legend.position='none',
      panel.background = element_blank()
    )
  
  return(p)
}



#############################################################
## Internal methods for enrichment analysis (not exported) ##
#############################################################

# This is a modified version of the PCGSE module
.pcgse = function(data, prcomp.output, feature.sets, pc.indexes, 
                  set.statistic, set.test) {
  
  # Sanity checks
  if (is.null(feature.sets))
    stop("'feature.sets' must be specified!")
  if (!(set.statistic %in% c("mean.diff", "rank.sum")))
    stop("set.statistic must be 'mean.diff' or 'rank.sum'")
  if (!(set.test %in% c("parametric", "cor.adj.parametric", "permutation")))
    stop("set.test must be one of 'parametric', 'cor.adj.parametric', 'permutation'")
  
  
  # Turn the feature set matrix into list form
  set.indexes <- feature.sets  
  if (is.matrix(feature.sets)) {
    set.indexes <- .createVarGroupList(var.groups=feature.sets)  
  }
  
  # Compute the feature statistics.
  feature.statistics <- matrix(0, nrow=ncol(data), ncol=length(pc.indexes))
  for (i in seq_along(pc.indexes)) {
    feature.statistics[,i] <- .compute_feature_statistics(
      data = data,
      prcomp.output = prcomp.output,
      pc.index = pc.indexes[i]
    )
  }
  
  # Compute the set statistics.
  if (set.test == "parametric" || set.test == "cor.adj.parametric") {
    if (set.statistic == "mean.diff") {
      results <- .pcgse_ttest(
        data = data, 
        prcomp.output = prcomp.output,
        pc.indexes = pc.indexes,
        set.indexes = set.indexes,
        feature.statistics = feature.statistics,
        cor.adjustment = (set.test == "cor.adj.parametric")
      )
    } else if (set.statistic == "rank.sum") {
      results <- .pcgse_wmw(
        data = data, 
        prcomp.output = prcomp.output,
        pc.indexes = pc.indexes,
        set.indexes = set.indexes,
        feature.statistics = feature.statistics,
        cor.adjustment = (set.test == "cor.adj.parametric")
      )
    }
  }
  
  # Add feature.statistics to the results
  results$feature.statistics <- feature.statistics
  
  return (results) 
}




# Turn the annotation matrix into a list of var group indexes for the valid sized var groups
.createVarGroupList <- function(var.groups) {
  var.group.indexes <- list()  
  for (i in seq_len(nrow(var.groups))) {
    member.indexes <- which(var.groups[i,]==1)
    var.group.indexes[[i]] <- member.indexes    
  }
  names(var.group.indexes) <- rownames(var.groups)    
  return (var.group.indexes)
}

# Computes the feature-level statistics
.compute_feature_statistics <- function(data, prcomp.output, pc.index) {
  feature.statistics <- prcomp.output$rotation[,pc.index]
  feature.statistics <- vapply(feature.statistics, abs, numeric(1))
  return (feature.statistics)
}

# Compute enrichment via t-test
#' @importFrom stats pt var
.pcgse_ttest <- function(data, prcomp.output, pc.indexes,
                         set.indexes, feature.statistics, cor.adjustment) {
  
  num.feature.sets <- length(set.indexes)
  
  # Create matrix for p-values
  p.values <- matrix(0, nrow=num.feature.sets, ncol=length(pc.indexes))  
  rownames(p.values) <- names(set.indexes)
  
  # Create matrix for set statistics
  set.statistics <- matrix(TRUE, nrow=num.feature.sets, ncol=length(pc.indexes))    
  rownames(set.statistics) <- names(set.indexes)    
  
  for (i in seq_len(num.feature.sets)) {
    indexes.for.feature.set <- set.indexes[[i]]
    m1 <- length(indexes.for.feature.set)
    not.set.indexes <- which(!(seq_len(ncol(data)) %in% indexes.for.feature.set))
    m2 <- length(not.set.indexes)
    
    if (cor.adjustment) {      
      # compute sample correlation matrix for members of feature set
      cor.mat <- cor(data[,indexes.for.feature.set], use = "complete.obs")
      # compute the mean pair-wise correlation 
      mean.cor <- (sum(cor.mat) - m1)/(m1*(m1-1))    
      # compute the VIF, using CAMERA formula from Wu et al., based on Barry et al.
      vif <- 1 + (m1 -1)*mean.cor
    }
    
    for (j in seq_along(pc.indexes)) {
      # get the feature statistics for this PC
      pc.feature.stats <- feature.statistics[,j]
      # compute the mean difference of the feature-level statistics
      mean.diff <- mean(pc.feature.stats[indexes.for.feature.set],na.rm=TRUE) - mean(pc.feature.stats[not.set.indexes], na.rm=TRUE)
      # compute the pooled standard deviation
      pooled.sd <- sqrt(((m1-1)*var(pc.feature.stats[indexes.for.feature.set], na.rm=TRUE) + (m2-1)*var(pc.feature.stats[not.set.indexes], na.rm=TRUE))/(m1+m2-2))
      # compute the t-statistic
      if (cor.adjustment) {
        t.stat <- mean.diff/(pooled.sd*sqrt(vif/m1 + 1/m2))
        df <- nrow(data)-2
      } else {
        t.stat <- mean.diff/(pooled.sd*sqrt(1/m1 + 1/m2))
        df <- m1+m2-2
      }
      set.statistics[i,j] <- t.stat      
      # compute the p-value via a two-sided test
      lower.p <- pt(t.stat, df=df, lower.tail=TRUE)
      upper.p <- pt(t.stat, df=df, lower.tail=FALSE)        
      p.values[i,j] <- 2*min(lower.p, upper.p)      
    }
  } 
  
  # Build the result list
  results <- list()
  results$p.values <- p.values
  results$statistics <- set.statistics  
  
  return (results)
}

# Compute enrichment via Wilcoxon Mann Whitney 
#' @importFrom stats wilcox.test pnorm
.pcgse_wmw <- function(data, prcomp.output, pc.indexes,
                       set.indexes, feature.statistics, cor.adjustment) {
  
  num.feature.sets <- length(set.indexes)
  
  # Create matrix for p-values
  p.values <- matrix(0, nrow=num.feature.sets, ncol=length(pc.indexes))  
  rownames(p.values) <- names(set.indexes)
  
  # Create matrix for set statistics
  set.statistics <- matrix(TRUE, nrow=num.feature.sets, ncol=length(pc.indexes))    
  rownames(set.statistics) <- names(set.indexes)    
  
  for (i in seq_len(num.feature.sets)) {
    indexes.for.feature.set <- set.indexes[[i]]
    m1 <- length(indexes.for.feature.set)
    not.set.indexes <- which(!(seq_len(ncol(data)) %in% indexes.for.feature.set))
    m2 <- length(not.set.indexes)
    
    if (cor.adjustment) {            
      # compute sample correlation matrix for members of feature set
      cor.mat <- cor(data[,indexes.for.feature.set], use="complete.obs")
      # compute the mean pair-wise correlation 
      mean.cor <- (sum(cor.mat) - m1)/(m1*(m1-1))    
    }
    
    for (j in seq_along(pc.indexes)) {
      # get the feature-level statistics for this PC
      pc.feature.stats <- feature.statistics[,j]
      # compute the rank sum statistic feature-level statistics
      wilcox.results <- wilcox.test(x=pc.feature.stats[indexes.for.feature.set],
                                    y=pc.feature.stats[not.set.indexes],
                                    alternative="two.sided", exact=FALSE, correct=FALSE)
      rank.sum = wilcox.results$statistic                
      if (cor.adjustment) {
        # Using correlation-adjusted formula from Wu et al.
        var.rank.sum <- ((m1*m2)/(2*pi))*
          (asin(1) + (m2 - 1)*asin(.5) + (m1-1)*(m2-1)*asin(mean.cor/2) +(m1-1)*asin((mean.cor+1)/2))
      } else {        
        var.rank.sum <- m1*m2*(m1+m2+1)/12
      }
      z.stat <- (rank.sum - (m1*m2)/2)/sqrt(var.rank.sum)
      set.statistics[i,j] <- z.stat
      
      # compute the p-value via a two-sided z-test
      lower.p <- pnorm(z.stat, lower.tail=TRUE)
      upper.p <- pnorm(z.stat, lower.tail=FALSE)        
      p.values[i,j] <- 2*min(lower.p, upper.p)
    }
  } 
  
  # Build the result list
  results <- list()
  results$p.values <- p.values
  results$statistics <- set.statistics  
  
  return (results)
}
