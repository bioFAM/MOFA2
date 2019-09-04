##########################################################
## Functions to perform Feature Set Enrichment Analysis ##
##########################################################

#' @title Run feature set Enrichment Analysis
#' @name run_enrichment 
#' @description Method to perform feature set enrichment analysis. Here we use a slightly modified version of the \link[PCGSE]{pcgse} function.
#' @param object a \code{\link{MOFA}} object.
#' @param view a character with the view name, or a numeric vector with the index of the view to use.
#' @param feature.sets data structure that holds feature set membership information. 
#' Must be a binary membership matrix (rows are feature sets and columns are features).
#' @param factors character vector with the factor names, or numeric vector with the index of the factors for which to perform the enrichment.
#' @param local.statistic the feature statistic used to quantify the association between each feature and each factor. Must be one of the following: loading (default), cor, z.
#' @param global.statistic the feature set statisic computed from the feature statistics. Must be one of the following: "mean.diff" (default) or "rank.sum".
#' @param statistical.test the statistical test used to compute the significance of the feature set statistics under a competitive null hypothesis.
#' Must be one of the following: "parametric" (default), "cor.adj.parametric", "permutation".
#' @param transformation optional transformation to apply to the feature-level statistics. Must be one of the following "none" or "abs.value" (default).
#' @param min.size Minimum size of a feature set (default is 10).
#' @param nperm number of permutations. Only relevant if statistical.test is set to "permutation". Default is 1000
#' @param cores number of cores to run the permutation analysis in parallel. Only relevant if statistical.test is set to "permutation". Default is 1
#' @param p.adj.method Method to adjust p-values factor-wise for multiple testing. Can be any method in p.adjust.methods(). Default uses Benjamini-Hochberg procedure.
#' @param alpha FDR threshold to generate lists of significant pathways. Default is 0.1
#' @details TO-DO
#' @return a list with three components: 
#' pval and pval.adj contain matrices with p-values and adjusted p-values, repectively. 
#' sigPathways contains a list with significant pathwayd at FDR alpha per factor.
#' @import foreach doParallel
#' @importFrom stats p.adjust var p.adjust.methods
#' @export

run_enrichment <- function(object, view, feature.sets, factors = "all",
                           local.statistic = c("loading", "cor", "z"), global.statistic = c("mean.diff", "rank.sum"),
                           statistical.test = c("parametric", "cor.adj.parametric", "permutation"), transformation = c("abs.value", "none"),
                           min.size = 10, nperm = 1000, cores = 1, p.adj.method = "BH", alpha = 0.1) {
  
  # Quality control
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (!(is(feature.sets, "matrix") & all(feature.sets %in% c(0,1)))) stop("feature.sets has to be a list or a binary matrix.")
  
  # Define views
  view  <- .check_and_get_views(object, view)
  
  # Define factors
  factors  <- .check_and_get_factors(object, factors)

  # Parse inputs
  local.statistic <- match.arg(local.statistic)
  transformation <- match.arg(transformation)
  global.statistic <- match.arg(global.statistic)
  statistical.test <- match.arg(statistical.test)
  
  # Collect observed data
  data <- get_data(object, views = view, as.data.frame = FALSE)[[1]]
  if(is(data, "list")) data <- Reduce(cbind, data) # concatenate groups
  data <- t(data)

  # Collect relevant expectations
  W <- get_weights(object, views=view, factors=factors)[[view]]
  Z <- get_factors(object, factors=factors)
  if(is(Z, "list")) Z <- Reduce(rbind, Z)
  stopifnot(rownames(data) == rownames(Z))
  
  # Check that there is no empty factor
  stopifnot( all(apply(Z,2,var, na.rm=TRUE)>0) )
  
  # Check if some features do not intersect between the feature sets and the observed data and remove them
  features <- intersect(colnames(data),colnames(feature.sets))
  if(length(features) == 0 ) stop("Feature names in feature.sets do not match feature names in model.")
  data <- data[,features]
  W <- W[features,]
  feature.sets <- feature.sets[,features]
  
  # Filter feature sets with small number of features
  feature.sets <- feature.sets[rowSums(feature.sets)>=min.size,]
    
  # Print options
  message("Running feature set Enrichment Analysis with the following options...")
  message(sprintf("View: %s", view))
  message(sprintf("Number of feature sets: %d", nrow(feature.sets)))
  message(sprintf("Local statistic: %s", local.statistic))
  message(sprintf("Transformation: %s", transformation))
  message(sprintf("Global statistic: %s", global.statistic))
  message(sprintf("Statistical test: %s", statistical.test))
  if (statistical.test=="permutation") {
    message(sprintf("Cores: %d", cores))
    message(sprintf("Number of permutations: %d", nperm))
  }
  
  # Non-parametric permutation test
  if (statistical.test == "permutation") {
    doParallel::registerDoParallel(cores=cores)
    
    null_dist_tmp <- foreach(rnd= seq_len(nperm)) %dopar% {
      perm <- sample(ncol(data))
      # Permute rows of the weight matrix to obtain a null distribution
      W_null <- W[perm,]
      rownames(W_null) <- rownames(W); colnames(W_null) <- colnames(W)
      # Permute columns of the data matrix correspondingly (only matters for cor.adjusted test)
      data_null <- data[,perm]
      rownames(data_null) <- rownames(data)
      
      # Compute null statistic
      s.null <- .pcgse(
        data=data_null, 
        prcomp.output = list(rotation=W_null, x=Z),
        pc.indexes = seq_along(factors), 
        feature.sets = feature.sets,
        feature.statistic = local.statistic,
        transformation = transformation,
        feature.set.statistic = global.statistic,
        feature.set.test = "parametric", nperm=NA)$statistic
      abs(s.null)
    }
    null_dist <- do.call("rbind", null_dist_tmp)
    colnames(null_dist) <- factors
    
    # Compute true statistics
    s.true <- .pcgse(
      data = data, 
      prcomp.output = list(rotation=W, x=Z),
      pc.indexes = seq_along(factors), 
      feature.sets = feature.sets,
      feature.statistic = local.statistic,
      transformation = transformation,
      feature.set.statistic = global.statistic,
      feature.set.test = "parametric", nperm=NA)$statistic
    colnames(s.true) <- factors
    rownames(s.true) <- rownames(feature.sets)
    
    # Calculate p-values based on fraction true statistic per factor and feature set is larger than permuted
    warning("A large number of permutations is required for the permutation approach!")
    xx <- array(unlist(null_dist_tmp),
                dim = c(nrow(null_dist_tmp[[1]]), ncol(null_dist_tmp[[1]]), length(null_dist_tmp)))
    ll <- lapply(seq_len(nperm), function(i) xx[,,i] > abs(s.true))
    values <- Reduce("+",ll)/nperm
    rownames(p.values) <- rownames(s.true); colnames(p.values) <- factors
    
    # Parametric test
  } else {
    results <- .pcgse(
      data = data,
      prcomp.output = list(rotation=W, x=Z),
      pc.indexes = seq_along(factors),
      feature.sets = feature.sets,
      feature.statistic = local.statistic,
      transformation = transformation,
      feature.set.statistic = global.statistic,
      feature.set.test = statistical.test, nperm=nperm)
  }
  
  # Parse results
  pathways <- rownames(feature.sets)
  colnames(results[["p.values"]]) <- colnames(results[["statistics"]]) <- colnames(results[["feature.statistics"]]) <- factors
  rownames(results[["p.values"]]) <- rownames(results[["statistics"]]) <- pathways
  rownames(results[["feature.statistics"]]) <- colnames(data)
  
  # adjust for multiple testing
  if(!p.adj.method %in%  p.adjust.methods) 
    stop("p.adj.method needs to be an element of p.adjust.methods")
  adj.p.values <- apply(results[["p.values"]], 2,function(lfw) p.adjust(lfw, method = p.adj.method))
  
  # obtain list of significant pathways
  sigPathways <- lapply(factors, function(j) rownames(adj.p.values)[adj.p.values[,j] <= alpha])
  
  # prepare output
  output <- list(
    pval = results[["p.values"]], 
    pval.adj = adj.p.values, 
    feature.statistics = results[["feature.statistics"]],
    set.statistics = results[["statistics"]],
    # pathway.activity = tmp,
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
#' @param object \code{\link{MOFAmodel}} object on which run_enrichment was performed
#' @param enrichment.results output of \link{runEnrichmentAnalysis} function
#' @param factor a string with the factor name or an integer with the factor index
#' @param alpha p.value threshold to filter out gene sets
#' @param max.pathways maximum number of enriched pathways to display
#' @param adjust use adjusted p-values?
#' @details it requires \code{\link{run_enrichment}} to be run beforehand.
#' @return a \code{ggplot2} object
#' @import ggplot2
#' @importFrom utils head
#' @export
plot_enrichment <- function(enrichment.results, factor, alpha = 0.1, max.pathways = 25, adjust = TRUE, 
                           text_size = 1.0, dot_size = 5.0) {
  
  # Sanity checks
  stopifnot(is.numeric(alpha)) 
  stopifnot(length(factor)==1) 
  if (is.numeric(factor)) factor <- colnames(enrichment.results$pval.adj)[factor]
  if(!factor %in% colnames(enrichment.results$pval)) 
    stop(paste0("No gene set enrichment calculated for factor ", factor))
  
  # get p-values
  if(adjust) p.values <- enrichment.results$pval.adj else p.values <- enrichment.results$pval
  
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
  tmp$logp <- -log10(tmp$pvalue)
  
  # Annotate significcant pathways
  # tmp$sig <- factor(tmp$pvalue<alpha)
  
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
#' @param enrichment.results output of \link{runEnrichmentAnalysis} function
#' @param alpha FDR threshold to filter out unsignificant feature sets which are
#'  not represented in the heatmap. Default is 0.05.
#' @param log_scale boolean indicating whether to plot the log of the p.values.
#' @param ... extra arguments to be passed to \link{pheatmap} function
#' @details it requires \code{\link{run_enrichment}} to be run beforehand.
#' @return produces a heatmap
#' @import pheatmap
#' @importFrom grDevices colorRampPalette
#' @export
plot_enrichment_heatmap <- function(enrichment.results, alpha = 0.05, log_scale = TRUE, ...) {
  
  # get p-values
  p.values <- enrichment.results$pval.adj
  p.values <- p.values[!apply(p.values, 1, function(x) sum(x>=alpha)) == ncol(p.values),, drop=FALSE]
  
  # Apply Log transform
  if (log_scale) {
    p.values <- -log10(p.values)
    alpha <- -log10(alpha)
    col <- colorRampPalette(c("lightgrey", "red"))(n=10)
  } else {
    col <- colorRampPalette(c("red", "lightgrey"))(n=10)
  }
  
  # Generate heatmap
  pheatmap::pheatmap(p.values, color = col, ...)
}


#' @title Plot detailed output of the Feature Set Enrichment Analysis
#' @name plotEnrichmentDetailed
#' @description Method to plot a detailed output of the Feature Set Enrichment Analyisis (FSEA). \cr
#' Each row corresponds to a significant pathway, sorted by statistical significance, and each dot corresponds to a gene. \cr
#' For each pathway, we display the top genes of the pathway sorted by the corresponding feature statistic (by default, the absolute value of the loading)
#' The top genes with the highest statistic (max.genes argument) are displayed and labeled in black. The remaining genes are colored in grey.
#' @param object \code{\link{MOFAmodel}} object on which FSEA was performed
#' @param factor string with factor name or numeric with factor index
#' @param enrichment.results output of \link{runEnrichmentAnalysis} function
#' @param feature.sets data structure that holds feature set membership information, as used in the \link{runEnrichmentAnalysis} function.
#' @param alpha p.value threshold to filter out feature sets
#' @param max.pathways maximum number of enriched pathways to display
#' @param max.genes maximum number of genes to display, per pathway
#' @param text_size size of the text to label the top genes
#' @param adjust use adjusted p-values?
#' @return a \code{ggplot2} object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @importFrom dplyr top_n
#' @importFrom ggrepel geom_text_repel
#' @export
plotEnrichmentDetailed <- function(object, factor, feature.sets, enrichment.results, 
                                   adjust = TRUE, alpha = 0.1, max.genes = 5, max.pathways = 10, text_size = 3) {
  
  # Sanity checks
  stopifnot(length(factor)==1) 
  if(is.numeric(factor)) factor <- factorNames(object)[factor]
  if(!factor %in% colnames(enrichment.results$pval)) 
    stop(paste0("No feature set enrichment calculated for factor ", factor, ".\n Use runEnrichmentAnalysis first."))
  
  # Fetch and prepare data  
  
  # foo
  foo <- melt(enrichment.results$feature.statistics, factorsAsStrings=TRUE)
  colnames(foo) <- c("feature","factor","feature.statistic")
  foo <- foo[foo$factor==factor,]
  foo$feature <- as.character(foo$feature)
  
  # bar
  bar <- melt(feature.sets)
  bar <- bar[bar$value==1,c(1,2)]
  colnames(bar) <- c("pathway","feature")
  bar$pathway <- as.character(bar$pathway)
  bar$feature <- as.character(bar$feature)
  
  # baz
  if (adjust) {
    baz <- melt(enrichment.results$pval.adj)
  } else {
    baz <- melt(enrichment.results$pval)
  }
  colnames(baz) <- c("pathway","factor","pvalue")
  baz$pathway <- as.character(baz$pathway)
  baz <- baz[baz$factor==factor,]
  
  # Filter out pathways by p-values
  baz <- baz[baz$pvalue<=alpha,,drop=FALSE]
  if(nrow(baz)==0) {
    stop("No siginificant pathways at the specified alpha threshold. \n
         For an overview use plotEnrichmentHeatmap() or plotEnrichmentBars().")
  }
  if (nrow(baz) > max.pathways)
    baz <- head(baz[order(baz$pvalue),],n=max.pathways)
  
  # order pathways according to significance
  baz$pathway <- factor(baz$pathway, levels = baz$pathway[order(baz$pvalue, decreasing = TRUE)])
  
  # Merge
  foobar <- merge(foo, bar, by="feature")
  tmp <- merge(foobar, baz, by=c("pathway","factor"))
  
  # Select the top N features with the largest feature.statistic (per pathway)
  tmp_filt <- top_n(group_by(tmp, pathway), n=max.genes, abs(feature.statistic))
  
  # Add number of features and p-value per pathway
  pathways <- unique(tmp_filt$pathway)
  
  # Add Ngenes and p-values to the pathway name
  df <- data.frame(pathway=pathways, nfeatures=rowSums(feature.sets)[pathways])
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
    labs(x="", y="Feature statistic", title="") +
    coord_flip() +
    theme(
      axis.line = element_line(color="black"),
      axis.text.y = element_text(size=rel(1.2), hjust=1, color='black'),
      axis.text.x = element_text(size=rel(1.2), vjust=0.5, color='black'),
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
.pcgse = function(data, prcomp.output, pc.indexes=1, feature.sets, feature.statistic="z", transformation="none", 
                  feature.set.statistic="mean.diff", feature.set.test="cor.adj.parametric", nperm=9999) {
  
  current.warn = getOption("warn")
  options(warn=-1)
  
  # Sanity checks
  if (is.null(feature.sets)) {
    stop("'feature.sets' must be specified!")
  }   
  options(warn=current.warn) 
  if (!(feature.statistic %in% c("loading", "cor", "z"))) {
    stop("feature.statistic must be 'loading', 'cor' or 'z'")
  }  
  if (!(transformation %in% c("none", "abs.value"))) {
    stop("transformation must be 'none' or 'abs.value'")
  }  
  if (!(feature.set.statistic %in% c("mean.diff", "rank.sum"))) {
    stop("feature.set.statistic must be 'mean.diff' or 'rank.sum'")
  }    
  if (!(feature.set.test %in% c("parametric", "cor.adj.parametric", "permutation"))) {
    stop("feature.set.test must be one of 'parametric', 'cor.adj.parametric', 'permutation'")
  }
  if (feature.set.test == "permutation" & feature.statistic == "loading") { 
    stop("feature.statistic cannot be set to 'loading' if feature.set.test is 'permutation'")
  }
  if (!is.matrix(feature.sets) & feature.set.test == "permutation") {
    stop("feature.sets must be specified as a binary membership matrix if
         feature.set.test is set to 'permutation'") 
  }  
  # if (feature.set.test == "parametric") {
  # warning("The 'parametric' test option ignores the correlation between feature-level
  #         test statistics and therefore has an inflated type I error rate.\n ",
  #         "This option should only be used for evaluation purposes.")
  # }  
  # if (feature.set.test == "permutation") {
  # warning("The 'permutation' test option can be extremely computationally expensive
  #         given the required modifications to the safe() function. ",
  #         "For most applications, it is recommended that feature.set.test
  #         is set to 'cor.adj.parametric'.")
  # }
  
  # Turn the feature set matrix into list form if feature.set.test is not "permutation"
  feature.set.indexes = feature.sets  
  if (is.matrix(feature.sets)) {
    feature.set.indexes = .createVarGroupList(var.groups=feature.sets)  
  }
  
  n = nrow(data)
  p = ncol(data)
  
  # Compute the feature statistics.
  feature.statistics = matrix(0, nrow=p, ncol=length(pc.indexes))
  for (i in seq_along(pc.indexes)) {
    pc.index = pc.indexes[i]
    feature.statistics[,i] = .computefeatureStatistics(
      data = data,
      prcomp.output = prcomp.output,
      pc.index = pc.index,
      feature.statistic = feature.statistic,
      transformation = transformation
    )
  }
  
  # Compute the feature-set statistics.
  if (feature.set.test == "parametric" | feature.set.test == "cor.adj.parametric") {
    if (feature.set.statistic == "mean.diff") {
      results = .pcgseViaTTest(
        data = data, 
        prcomp.output = prcomp.output,
        pc.indexes = pc.indexes,
        feature.set.indexes = feature.set.indexes,
        feature.statistics = feature.statistics,
        cor.adjustment = (feature.set.test == "cor.adj.parametric")
      )
    } else if (feature.set.statistic == "rank.sum") {
      results = .pcgseViaWMW(
        data = data, 
        prcomp.output = prcomp.output,
        pc.indexes = pc.indexes,
        feature.set.indexes = feature.set.indexes,
        feature.statistics = feature.statistics,
        cor.adjustment = (feature.set.test == "cor.adj.parametric")
      )
    }
  }
  
  # Add feature.statistics to the results
  results[["feature.statistics"]] <- feature.statistics
  
  return (results) 
}




# Turn the annotation matrix into a list of var group indexes for the valid sized var groups
.createVarGroupList = function(var.groups) {
  var.group.indexes = list()  
  for (i in seq_len(nrow(var.groups))) {
    member.indexes = which(var.groups[i,]==1)
    var.group.indexes[[i]] = member.indexes    
  }
  names(var.group.indexes) = rownames(var.groups)    
  return (var.group.indexes)
}

# Computes the feature-level statistics
.computefeatureStatistics = function(data, prcomp.output, pc.index, feature.statistic, transformation) {
  p = ncol(data)
  n = nrow(data)
  feature.statistics = rep(0, p)
  if (feature.statistic == "loading") {
    # get the PC loadings for the selected PCs
    feature.statistics = prcomp.output$rotation[,pc.index]
  } else {
    # compute the Pearson correlation between the selected PCs and the data
    feature.statistics = cor(data, prcomp.output$x[,pc.index], use = "complete.obs") 
    if (feature.statistic == "z") {
      # use Fisher's Z transformation to convert to Z-statisics
      feature.statistics = vapply(feature.statistics, function(x) {
        return (sqrt(n-3)*atanh(x))}, numeric(1))      
    }    
  }
  
  # Absolute value transformation of the feature-level statistics if requested
  if (transformation == "abs.value") {
    feature.statistics = vapply(feature.statistics, abs, numeric(1))
  }  
  
  return (feature.statistics)
}

# Compute enrichment via t-test
#' @importFrom stats pt var
.pcgseViaTTest = function(data, prcomp.output, pc.indexes,
                          feature.set.indexes, feature.statistics, cor.adjustment) {
  
  num.feature.sets = length(feature.set.indexes)
  n= nrow(data)
  p.values = matrix(0, nrow=num.feature.sets, ncol=length(pc.indexes))  
  rownames(p.values) = names(feature.set.indexes)
  feature.set.statistics = matrix(TRUE, nrow=num.feature.sets, ncol=length(pc.indexes))    
  rownames(feature.set.statistics) = names(feature.set.indexes)    
  
  for (i in seq_len(num.feature.sets)) {
    indexes.for.feature.set = feature.set.indexes[[i]]
    m1 = length(indexes.for.feature.set)
    not.feature.set.indexes = which(!(seq_len(ncol(data)) %in% indexes.for.feature.set))
    m2 = length(not.feature.set.indexes)
    
    if (cor.adjustment) {      
      # compute sample correlation matrix for members of feature set
      cor.mat = cor(data[,indexes.for.feature.set], use = "complete.obs")
      # compute the mean pair-wise correlation 
      mean.cor = (sum(cor.mat) - m1)/(m1*(m1-1))    
      # compute the VIF, using CAMERA formula from Wu et al., based on Barry et al.
      vif = 1 + (m1 -1)*mean.cor
    }
    
    for (j in seq_along(pc.indexes)) {
      # get the feature-level statistics for this PC
      pc.feature.stats = feature.statistics[,j]
      # compute the mean difference of the feature-level statistics
      mean.diff = mean(pc.feature.stats[indexes.for.feature.set]) -
        mean(pc.feature.stats[not.feature.set.indexes])
      # compute the pooled standard deviation
      pooled.sd = sqrt(((m1-1)*var(pc.feature.stats[indexes.for.feature.set]) +
                          (m2-1)*var(pc.feature.stats[not.feature.set.indexes]))/(m1+m2-2))      
      # compute the t-statistic
      if (cor.adjustment) {
        t.stat = mean.diff/(pooled.sd*sqrt(vif/m1 + 1/m2))
        df = n-2
      } else {
        t.stat = mean.diff/(pooled.sd*sqrt(1/m1 + 1/m2))
        df = m1+m2-2
      }
      feature.set.statistics[i,j] = t.stat      
      # compute the p-value via a two-sided test
      lower.p = pt(t.stat, df=df, lower.tail=TRUE)
      upper.p = pt(t.stat, df=df, lower.tail=FALSE)        
      p.values[i,j] = 2*min(lower.p, upper.p)      
    }
  } 
  
  # Build the result list
  results = list()
  results$p.values = p.values
  results$statistics = feature.set.statistics  
  
  return (results)
}

# Compute enrichment via Wilcoxon Mann Whitney 
#' @importFrom stats wilcox.test pnorm
.pcgseViaWMW = function(data, prcomp.output, pc.indexes,
                        feature.set.indexes, feature.statistics, cor.adjustment) {
  
  num.feature.sets = length(feature.set.indexes)
  n= nrow(data)
  p.values = matrix(0, nrow=num.feature.sets, ncol=length(pc.indexes))  
  rownames(p.values) = names(feature.set.indexes)
  feature.set.statistics = matrix(TRUE, nrow=num.feature.sets, ncol=length(pc.indexes))    
  rownames(feature.set.statistics) = names(feature.set.indexes)    
  
  for (i in seq_len(num.feature.sets)) {
    indexes.for.feature.set = feature.set.indexes[[i]]
    m1 = length(indexes.for.feature.set)
    not.feature.set.indexes = which(!(seq_len(ncol(data)) %in% indexes.for.feature.set))
    m2 = length(not.feature.set.indexes)
    
    if (cor.adjustment) {            
      # compute sample correlation matrix for members of feature set
      cor.mat = cor(data[,indexes.for.feature.set])
      # compute the mean pair-wise correlation 
      mean.cor = (sum(cor.mat) - m1)/(m1*(m1-1))    
    }
    
    for (j in seq_along(pc.indexes)) {
      # get the feature-level statistics for this PC
      pc.feature.stats = feature.statistics[,j]
      # compute the rank sum statistic feature-level statistics
      wilcox.results = wilcox.test(x=pc.feature.stats[indexes.for.feature.set],
                                   y=pc.feature.stats[not.feature.set.indexes],
                                   alternative="two.sided", exact=FALSE, correct=FALSE)
      rank.sum = wilcox.results$statistic                
      if (cor.adjustment) {
        # Using correlation-adjusted formula from Wu et al.
        var.rank.sum = ((m1*m2)/(2*pi))*
          (asin(1) + (m2 - 1)*asin(.5) + (m1-1)*(m2-1)*asin(mean.cor/2) +(m1-1)*asin((mean.cor+1)/2))
      } else {        
        var.rank.sum = m1*m2*(m1+m2+1)/12
      }
      z.stat = (rank.sum - (m1*m2)/2)/sqrt(var.rank.sum)
      feature.set.statistics[i,j] = z.stat
      
      # compute the p-value via a two-sided z-test
      lower.p = pnorm(z.stat, lower.tail=TRUE)
      upper.p = pnorm(z.stat, lower.tail=FALSE)        
      p.values[i,j] = 2*min(lower.p, upper.p)
    }
  } 
  
  # Build the result list
  results = list()
  results$p.values = p.values
  results$statistics = feature.set.statistics  
  
  return (results)
}
