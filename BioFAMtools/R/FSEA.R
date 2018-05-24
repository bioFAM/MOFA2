##########################################################
## Functions to perform Feature Set Enrichment Analysis ##
##########################################################

#' @title Feature Set Enrichment Analysis
#' @name feature_set_enrichment_analysis 
#' @description Method to perform feature set enrichment analysis. Here we use a slightly modified version of the \link[PCGSE]{pcgse} function.
#' @param object a \code{\link{BioFAModel}} object.
#' @param view name of the view
#' @param feature.sets data structure that holds feature set membership information. Must be either a binary membership matrix (rows are feature sets and columns are features) or a list of feature set indexes (see vignette for details).
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
#' @return a list with three components: pval and pval.adj contain matrices with p-values and adjusted p-values, repectively. sigPathways contains a list with significant pathwayd at FDR alpha per factor.
#' @import foreach doParallel
#' @importFrom stats p.adjust
#' @export

FSEA <- function(object, view, feature.sets, factors = "all", local.statistic = c("loading", "cor", "z"),
                                         global.statistic = c("mean.diff", "rank.sum"), statistical.test = c("parametric", "cor.adj.parametric", "permutation"),
                                         transformation = c("abs.value", "none"), min.size = 10, nperm = 1000, cores = 1, p.adj.method = "BH", alpha=0.1) {
  
  # Parse inputs
  local.statistic <- match.arg(local.statistic)
  transformation <- match.arg(transformation)
  global.statistic <- match.arg(global.statistic)
  statistical.test <- match.arg(statistical.test)


  # Define factors
  if (paste0(factors,collapse="") == "all") { factors <- factors_names(object) } 
    else if(is.numeric(factors)) {
      if (object@model_options$learn_intercept) factors <- factors_names(object)[factors+1]
      else factors <- factors_names(object)[factors]
    }
      else{ stopifnot(all(factors %in% factors_names(object))) }

  # remove intercept factors
  factors <- factors[factors!="intercept"]
  
  # Collect observed data
  data <- object@training_data[[view]]
  if(class(data)=="list") data <- Reduce(cbind, data)
  data <- t(data)

  # Collect relevant expectations
  W <- get_weights(object, views=view, factors=factors)[[view]]
  Z <- get_factors(object, factors=factors)
  if(class(Z)=="list") Z <- Reduce(rbind, Z)
  stopifnot(rownames(data) == rownames(Z))
  
  # Check that there is no constant factor
  stopifnot( all(apply(Z,2,var, na.rm=T)>0) )
    
  # turn feature.sets into binary membership matrices if provided as list
  if(class(feature.sets) == "list") {
    features <- Reduce(union, feature.sets)
    feature.sets <- sapply(names(feature.sets), function(nm) {
      tmp <- features %in% feature.sets[[nm]]
      names(tmp) <- features
      tmp
    })
    feature.sets <-t(feature.sets)*1
  }

  if(!(class(feature.sets)=="matrix" & all(feature.sets %in% c(0,1)))) stop("feature.sets has to be a list or a binary matrix.")
  
  # Check if some features do not intersect between the feature sets and the observed data and remove them
  features <- intersect(colnames(data),colnames(feature.sets))
  if(length(features) == 0 ) stop("Feautre names in feature.sets do not match feature names in model.")
  data <- data[,features]
  W <- W[features,]
  feature.sets <- feature.sets[,features]
  
  # Filter feature sets with small number of features
  feature.sets <- feature.sets[rowSums(feature.sets)>=min.size,]
    
  # Print options
  # message("Doing feature Ontology Enrichment Analysis with the following options...")
  # message(sprintf("View: %s", view))
  # message(sprintf("Latent variables: %s", paste(as.character(factors),collapse=" ")))
  # message(sprintf("Number of feature sets: %d", nrow(feature.sets)))
  # message(sprintf("Local statistic: %s", local.statistic))
  # message(sprintf("Transformation: %s", transformation))
  # message(sprintf("Global statistic: %s", global.statistic))
  # message(sprintf("Statistical test: %s", statistical.test))
  # if (statistical.test=="permutation") {
  #   message(sprintf("Cores: %d", cores))
  #   message(sprintf("Number of permutations: %d", nperm))
  # }
  
  # use own version for permutation test because of bugs in PCGSE package
  if (statistical.test == "permutation") {
    doParallel::registerDoParallel(cores=cores)
    	
    null_dist_tmp <- foreach(rnd=1:nperm) %dopar% {
      perm <- sample(ncol(data))
      # Permute rows of the weight matrix to obtain a null distribution
      W_null <- W[perm,]
      rownames(W_null) <- rownames(W); colnames(W_null) <- colnames(W)
      # Permute columns of the data matrix correspondingly (only matters for cor.adjusted test)
	  data_null <- data[,perm]
      rownames(data_null) <- rownames(data)
      
      # Compute null statistic
      s.null <- pcgse(data=data_null, prcomp.output=list(rotation=W_null, x=Z), pc.indexes=1:length(factors), feature.sets=feature.sets, feature.statistic=local.statistic,
                      transformation=transformation, feature.set.statistic=global.statistic, feature.set.test="parametric", nperm=NA)$statistic
      abs(s.null)
    }
    null_dist <- do.call("rbind", null_dist_tmp)
    colnames(null_dist) <- factors
    
    # Compute true statistics
    s.true <- pcgse(data=data, prcomp.output=list(rotation=W, x=Z), pc.indexes=1:length(factors), feature.sets=feature.sets, feature.statistic=local.statistic,
                    transformation=transformation, feature.set.statistic=global.statistic, feature.set.test="parametric", nperm=NA)$statistic
    colnames(s.true) <- factors
    rownames(s.true) <- rownames(feature.sets)
    
    # Calculate p-values based on fraction true statistic per factor and gene set is larger than permuted
    warning("A large number of permutations is required for the permutation approach!")
    xx <- array(unlist(null_dist_tmp), dim = c(nrow(null_dist_tmp[[1]]), ncol(null_dist_tmp[[1]]), length(null_dist_tmp)))
	ll <- lapply(1:nperm, function(i) xx[,,i] > abs(s.true))
	p.values <- Reduce("+",ll)/nperm
	rownames(p.values) <- rownames(s.true); colnames(p.values) <- factors

# parametric version
  } else {
    p.values <- pcgse(data=data, prcomp.output=list(rotation=W, x=Z), pc.indexes=1:length(factors), feature.sets=feature.sets, feature.statistic=local.statistic,
                      transformation=transformation, feature.set.statistic=global.statistic, feature.set.test=statistical.test, nperm=nperm)$p.values
    colnames(p.values) <- factors
    rownames(p.values) <- rownames(feature.sets)
    }
  
  # adjust for multiple testing per factor
  if(!p.adj.method %in%  p.adjust.methods) stop("p.adj.method needs to be an element of p.adjust.methods")
  adj.p.values <- apply(p.values, 2,function(lfw) p.adjust(lfw, method = p.adj.method))

  # list of significant pathwasy at level alpha
  sigPathways <- lapply(factors, function(j) rownames(adj.p.values)[adj.p.values[,j] <= alpha])

  return(list(pval = p.values, pval.adj = adj.p.values, sigPathways=sigPathways))
}


########################
## Plotting functions ##
########################


#' @title Line plot of Feature Set Enrichment Analysis results
#' @name lineplot_FSEA
#' @description Line plot of the Feature Set Enrichment Analyisis results for a specific latent variable
#' @param fsea.out output of \link{FeatureSetEnrichmentAnalysis} function
#' @param factor Factor for which to show enriched pathways in the lineplot
#' @param threshold p.value threshold to filter out feature sets
#' @param max.pathways maximum number of enriched pathways to display
#' @param adjust use multiple testing correction
#' @details fill this
#' @return nothing
#' @import ggplot2
#' @export
lineplot_FSEA <- function(fsea.out, factor, threshold=0.1, max.pathways=25, adjust=T) {
  
  # Sanity checks
  # (...)
  
  # get p-values
  if(adjust) p.values <- fsea.out$pval.adj else p.values <- fsea.out$pval

  # Get data  
  tmp <- as.data.frame(p.values[,factor, drop=F])
  tmp$pathway <- rownames(tmp)
  colnames(tmp) <- c("pvalue")
  
  # Filter out pathways
  tmp <- tmp[tmp$pvalue<=threshold,,drop=F]
  if(nrow(tmp)==0) {
    warning("No siginificant pathways at the specified threshold. For an overview use Heatmap_FeatureSetEnrichmentAnalysis().")
    return()
  }
  
  # If there are too many pathways enriched, just keep the 'max_pathways' more significant
  if (nrow(tmp) > max.pathways)
    tmp <- head(tmp[order(tmp$pvalue),],n=max.pathways)
  
  # Convert pvalues to log scale (add a small regulariser to avoid numerical errors)
  tmp$log <- -log10(tmp$pvalue)
  
  # Annotate significcant pathways
  # tmp$sig <- factor(tmp$pvalue<threshold)
  
  #order according to significance
  tmp$pathway <- factor(tmp$pathway <- rownames(tmp), levels = tmp$pathway[order(tmp$pvalue, decreasing = T)])
  
  p <- ggplot(tmp, aes(x=pathway, y=log)) +
    # ggtitle(paste("Enriched sets in factor", factor)) +
    geom_point(size=5) +
    geom_hline(yintercept=-log10(threshold), linetype="longdash") +
    # scale_y_continuous(limits=c(0,7)) +
    scale_color_manual(values=c("black","red")) +
    geom_segment(aes(xend=pathway, yend=0)) +
    ylab("-log pvalue") +
    coord_flip() +
    theme(
      axis.text.y = element_text(size=rel(1.2), hjust=1, color='black'),
      axis.text.x = element_text(size=rel(1.2), vjust=0.5, color='black'),
      axis.title.y=element_blank(),
      legend.position='none',
      panel.background = element_blank()
    )
  
  return(p)
}

#' @title Heatmap of Feature Set Enrichment Analysis results
#' @name heatmap_FSEA
#' @description this method generates a heatmap with the adjusted p.values that result from the the feature set enrichment analysis. Rows are feature sets and columns are factors.
#' @param fsea.out output of \link{FeatureSetEnrichmentAnalysis} function
#' @param factors : to plot only a subset of given factors
#' @param threshold p.value threshold to filter out unsignificant feature sets. Default is 0.05.
#' @param log boolean indicating whether to plot the log of the p.values.
#' @param ... extra arguments to be passed to \link{pheatmap} function
#' @details fill this
#' @import pheatmap
#' @importFrom grDevices colorRampPalette
#' @export

heatmap_FSEA <- function(fsea.out, factors="all", threshold = 0.05, log = TRUE, ...) {

  # get p-values
  p.values <- fsea.out$pval.adj
  if (factors!="all"){
    p.values <- p.values[,factors]
  }
  p.values <- p.values[!apply(p.values, 1, function(x) sum(x>=threshold)) == ncol(p.values),]

  # Apply Log transform
  if (log==T) {
    p.values <- -log10(p.values)
    threshold <- -log10(threshold)
    col <- colorRampPalette(c("lightgrey", "red"))(n=10)
  } else {
    col <- colorRampPalette(c("red", "lightgrey"))(n=10)
  }
  
  # Generate heatmap
  pheatmap::pheatmap(p.values, color = col, ...)
}


#' @title Barplot of Feature Set Enrichment Analysis results
#' @name barplot_FSEA
#' @description this method generates a barplot with the number of enriched feature sets per factor
#' @param fsea.out output of \link{FeatureSetEnrichmentAnalysis} function
#' @param alpha FDR threshold for calling enriched feature sets. Default is 0.05
#' @details TO-DO: IT NEEDS A BIT MORE WORK ON THE THEME
#' @return \link{ggplot} object
#' @import ggplot2
#' @importFrom grDevices colorRampPalette
#' @export
barplot_FSEA <- function(fsea.out, alpha = 0.05) {

  # Get enriched pathways at FDR of alpha
  pathwayList <- apply(fsea.out$pval.adj, 2, function(f) names(f)[f<=alpha])
  pathwaysDF <- reshape2::melt(pathwayList, value.name="pathway")
  colnames(pathwaysDF) <- c("pathway", "factor")
  
  # Count enriched gene sets per pathway
  # ERROR. no visible binding for global variable ‘n_enriched’
  pathwaysSummary <- dplyr::group_by(pathwaysDF,factor)
  pathwaysSummary <- dplyr::summarise(pathwaysSummary, n_enriched=length(pathway)) 
  pathwaysSummary <- rbind(pathwaysSummary, data.frame(factor = colnames(fsea.out$pval)[!colnames(fsea.out$pval) %in% pathwaysSummary$factor],
                                                     n_enriched=0))
  
  # Generate plot (TO-DO: IT NEEDS A BIT MORE THEME)
  ggplot(pathwaysSummary, aes(x=factor, y=n_enriched)) +
    geom_bar(stat="identity") + coord_flip() + 
    ylab(paste0("Enriched gene sets at FDR", alpha*100,"%")) +
    xlab("Factor") + 
    theme(
      legend.position = "bottom"
    )
}



##############################################
## From here downwards are internal methods ##
##############################################

# This is a modified version of the PCGSE module
pcgse = function(data, 
                 prcomp.output, 
                 pc.indexes=1, 
                 feature.sets,
                 feature.statistic="z",
                 transformation="none",
                 feature.set.statistic="mean.diff",
                 feature.set.test="cor.adj.parametric",
                 nperm=9999 # for feature.set.test value of "permutation"
) {
  current.warn = getOption("warn")
  options(warn=-1)
  # if (is.na(data)) {
  #   stop("'data must' be specified!")
  # }  
  if (is.na(feature.sets)) {
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
    stop("feature.sets must be specified as a binary membership matrix if feature.set.test is set to 'permutation'") 
  }  
  # if (feature.set.test == "parametric") {
  #   warning("The 'parametric' test option ignores the correlation between feature-level test statistics and therefore has an inflated type I error rate. ",
  #           "This option should only be used for evaluation purposes.")    
  # }  
  # if (feature.set.test == "permutation") {
  #   warning("The 'permutation' test option can be extremely computationally expensive given the required modifications to the safe() function. ",
  #           "For most applications, it is recommended that feature.set.test is set to 'cor.adj.parametric'.")
  # }
  
  # Turn the feature set matrix into list form if feature.set.test is not "permutation"
  feature.set.indexes = feature.sets  
  if (is.matrix(feature.sets)) {
    feature.set.indexes = create_var_group_list(var.groups=feature.sets)  
  }
  
  n = nrow(data)
  p = ncol(data)
  
  # Compute the feature-level statistics.
  feature.statistics = matrix(0, nrow=p, ncol=length(pc.indexes))
  for (i in 1:length(pc.indexes)) {
    pc.index = pc.indexes[i]
    feature.statistics[,i] = compute_feature_statistics(data=data, prcomp.output=prcomp.output, pc.index=pc.index, feature.statistic, transformation)
  }
  
  # Perform the specified feature set test for each feature set on each specified PC using the feature-level statistics
  if (feature.set.test == "parametric" | feature.set.test == "cor.adj.parametric") {
    if (feature.set.statistic == "mean.diff") {
      results = pcgse_via_ttest(data=data, prcomp.output=prcomp.output, pc.indexes=pc.indexes, feature.set.indexes=feature.set.indexes,
                              feature.statistics=feature.statistics, cor.adjustment=(feature.set.test == "cor.adj.parametric"))      
    } else if (feature.set.statistic == "rank.sum") {
      results = pcgse_via_WMW(data=data, prcomp.output=prcomp.output, pc.indexes=pc.indexes, feature.set.indexes=feature.set.indexes,
                            feature.statistics=feature.statistics, cor.adjustment=(feature.set.test == "cor.adj.parametric"))
    }     
  } else if (feature.set.test == "permutation") {
    # results = pcgseViaSAFE(data=data, prcomp.output=prcomp.output, pc.indexes=pc.indexes, feature.set.indexes=feature.set.indexes, 
    #                        feature.statistic=feature.statistic, transformation=transformation, feature.set.statistic=feature.set.statistic, nperm=nperm)
    results = pcgseViaPermutation(data=data, prcomp.output=prcomp.output, pc.indexes=pc.indexes, feature.set.indexes=feature.set.indexes, 
                                  feature.statistics=feature.statistics, feature.set.statistic=feature.set.statistic, nperm=nperm)        
  }
  
  return (results) 
}




# Turn the annotation matrix into a list of var group indexes for the valid sized var groups
create_var_group_list = function(var.groups) {
  var.group.indexes = list()  
  for (i in 1:nrow(var.groups)) {
    member.indexes = which(var.groups[i,]==1)
    var.group.indexes[[i]] = member.indexes    
  }
  names(var.group.indexes) = rownames(var.groups)    
  return (var.group.indexes)
}

# Computes the feature-level statistics
compute_feature_statistics = function(data, prcomp.output, pc.index, feature.statistic, transformation) {
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
      feature.statistics = sapply(feature.statistics, function(x) {
        return (sqrt(n-3)*atanh(x))})      
    }    
  }
  
  # Absolute value transformation of the feature-level statistics if requested
  if (transformation == "abs.value") {
    feature.statistics = sapply(feature.statistics, abs)
  }  
  
  return (feature.statistics)
}

# Compute enrichment via t-test
pcgse_via_ttest = function(data, prcomp.output, pc.indexes, feature.set.indexes, feature.statistics, cor.adjustment) {
  
  num.feature.sets = length(feature.set.indexes)
  n= nrow(data)
  p.values = matrix(0, nrow=num.feature.sets, ncol=length(pc.indexes))  
  rownames(p.values) = names(feature.set.indexes)
  feature.set.statistics = matrix(T, nrow=num.feature.sets, ncol=length(pc.indexes))    
  rownames(feature.set.statistics) = names(feature.set.indexes)    
  
  for (i in 1:num.feature.sets) {
    indexes.for.feature.set = feature.set.indexes[[i]]
    m1 = length(indexes.for.feature.set)
    not.feature.set.indexes = which(!(1:ncol(data) %in% indexes.for.feature.set))
    m2 = length(not.feature.set.indexes)
    
    if (cor.adjustment) {      
      # compute sample correlation matrix for members of feature set
      cor.mat = cor(data[,indexes.for.feature.set], use = "complete.obs")
      # compute the mean pair-wise correlation 
      mean.cor = (sum(cor.mat) - m1)/(m1*(m1-1))    
      # compute the VIF, using CAMERA formula from Wu et al., based on Barry et al.
      vif = 1 + (m1 -1)*mean.cor
    }
    
    for (j in 1:length(pc.indexes)) {
      # get the feature-level statistics for this PC
      pc.feature.stats = feature.statistics[,j]
      # compute the mean difference of the feature-level statistics
      mean.diff = mean(pc.feature.stats[indexes.for.feature.set]) - mean(pc.feature.stats[not.feature.set.indexes])
      # compute the pooled standard deviation
      pooled.sd = sqrt(((m1-1)*var(pc.feature.stats[indexes.for.feature.set]) + (m2-1)*var(pc.feature.stats[not.feature.set.indexes]))/(m1+m2-2))      
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
      lower.p = pt(t.stat, df=df, lower.tail=T)
      upper.p = pt(t.stat, df=df, lower.tail=F)        
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
pcgse_via_WMW = function(data, prcomp.output, pc.indexes, feature.set.indexes, feature.statistics, cor.adjustment) {
  
  num.feature.sets = length(feature.set.indexes)
  n= nrow(data)
  p.values = matrix(0, nrow=num.feature.sets, ncol=length(pc.indexes))  
  rownames(p.values) = names(feature.set.indexes)
  feature.set.statistics = matrix(T, nrow=num.feature.sets, ncol=length(pc.indexes))    
  rownames(feature.set.statistics) = names(feature.set.indexes)    
  
  for (i in 1:num.feature.sets) {
    indexes.for.feature.set = feature.set.indexes[[i]]
    m1 = length(indexes.for.feature.set)
    not.feature.set.indexes = which(!(1:ncol(data) %in% indexes.for.feature.set))
    m2 = length(not.feature.set.indexes)
    
    if (cor.adjustment) {            
      # compute sample correlation matrix for members of feature set
      cor.mat = cor(data[,indexes.for.feature.set])
      # compute the mean pair-wise correlation 
      mean.cor = (sum(cor.mat) - m1)/(m1*(m1-1))    
    }
    
    for (j in 1:length(pc.indexes)) {
      # get the feature-level statistics for this PC
      pc.feature.stats = feature.statistics[,j]
      # compute the rank sum statistic feature-level statistics
      wilcox.results = wilcox.test(x=pc.feature.stats[indexes.for.feature.set], y=pc.feature.stats[not.feature.set.indexes],
                                   alternative="two.sided", exact=F, correct=F)
      rank.sum = wilcox.results$statistic                
      if (cor.adjustment) {
        # Using correlation-adjusted formula from Wu et al.
        var.rank.sum = ((m1*m2)/(2*pi))*(asin(1) + (m2 - 1)*asin(.5) + (m1-1)*(m2-1)*asin(mean.cor/2) +(m1-1)*asin((mean.cor+1)/2))
      } else {        
        var.rank.sum = m1*m2*(m1+m2+1)/12
      }
      z.stat = (rank.sum - (m1*m2)/2)/sqrt(var.rank.sum)
      feature.set.statistics[i,j] = z.stat      
      # compute the p-value via a two-sided z-test
      lower.p = pnorm(z.stat, lower.tail=T)
      upper.p = pnorm(z.stat, lower.tail=F)        
      p.values[i,j] = 2*min(lower.p, upper.p)
    }
  } 
  
  # Build the result list
  results = list()
  results$p.values = p.values
  results$statistics = feature.set.statistics  
  
  return (results)
}
