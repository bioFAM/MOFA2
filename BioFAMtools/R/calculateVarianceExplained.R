
#' @title Calculate variance explained by the model
#' @name calculateVarianceExplained
#' @description Method to calculate variance explained by the BioFAModel for each view and latent factor. \cr
#' As a measure of variance explained for gaussian data we adopt the coefficient of determination (R2). \cr
#' For non-gaussian views the calculations are based on the normally-distributed pseudo-data 
#' (for more information on the non-gaussian model see Supplementary Methods of the MOFA paper or Seeger & Bouchard, 2012).
#' @param object a \code{\link{BioFAModel}} object.
#' @param views character vector with the view names, or numeric vector with view indexes. Default is 'all'
#' @param factors character vector with the factor names, or numeric vector with the factor indexes. Default is 'all'
#' @param include_intercept include the intercept factor for calculation of variance explained (only used when an intercept was learned)
#' @param batches character vector with the batch names, or numeric vector with batch indexes. Default is 'all'
#' @details This function takes a trained BioFAModel as input and calculates for each view the coefficient of determination (R2),
#' i.e. the proportion of variance in the data explained by the BioFAM factor(s) (both jointly and for each individual factor). 
#' In case of non-Gaussian data the variance explained on the Gaussian pseudo-data is calculated. 
#' @return a list with matrices with the amount of variation explained per factor and view, and optionally total variance explained per view and variance explained by each feature alone
#' @export
calculateVarianceExplained <- function(object, views = "all", factors = "all", include_intercept = TRUE, batches = "all", ...) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  
  # check whether the intercept was learned
  if(!object@ModelOpts$learnIntercept & include_intercept) {
    include_intercept <- FALSE
    # warning("No intercept was learned in BioFAM.\n Intercept is not included in the model prediction.")
  }
  
  # Define views
  if (paste0(views, sep="", collapse="") == "all") { 
    views <- viewNames(object) 
  } else {
    stopifnot(all(views %in% viewNames(object)))  
  }
  M <- length(views)

  # Define batches
  if (paste0(batches, sep="", collapse="") == "all") { 
    batches <- batchNames(object) 
  } else {
    stopifnot(all(batches %in% batchNames(object)))  
  }
  H <- length(batches)

  # Define factors
  if (paste0(factors, collapse="") == "all") { 
    factors <- factorNames(object) 
  } else if (is.numeric(factors)) {
    if (include_intercept == T) {
      factors <- factorNames(object)[factors+1] 
    } else {
      factors <- factorNames(object)[factors]
    }
  } else { 
    stopifnot(all(factors %in% factorNames(object))) 
  }
  factors <- factors[factors!="intercept"]
  K <- length(factors)

  # Collect relevant expectations
  W <- getWeights(object, views, factors)
  Z <- getFactors(object, batches, factors)
  Y <- getExpectations(object, "Y")  # for non-Gaussian likelihoods the pseudodata is considered
  
  # Calulcate feature-wise means as null model
  FeatureMean <- lapply(views, function(m) {
    lapply(batches, function(h) {
      apply(Y[[m]][[h]], 2, mean, na.rm=T) 
    })
  })
  names(FeatureMean) <- views
  for (view in views) { names(FeatureMean[[view]]) <- batches }

  # Sweep out the feature-wise mean to calculate null model residuals
  resNullModel <- lapply(views, function(m) {
    lapply(batches, function(h) {
      sweep(Y[[m]][[h]], 2, FeatureMean[[m]][[h]], "-")
    })
  })
  names(resNullModel) <- views
  for (view in views) { names(resNullModel[[view]]) <- batches }
  browser()
  # replace masked values on Z by 0 (so that they do not contribute to predictions)
  for (batch in batches) { Z[[batch]][is.na(Z[[batch]])] <- 0 }
    
  # Calculate predictions under the MOFA model using all (non-intercept) factors
  Ypred_m <- lapply(views, function(m) Z%*%t(W[[m]])); names(Ypred_m) <- views

  # Calculate predictions under the MOFA model using each (non-intercept) factors on its own
  Ypred_mk <- lapply(views, function(m) {
                      ltmp <- lapply(factors, function(k) Z[,k]%*%t(W[[m]][,k]) )
                      names(ltmp) <- factors
                      ltmp
                    })
  names(Ypred_mk) <- views
  
  # If an intercept is included, regress out the intercept from the data
  if (include_intercept==T) {
      intercept <- getWeights(object,views,"intercept")
      Y <- lapply(views, function(m) sweep(Y[[m]],2,intercept[[m]],"-"))
      names(Y) <- views
  }

  # Calculate coefficient of determination per view
  fvar_m <- sapply(views, function(m) 1 - sum((Y[[m]]-Ypred_m[[m]])**2, na.rm=T) / sum(resNullModel[[m]]**2, na.rm=T))
   
  # Calculate coefficient of determination per factor and view
  tmp <- sapply(views, function(m) {
    sapply(factors, function(k) {
      1 - sum((Y[[m]]-Ypred_mk[[m]][[k]])**2, na.rm=T) / sum(resNullModel[[m]]**2, na.rm=T) 
    })
  }) 
  fvar_mk <- matrix(tmp, ncol=length(views), nrow=length(factors))
  
  # Set names
  names(fvar_m) <- views
  colnames(fvar_mk) <- views
  rownames(fvar_mk) <- factors

  # Store results
  R2_list <- list(R2Total = fvar_m, R2PerFactor = fvar_mk)

  return(R2_list)
 
}

#' @title Plot variance explained by the model
#' @name plotVarianceExplained
#' @description Method to plot variance explained (R-squared) by the MOFA model for each view and latent factor. \cr
#' As a measure of variance explained for gaussian data we adopt the coefficient of determination (R2). \cr
#' For details on the computation see the help of the \code{\link{calculateVarianceExplained}} function
#' @param object a \code{\link{MOFAmodel}} object.
#' @param cluster logical indicating whether to do hierarchical clustering on the plot
#' @param ... extra arguments to be passed to \code{\link{calculateVarianceExplained}}
#' @return ggplot object
#' @import pheatmap ggplot2 reshape2
#' @importFrom cowplot plot_grid
#' @export
plotVarianceExplained <- function(object, cluster = T, ...) {
  
  # Calculate Variance Explained
  R2_list <- calculateVarianceExplained(object, ...)
  fvar_m <- R2_list$R2Total
  fvar_mk <- R2_list$R2PerFactor
  
  ## Plot variance explained by factor ##
  
  # convert matrix to data frame for ggplot2  
  fvar_mk_df <- reshape2::melt(fvar_mk, varnames=c("factor","view"))
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
  fvar_m_df <- data.frame(view=names(fvar_m), R2=fvar_m)
  
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


