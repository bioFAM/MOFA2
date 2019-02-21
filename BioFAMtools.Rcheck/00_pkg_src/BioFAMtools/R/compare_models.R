
################################################
## Functions to compare different BioFAModels ##
################################################


#' @title Plot the robustness of the latent factors across diferent trials
#' @name compare_factors
#' @description Different objects of \code{\link{BioFAModel}} are compared in terms of correlation between
#' their latent factors. The correlation is calculated only on those samples which are present in all models.
#' Ideally, the output should look like a block diagonal matrix, suggesting that all detected factors are robust under different initialisations.
#' If not, it suggests that some factors are weak and not captured by all models.
#' @param models a list containing \code{\link{BioFAModel}} objects.
#' @param comparison tye of comparison, either 'pairwise' or 'all'
#' @param ... extra arguments passed to pheatmap
#' @details TO-FILL
#' @return Plots a heatmap of correlation of Latent Factors in all models when 'comparison' is 'all'.
#' Otherwise, for each pair of models, a seperate heatmap is produced comparing one model againt the other.
#' The corresponding correlation matrix or list or pairwise correlation matrices is returned
#' @references fill this
#' @importFrom stats cor
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @export

compare_factors <- function(models, comparison = "all", show_rownames=FALSE, show_colnames=FALSE,...) {
  breaksList = seq(0, 1, by = 0.01)  # define a break list to fix the colour range irrespective of data range
  # Sanity checks
  if(!is.list(models))
    stop("'models' has to be a list")
  if (!all(sapply(models, function (l) class(l)=="BioFAModel")))
    stop("Each element of the the list 'models' has to be an instance of BioFAModel")
  if (!comparison %in% c("all", "pairwise"))
    stop("'comparison' has to be either 'all' or 'pairwise'")

  # give generic names if no names present
  if(is.null(names(models))) names(models) <- paste("model", 1: length(models), sep="")

  # TODO maybe find a better way to handle that
  # if(node=='Z'){
  #     Z <- getFactors(models)
  #   }
  #   else{
  #     SW <- getWeights(model)
  #     Z  <- do.call('rbind', SW)
  #     # concatenate weight matrices across views
  #   }

  # get latent factors
  LFs <- lapply(seq_along(models), function(modelidx){
    model <- models[[modelidx]]
    Z <- get_factors(model)
    # NOTE: concatenate all groups by default
    Z <- do.call(rbind, Z)
    Z
    })

  for (i in seq_along(LFs))
    colnames(LFs[[i]]) <- paste(names(models)[i], colnames(LFs[[i]]), sep="_")

  if (comparison=="all") {
    #get common samples between models
    common_samples <- Reduce(intersect,lapply(LFs, rownames))
    if(is.null(common_samples))
      stop("No common samples in all models for comparison")

    #subset LFs to common samples
    LFscommon <- Reduce(cbind, lapply(LFs, function(Z) Z[common_samples,, drop=FALSE]))

    # calculate correlation
    corLFs <- cor(LFscommon, use="complete.obs")

    #annotation by model
    # modelAnnot <- data.frame(model = rep(names(models), times=sapply(LFs, ncol)))
    # rownames(modelAnnot) <- colnames(LFscommon)

    #plot heatmap
    # TODO find a better way to plot heatmap with NAs
    corLFs[is.na(corLFs)]  = 0
    # if(is.null(main)) main <- "Absolute correlation between latent factors"
    if(length(unique(as.numeric(abs(corLFs))))>1){
    pheatmap(abs(corLFs), show_rownames = show_rownames,show_colnames = show_colnames,
             color = colorRampPalette(c("white",RColorBrewer::brewer.pal(9,name="YlOrRd")))(length(breaksList)), breaks=breaksList,
             # color=colorRampPalette(c("white", "orange" ,"red"))(100),
             # annotation_col = modelAnnot, main= main , ...)
             ...)
    } else warning("No plot produced as correlations consist of only one value")
    # return(corLFs)
  }

  if(comparison=="pairwise"){
    PairWiseCor <- lapply(seq_along(LFs[-length(LFs)]), function(i){
      LFs1<-LFs[[i]]
        sublist <- lapply((i+1):length(LFs), function(j){
          LFs2<-LFs[[j]]
          common_pairwise <- intersect(rownames(LFs1), rownames(LFs2))
          if(is.null(common_pairwise)) {
            warning(paste("No common samples between models",i,"and",j,"- No comparison possible"))
            NA
          }
          else{
          # if(is.null(main)) main <- paste("Absolute correlation between factors in model", i,"and",j)
          corLFs_pairs <- cor(LFs1[common_pairwise,], LFs2[common_pairwise,], use="complete.obs")
          corLFs_pairs[is.na(corLFs_pairs)]  = 0
          if(length(unique(abs(corLFs_pairs)))>1){
          pheatmap(abs(corLFs_pairs),color=colorRampPalette(c("white", "orange" ,"red"))(length(breaksList)), breaks=breaksList,show_rownames= show_rownames, show_colnames = show_colnames, ...)
          } else warning("No plot produced as correlations consist of only one value")
          corLFs_pairs
          }
        })
        names(sublist) <- names(models)[(i+1):length(LFs)]
        sublist
    })

    names(PairWiseCor) <- names(models[-length(models)])
    return(NULL)
    #return(PairWiseCor)
  }
}



#' @title Plot the robustness of the loading matrices, ie. the definition of the factors (how they impact features), across diferent trials
#' @name compare_weights
#' @description Different objects of \code{\link{BioFAModel}} are compared in terms of correlation between
#' their weights matrices. The correlations are computed between the columns of the concatenation of the weights matrices, only on those features which are present in all models.
#' Ideally, the output should look like a block diagonal matrix, suggesting that the definition of factors (how they impact features) is robust under different initialisations.
#' If not, it suggests that some factors are weak and not captured by all models.
#' @param models a list containing \code{\link{BioFAModel}} objects.
#' @param views : list of views or 'all'. Correlations computed on the concatenation of the loading matrices for the given views
#' @param comparison tye of comparison, either 'pairwise' or 'all'
#' @param ... extra arguments passed to pheatmap
#' @details TO-FILL
#' @return Plots a heatmap of correlation of weights columns in all models when 'comparison' is 'all'.
#' Otherwise, for each pair of models, a seperate heatmap is produced comparing one model againt the other.
#' The corresponding correlation matrix or list or pairwise correlation matrices is returned
#' @references fill this
#' @importFrom stats cor
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @export

compare_weights <- function(models, views = "all", comparison = "all", show_rownames=FALSE, show_colnames=FALSE,...) {
  breaksList = seq(0, 1, by = 0.001)
  # Sanity checks
  if(!is.list(models))
    stop("'models' has to be a list")
  if (!all(sapply(models, function (l) class(l)=="BioFAModel")))
    stop("Each element of the the list 'models' has to be an instance of BioFAModel")
  if (!comparison %in% c("all", "pairwise"))
    stop("'comparison' has to be either 'all' or 'pairwise'")

  # give generic names if no names present
  if(is.null(names(models))) names(models) <- paste("model", 1: length(models), sep="")

  # get latent factors
  LFs <- lapply(seq_along(models), function(modelidx){
    model <- models[[modelidx]]

    W <- get_weights(model)
    if (views == "all"){
      W <- do.call(rbind, W)
    }
    else{
      W <- do.call(rbind, W[views])
    }

    # NOTE: concatenate all groups by default
  })

  for (i in seq_along(LFs))
    colnames(LFs[[i]]) <- paste(names(models)[i], colnames(LFs[[i]]), sep="_")

  if (comparison=="all") {
    #get common samples between models
    common_features <- Reduce(intersect,lapply(LFs, rownames))
    if(is.null(common_features))
      stop("No common features in all models for comparison for the specified view")

    #subset LFs to common samples
    LFscommon <- Reduce(cbind, lapply(LFs, function(W) W[common_features,, drop=FALSE]))

    # calculate correlation
    corLFs <- cor(LFscommon, use="complete.obs")

    #annotation by model
    # modelAnnot <- data.frame(model = rep(names(models), times=sapply(LFs, ncol)))
    # rownames(modelAnnot) <- colnames(LFscommon)

    #plot heatmap
    # TODO find a better way to plot heatmap with NAs
    corLFs[is.na(corLFs)]  = 0
    # if(is.null(main)) main <- "Absolute correlation between latent factors"
    if(length(unique(as.numeric(abs(corLFs))))>1){
      pheatmap(abs(corLFs), show_rownames = show_rownames,show_colnames = show_colnames,
               color = colorRampPalette(c("white",RColorBrewer::brewer.pal(n = 7, name = "YlOrRd")))(length(breaksList)),
               # color=colorRampPalette(c("white", "orange" ,"red"))(100),
               # annotation_col = modelAnnot, main= main , ...)
               ...)
    } else warning("No plot produced as correlations consist of only one value")
    # return(corLFs)
  }

  if(comparison=="pairwise"){
    PairWiseCor <- lapply(seq_along(LFs[-length(LFs)]), function(i){
      LFs1<-LFs[[i]]
      sublist <- lapply((i+1):length(LFs), function(j){
        LFs2<-LFs[[j]]
        common_pairwise <- intersect(rownames(LFs1), rownames(LFs2))
        if(is.null(common_pairwise)) {
          warning(paste("No common samples between models",i,"and",j,"- No comparison possible"))
          NA
        }
        else{
          # if(is.null(main)) main <- paste("Absolute correlation between factors in model", i,"and",j)
          # corLFs_pairs <- cor(LFs1[common_pairwise,], LFs2[common_pairwise,], use="complete.obs")
          corLFs_pairs <- cor(LFs1, LFs2, use="complete.obs")
          corLFs_pairs[is.na(corLFs_pairs)]  = 0
          if(length(unique(abs(corLFs_pairs)))>1){
            pheatmap(abs(corLFs_pairs),
                     color = colorRampPalette(c("white",RColorBrewer::brewer.pal(n = 7, name = "YlOrRd")))(length(breaksList)),
                     show_rownames=F, show_colnames = F, ...)
          } else warning("No plot produced as correlations consist of only one value")
          corLFs_pairs
        }
      })
      names(sublist) <- names(models)[(i+1):length(LFs)]
      sublist
    })

    names(PairWiseCor) <- names(models[-length(models)])
    return(NULL)
    #return(PairWiseCor)
  }
}

#' @title Compare different trained \code{\link{BioFAModel}} objects in terms of the final value of the ELBO statistics and number of inferred factors
#' @name compareModels
#' @description Different objects of \code{\link{BioFAModel}} are compared in terms of the final value of the ELBO statistics.
#' For model selection the model with the highest ELBO value is selected.
#' @param models a list containing \code{\link{BioFAModel}} objects.
#' @export

compare_models <- function(models, show_modelnames = FALSE) {
  # Sanity checks
  if(!is.list(models))
    stop("'models' has to be a list")
  if (!all(sapply(models, function (l) class(l)=="BioFAModel")))
    stop("Each element of the the list 'models' has to be an instance of BioFAModel")

  elbo_vals <- sapply(models, get_elbo)
  n_factors <- sapply(models, function(m) get_dimensions(m)$K)
  if(is.null(names(models))) names(models) <- paste0("model_", seq_along(models))
  df <- data.frame(ELBO=elbo_vals, model = names(models))
  gg <- ggplot(df, aes(x=model,y=n_factors, fill=ELBO)) + geom_bar(stat="identity")+
    ylab("Number of inferred factors")
  if(show_modelnames) gg <- gg + theme(axis.text.x=element_text(angle=60, vjust=1, hjust=1))
  else gg <- gg + theme(axis.text.x=element_blank())
  return(gg)
}

#' @title Select a model from a list of trained \code{\link{BioFAModel}} objects based on the best ELBO value
#' @name select_model
#' @description Different objects of \code{\link{BioFAModel}} are compared in terms of the final value of the ELBO statistics

#' and the model with the highest ELBO value is selected.
#' @param models a list containing \code{\link{BioFAModel}} objects.
#' @export

select_model <- function(models, plotit = TRUE) {
  # Sanity checks
  if(!is.list(models))
    stop("'models' has to be a list")
  if (!all(sapply(models, function (l) class(l)=="BioFAModel")))
    stop("Each element of the the list 'models' has to be an instance of BioFAModel")

  elbo_vals <- sapply(models, get_elbo)
  if(plotit) compare_models(models)
  models[[which.max(elbo_vals)]]
}
