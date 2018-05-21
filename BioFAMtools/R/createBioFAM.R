
#' @title Initialize a BioFAM object
#' @name createBioFAMbject
#' @description Method to initialize a \code{\link{BioFAModel}} object with a multi-omics data set.
#' @param data either a \code{MultiAssayExperiment} or a list of matrices with features as rows and samples as columns.
#' @details
#' If the multi-omics data is provided as a list of matrices, please make sure that features 
#' are stored as rows and samples are stored as columns. \cr 
#' If the matrices have sample names, we will use them to match the different matrices, filling the corresponding missing values. \cr
#' If matrices have no column names, all matrices must have the same number of columns, and you are responsible for filling any missing values.
#' @return Returns an untrained \code{\link{BioFAModel}} object
#' @export
createBioFAMobject <- function(data) {
  
  if (is(data,"MultiAssayExperiment")) {
    stop("Not functional")
    # message("Creating BioFAM object from a MultiAssayExperiment object...")
    # object <- .createBioFAMobjectFromMAE(data)
    
  } else if (is(data,"SummarizedExperiment")) {
    stop("Not functional")
    
  } else if (is(data,"data.frame")) {
      message("Creating BioFAM object from a data.frame...")
    
      object <- .createBioFAMobjectFromDataFrame(data)
      
      # Set dimensionalities
      object@Dimensions[["M"]] <- length(unique(object$feature_groups))
      object@Dimensions[["N"]] <- length(unique(object$sample))
      object@Dimensions[["D"]] <- length(unique(object$feature))
      object@Dimensions[["K"]] <- 0
      
      # Set view names
      viewNames(data) <- unique(object$feature_group)
      
  } else {
    stop("Error: input data has to be provided either as a list of matrices or as a MultiAssayExperiment object")
  }
  
  print(object)
  
  return(object)
}


# (Hidden) function to initialise a BioFAModel object using a MultiAssayExperiment
#' @import MultiAssayExperiment
.createBioFAMobjectFromMAE <- function(data) {

  # Initialise BioFAM object
  object <- new("BioFAModel")
  object@Status <- "untrained"
  object@InputData <- data
  
  # Re-arrange data for training in BioFAM to matrices, fill in NAs and store in TrainData slot
  object@TrainData <- lapply(names(data), function(m) {
    
    # Extract general sample names
    primary <- unique(sampleMap(data)[,"primary"])
    
    # Extract view
    subdata <- assays(data)[[m]]
    
    # Rename view-specific sample IDs with the general sample names
    stopifnot(colnames(subdata)==sampleMap(data)[sampleMap(data)[,"assay"]==m,"colname"])
    colnames(subdata) <- sampleMap(data)[sampleMap(data)[,"assay"]==m,"primary"]
    
    # Fill subdata with NAs
    subdata_filled <- .subset_augment(subdata,primary)
    return(subdata_filled)
  })
  return(object)
}


# (Hidden) function to fill NAs for missing samples
.subset_augment<-function(mat, pats) {
  pats <- unique(pats)
  mat <- t(mat)
  aug_mat<-matrix(NA, ncol=ncol(mat), nrow=length(pats))
  aug_mat<-mat[match(pats,rownames(mat)),,drop=FALSE]
  rownames(aug_mat)<-pats
  colnames(aug_mat)<-colnames(mat)
  return(t(aug_mat))
}

# (Hidden) function to initialise a BioFAModel object using a Dataframe
.createBioFAMFromDataFrame <- function(data) {
  
  # Quality controls
  stopifnot(all(colnames(data) %in% (c("sample","feature","value","sample_group","feature_group"))))
  stopifnot(all(is.numeric(data$value)))
  
  # Convert 'sample' and'feature' columns to factors
  if (!is.factor(data$sample))
    data$sample <- as.factor(data$sample)
  if (!is.factor(data$feature))
    data$feature <- as.factor(data$feature)
  
  # Convert 'sample_group' columns to factors
  if (!"sample_group" %in% colnames(data)) {
    data$sample_group <- as.factor("1")
  } else {
    data$sample_group <- as.factor(data$sample_group)
  }
  
  # Convert 'feature_group' columns to factors
  if (!"feature_group" %in% colnames(data)) {
    data$feature_group <- as.factor("1")
  } else {
    data$feature_group <- as.factor(data$feature_group)
  }
  
  # Initialise BioFAM object
  object <- new("BioFAModel")
  object@Status <- "untrained"
  object@InputData <- data
  
  return(object)
}
