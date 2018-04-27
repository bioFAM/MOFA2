
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
    message("Creating BioFAM object from a MultiAssayExperiment object...")
    object <- .createBioFAMobjectFromMAE(data)
  } else if (is(data,"list")) {
    message("Creating BioFAM object from list of matrices, please make sure that samples are columns and features are rows...")
    object <- .createBioFAMobjectFromList(data)
  } else {
    stop("Error: input data has to be provided either as a list of matrices or as a MultiAssayExperiment object")
  }
  
  # Set dimensionalities
  object@Dimensions[["M"]] <- length(object@TrainData)
  object@Dimensions[["N"]] <- ncol(object@TrainData[[1]])
  object@Dimensions[["D"]] <- sapply(object@TrainData, nrow)
  object@Dimensions[["K"]] <- 0
  
  # Set view names
  if(!is.null(names(data))) {
    viewNames(object) <- names(data) 
  } else { 
    viewNames(object) <- paste("view", 1:length(object@TrainData), sep="_")
    warning(paste0("View names are not specified in data, renaming them to: ",paste("view",1:length(object@TrainData), collapse=" "), "\n"))
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


# (Hidden) function to initialise a BioFAModel object using a list of matrices
.createBioFAMobjectFromList <- function(data) {
  
  # Initialise BioFAM object
  object <- new("BioFAModel")
  object@Status <- "untrained"
  
  # Fetch or assign sample names
  samples <- Reduce(union, lapply(data, colnames))
  if (is.null(samples)) {
    N <- unique(sapply(data,ncol))
    if (length(N)>1) { 
      stop("If the matrices have no column (samples) names that can be used to match the different views, all matrices must have the same number of columns")
    }
    samples <- as.character(1:N)
    for (m in 1:length(data)) { colnames(data[[m]]) <- samples }
  }
  
  object@TrainData <- lapply(data, function(view) .subset_augment(view, samples))
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