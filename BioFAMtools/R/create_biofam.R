
#' @title create a BioFAM object
#' @name create_biofam
#' @description Method to create a \code{\link{BioFAModel}} object
#' @param data TO SPECIFY
#' @return Returns an untrained \code{\link{BioFAModel}} object
#' @export
create_biofam <- function(data) {
  
  if (is(data,"MultiAssayExperiment")) {
    stop("Not functional")
    # message("Creating BioFAM object from a MultiAssayExperiment object...")
    # object <- .createBioFAMobjectFromMAE(data)
    
  } else if (is(data,"SummarizedExperiment")) {
    stop("Not functional")
    
  } else if (is(data,"data.frame")) {
      message("Creating BioFAM object from a dataframe...")
    
      object <- .create_biofam_from_dataframe(data)
      
      # Set dimensionalities
      object@dimensions[["M"]] <- length(unique(data$feature_group))
      object@dimensions[["D"]] <- sapply(unique(data$feature_group), function(m) length(unique(data[data$feature_group==m,]$feature)))
      object@dimensions[["P"]] <- length(unique(data$sample_group))
      object@dimensions[["N"]] <- sapply(unique(data$sample_group), function(p) length(unique(data[data$sample_group==p,]$sample)))
      object@dimensions[["K"]] <- 0
      
      # Set view names
      views_names(object) <- unique(data$feature_group)
      
      # Set sample group names
      groups_names(object) <- unique(data$sample_group)
      
  } else {
    stop("Error: input data has to be provided either as .....")
  }
  
  print(object)
  
  return(object)
}


# (Hidden) function to initialise a BioFAModel object using a MultiAssayExperiment
# #' @import MultiAssayExperiment
# .create_biofam_from_MultiAssayExperiment <- function(data) {
# 
#   # Initialise BioFAM object
#   object <- new("BioFAModel")
#   object@status <- "untrained"
#   object@input_data <- data
#   
#   # Re-arrange data for training in BioFAM to matrices, fill in NAs and store in TrainData slot
#   object@TrainData <- lapply(names(data), function(m) {
#     
#     # Extract general sample names
#     primary <- unique(sampleMap(data)[,"primary"])
#     
#     # Extract view
#     subdata <- assays(data)[[m]]
#     
#     # Rename view-specific sample IDs with the general sample names
#     stopifnot(colnames(subdata)==sampleMap(data)[sampleMap(data)[,"assay"]==m,"colname"])
#     colnames(subdata) <- sampleMap(data)[sampleMap(data)[,"assay"]==m,"primary"]
#     
#     # Fill subdata with NAs
#     subdata_filled <- .subset_augment(subdata,primary)
#     return(subdata_filled)
#   })
#   return(object)
# }


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
.create_biofam_from_dataframe<- function(data) {
  
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
  object@status <- "untrained"
  object@input_data <- data
  
  return(object)
}
