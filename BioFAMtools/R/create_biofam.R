
#' @title Initialize a BioFAM object
#' @name create_biofam_object
#' @description Method to initialize a \code{\link{BioFAModel}} object with a multi-omics data set.
#' @param data either a \code{MultiAssayExperiment} or a list of matrices with features as rows and samples as columns.
#' @details
#' If the multi-omics data is provided as a list of matrices, please make sure that features 
#' are stored as rows and samples are stored as columns. \cr 
#' If the matrices have sample names, we will use them to match the different matrices, filling the corresponding missing values. \cr
#' If matrices have no column names, all matrices must have the same number of columns, and you are responsible for filling any missing values.
#' @return Returns an untrained \code{\link{BioFAModel}} object
#' @export
create_biofam_object <- function(data, samples_groups = NULL) {
  
  if (is(data, "MultiAssayExperiment")) {
    message("Creating BioFAM object from a MultiAssayExperiment object...")
  } else if (is(data, "list")) {
    message("Creating BioFAM object from list of matrices,\n please make sure that samples are columns and features are rows...\n")
    data <- .create_mae_object_from_list(data)
  } else {
    stop("Error: input data has to be provided either as a list of matrices or as a MultiAssayExperiment object \n")
  }
  object <- .create_biofam_object_from_mae(data)

  # ALl samples are coming from one group if samples_groups is not provided
  if (is.null(samples_groups)) {
    samples_groups <- rep("group_1", ncol(object@training_data[[1]]))
  }
  object <- .split_samples_by_group(object, samples_groups)
  
  # Set dimensionalities
  object@dimensions[["M"]] <- length(object@training_data)
  object@dimensions[["P"]] <- length(object@training_data[[1]])
  object@dimensions[["N"]] <- sapply(object@training_data[[1]], ncol)
  object@dimensions[["D"]] <- sapply(object@training_data, function(p) nrow(p[[1]]))
  object@dimensions[["K"]] <- 0
  
  # Set view names
  if(!is.null(names(data))) {
    views_names(object) <- names(data) 
  } else { 
    views_names(object) <- paste("view", 1:length(object@training_data), sep="_")
    warning(paste0("View names are not specified in the data, renaming them to: ",
                   paste0("view_", 1:length(object@training_data), collapse=", "), "\n"))
  }

  # Set group names
  groups_names(object) <- unique(samples_groups)
  
  # Set feature names
  for (m in 1:object@dimensions[["M"]]) {
    if (is.null(rownames(object@training_data[[m]][[1]]))) {
      warning(sprintf("Feature names are not specified for view %d, using default: feature1_v%d,feature2_v%d...\n",m,m,m))
      for (p in 1:object@dimensions[["P"]]) {
        rownames(object@training_data[[m]][[p]]) <- paste0("feature_",1:nrow(object@training_data[[m]]),"_v",m )
      }
    }
  }
  
  # Feature names are already set in .createBioFAMobjectFromList in case they are missing
  
  return(object)
}

.split_samples_by_group <- function(object, samples_groups) {

  groups_names <- unique(samples_groups)

  object@training_data <- lapply(names(object@training_data), function(m) {
    groups_list <- lapply(groups_names, function(p) {
      object@training_data[[m]][,samples_groups == p]
    })
    names(groups_list) <- groups_names
    return(groups_list)
  })

  return(object)

}

# (Hidden) function to initialise a BioFAModel object using a MultiAssayExperiment
#' @import MultiAssayExperiment
.create_biofam_object_from_mae <- function(data) {

  # Initialise BioFAM object
  object <- new("BioFAModel")
  object@status <- "untrained"
  object@input_data <- data
  
  # Re-arrange data for training in BioFAM to matrices, fill in NAs and store in training_data slot
  object@training_data <- lapply(names(data), function(m) {
    
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
  names(object@training_data) <- names(data)
  return(object)
}

#' @importFrom Biobase ExpressionSet
#' @import MultiAssayExperiment
.create_mae_object_from_list <- function(data) {
  
  # Fetch or assign sample names
  samples <- lapply(data[[1]], colnames)
  if (is.null(samples)) {
    N <- unique(sapply(data,ncol))
    if (length(N)>1) {
      stop("If the matrices have no column (samples) names that can be used to match the different views,
           all matrices must have the same number of columns")
    }
    warning("Sample names are not specified, using default: sample_1,sample_2...
             Make sure the columns match between data matrices or provide sample names! \n")
    samples <- paste0("sample_",1:N)
    for (m in 1:length(data)) { colnames(data[[m]]) <- samples }
  }
  
  # Create ExpressionSet for each assay
  expressionset_list <- list()
  for (i in 1:length(data)) {
    expressionset_list[[i]] <- ExpressionSet(assayData=data[[i]])
  }
  names(expressionset_list) <- names(data)
  
  # Create MultiAssayExperiment
  MAE <- MultiAssayExperiment(
    experiments = expressionset_list
  )
}
  
# (Hidden) function to initialise a BioFAModel object using a list of matrices
# .createBioFAMobjectFromList <- function(data) {
#   
#   # Initialise BioFAM object
#   object <- new("BioFAModel")
#   object@Status <- "untrained"
#   
#   # Fetch or assign sample names
#   samples <- Reduce(union, lapply(data, colnames))
#   if (is.null(samples)) {
#     N <- unique(sapply(data,ncol))
#     if (length(N)>1) { 
#       stop("If the matrices have no column (samples) names that can be used to match the different views,
#            all matrices must have the same number of columns")
#     }
#     warning("Sample names are not specified, using default: sample_1,sample_2...
#              Make sure the columns match between data matrices or provide sample names! \n")
#     samples <- paste0("sample_",1:N)
#     for (m in 1:length(data)) { colnames(data[[m]]) <- samples }
#   }
#   
#   object@training_data <- lapply(data, function(view) .subset_augment(view, samples))
#   
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