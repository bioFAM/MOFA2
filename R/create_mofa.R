
#' @title create a MOFA object
#' @name create_mofa
#' @description Method to create a \code{\link{MOFA}} object. Depending on the input data format, this method calls one of the following functions:
#' \itemize{
#'   \item{\strong{long data.frame}: }{\code{\link{create_mofa_from_df}}}
#'   \item{\strong{List of matrices}: }{\code{\link{create_mofa_from_matrix}}}
#'   \item{\strong{MultiAssayExperiment}: }{\code{\link{create_mofa_from_MultiAssayExperiment}}}
#'   \item{\strong{Seurat}: }{\code{\link{create_mofa_from_Seurat}}}
#'   \item{\strong{SingleCellExperiment}: }{\code{\link{create_mofa_from_SingleCellExperiment}}}
#'   }
#'  Please read the documentation of the corresponding function for more details on your specific data format.
#' @param data one of the formats above
#' @param groups group information, only relevant when using the multi-group framework. 
#' @param ... further arguments that can be passed to the function depending on the inout data format.
#' See the dpcumentation of above functions for details.
#' @return Returns an untrained \code{\link{MOFA}} object
#' @export
#' @examples
#' # Using an existing simulated data with two groups and two views
#' file <- system.file("extdata", "test_data.RData", package = "MOFA2")
#' 
#' # Load data (in long data.frame format)
#' load(file) 
#' MOFAmodel <- create_mofa(dt)
create_mofa <- function(data, groups = NULL, ...) {
  
  # Creating MOFA object from a Seurat object
  if (is(data, "Seurat")) {
    
    message("Creating MOFA object from a Seurat object...")
    object <- create_mofa_from_Seurat(data, groups, ...)

  # Creating MOFA object from a SingleCellExperiment object
  } else if (is(data, "SingleCellExperiment")) {
    
    message("Creating MOFA object from a SingleCellExperiment object...")
    object <- create_mofa_from_SingleCellExperiment(data, groups, ...)
    

    # Creating MOFA object from a data.frame object
  } else if (is(data, "data.frame")) {
    
    message("Creating MOFA object from a data.frame...")
    object <- create_mofa_from_df(data)
    
    # Creating MOFA object from a (sparse) matrix object
  } else if (is(data, "list") && (length(data) > 0) && 
             (all(sapply(data, function(x) is(x, "matrix"))) || 
              all(sapply(data, function(x) is(x, "dgCMatrix"))) || 
              all(sapply(data, function(x) is(x, "dgTMatrix"))))) {
    
    message("Creating MOFA object from a list of matrices (features as rows, sample as columns)...\n")
    object <- create_mofa_from_matrix(data, groups)
    
  } else if(is(data, "MultiAssayExperiment")){

    object <- create_mofa_from_MultiAssayExperiment(data, groups, ...)

  } else {
    stop("Error: input data has to be provided as a list of matrices, a data frame or a Seurat object. Please read the documentation for more details.")
  }
  
  # Create sample metadata
  foo <- lapply(object@data[[1]], colnames)
  tmp <- data.frame(
    sample = unname(unlist(foo)),
    group = unlist(lapply(names(foo), function(x) rep(x, length(foo[[x]])) )),
    stringsAsFactors = FALSE
  )
  if (.hasSlot(object, "samples_metadata") && (length(object@samples_metadata) > 0)) {
    object@samples_metadata <- cbind(tmp, object@samples_metadata[match(tmp$sample, rownames(object@samples_metadata)),])
  } else {
    object@samples_metadata <- tmp
  }
  
  # Create features metadata
  tmp <- data.frame(
    feature = unname(unlist(lapply(object@data, function(x) rownames(x[[1]])))),
    view = unlist(lapply(seq_len(object@dimensions$M), function(x) rep(views_names(object)[[x]], object@dimensions$D[[x]]) )),
    stringsAsFactors = FALSE
  )
  if (.hasSlot(object, "features_metadata") && (length(object@features_metadata) > 0)) {
    object@features_metadata <- cbind(tmp, object@features_metadata[match(tmp$feature, rownames(object@features_metadata)),])
  } else {
    object@features_metadata <- tmp
  }
  
  # Do quality control
  object <- .quality_control(object)
  
  # print verbose messages
  if (length(unique(object@samples_metadata$group))>1) {
    message("\n# Multi-group mode requested.")
    message("\nThis is an advanced option, if this is the first time that you are running MOFA, we suggest that you try do some exploration first without specifying groups. Two important remarks:")
    message("\n - The aim of the multi-group framework is to identify the sources of variability *within* the groups. If your aim is to find a factor that 'separates' the groups, you DO NOT want to use the multi-group framework. Please see the FAQ in (https://biofam.github.io/MOFA2)") 
    message("\n - It is important to account for the group effect before selecting highly variable features (HVFs). We suggest that either you calculate HVFs per group and then take the union, or regress out the group effect before HVF selection")
  }
  return(object)
}

#' @title create a MOFA object from a MultiAssayExperiment object
#' @name create_mofa_from_MultiAssayExperiment
#' @description Method to create a \code{\link{MOFA}} object from a MultiAssayExperiment object
#' @param mae a MultiAssayExperiment object
#' @param groups a string specifying column name of the colData to use it as a group variable. 
#' Alternatively, a character vector with group assignment for every sample.
#' Default is \code{NULL} (no group structure).
#' @param save_metadata logical indicating whether to incorporate the metadata from the MultiAssayExperiment object into the MOFA object
#' @return Returns an untrained \code{\link{MOFA}} object
#' @export
create_mofa_from_MultiAssayExperiment <- function(mae, groups = NULL, save_metadata = FALSE) {
  
  # Sanity check
  if(!requireNamespace("MultiAssayExperiment", quietly = TRUE)){
    stop("Package \"MultiAssayExperiment\" is required but is not installed.", call. = FALSE)
  } else {

    # Re-arrange data for training in MOFA to matrices, fill in NAs
    data_list <- lapply(names(mae), function(m) {
      
      # Extract general sample names
      primary <- unique(MultiAssayExperiment::sampleMap(mae)[,"primary"])
      
      # Extract view
      subdata <- MultiAssayExperiment::assays(mae)[[m]]
      
      # Rename view-specific sample IDs with the general sample names
      stopifnot(colnames(subdata)==MultiAssayExperiment::sampleMap(mae)[MultiAssayExperiment::sampleMap(mae)[,"assay"]==m,"colname"])
      colnames(subdata) <- MultiAssayExperiment::sampleMap(mae)[MultiAssayExperiment::sampleMap(mae)[,"assay"]==m,"primary"]
      
      # Fill subdata with NAs
      subdata_filled <- .subset_augment(subdata, primary)
      return(subdata_filled)
    })
    
    # Define groups
    if (is(groups, 'character') && (length(groups) == 1)) {
      if (!(groups %in% colnames(MultiAssayExperiment::colData(mae))))
        stop(paste0(groups, " is not found in the colData of the MultiAssayExperiment.\n",
                    "If you want to use groups information from MultiAssayExperiment,\n",
                    "please ensure to provide a column name that exists. The columns of colData are:\n",
                    paste0(colnames(MultiAssayExperiment::colData(mae)), collapse = ", ")))
      groups <- MultiAssayExperiment::colData(mae)[,groups]
    }
    
    # If no groups provided, treat all samples as coming from one group
    if (is.null(groups)) {
      # message("No groups provided as argument, we assume that all samples belong to the same group.\n")
      groups <- rep("group1",  length(unique(MultiAssayExperiment::sampleMap(mae)[,"primary"])))
    }
    
    # Initialise MOFA object
    object <- new("MOFA")
    object@status <- "untrained"
    object@data <- .split_data_into_groups(data_list, groups)
    
    # groups_nms <- unique(as.character(groups))
    groups_nms <- names(object@data[[1]])
    
    # Set dimensionalities
    object@dimensions[["M"]] <- length(data_list)
    object@dimensions[["G"]] <- length(groups_nms)
    object@dimensions[["D"]] <- sapply(data_list, nrow)
    object@dimensions[["N"]] <- sapply(groups_nms, function(x) sum(groups == x))
    object@dimensions[["K"]] <- 0
    
    # Set view names
    views_names(object) <- names(mae)
    
    # Set samples group names
    groups_names(object) <- groups_nms
    
    # Set metadata
    if (isTRUE(save_metadata)) {
      # Samples metadata
      if (ncol(MultiAssayExperiment::colData(mae)) > 0) {
        object@samples_metadata <- data.frame(MultiAssayExperiment::colData(mae))
      }
      # No features metadata is typically contained
    }
    
    return(object)
  }
}


#' @title create a MOFA object from a data.frame object
#' @name create_mofa_from_df
#' @description Method to create a \code{\link{MOFA}} object from a data.frame object
#' @param df \code{data.frame} object with at most 5 columns: \code{sample}, \code{group}, \code{feature}, \code{view}, \code{value}. 
#'   The \code{group} column (optional) indicates the group of each sample when using the multi-group framework.
#'   The \code{view} column (optional) indicates the view of each feature when having multi-view data.
#' @return Returns an untrained \code{\link{MOFA}} object
#' @export
#' @examples
#' # Using an existing simulated data with two groups and two views
#' file <- system.file("extdata", "test_data.RData", package = "MOFA2")
#' 
#' # Load data (in long data.frame format)
#' load(file) 
#' MOFAmodel <- create_mofa_from_df(dt)
create_mofa_from_df <- function(df) {
  
  # Quality controls
  df <- as.data.frame(df)
  if (!"group" %in% colnames(df)) {
    # message('No "group" column found in the data.frame, we will assume a common group for all samples')
    df$group <- "single_group"
  }
  if (!"view" %in% colnames(df)) {
    # message('No "view" column found in the data.frame, we will assume a common view for all features')
    df$view <- "single_view"
  }
  stopifnot(all(colnames(df) %in% (c("sample","feature","value","group","view"))))
  stopifnot(all(is.numeric(df$value)))
  
  # Convert 'sample' and 'feature' columns to factors
  if (!is.factor(df$sample))
    df$sample <- as.factor(df$sample)
  if (!is.factor(df$feature))
    df$feature <- as.factor(df$feature)
  
  # Convert 'group' columns to factors
  if (!"group" %in% colnames(df)) {
    df$group <- factor("group1")
  } else {
    df$group <- factor(df$group)
  }
  
  # Convert 'view' columns to factors
  if (!"view" %in% colnames(df)) {
    df$view <- factor("view1")
  } else {
    df$view <- factor(df$view)
  }
  
  data_matrix <- list()
  for (m in levels(df$view)) {
    data_matrix[[m]] <- list()
    features <- as.character( unique( df[df$view==m,"feature"] ) )
    for (g in levels(df$group)) {
      samples <- as.character( unique( df[df$group==g,"sample"] ) )
      Y <- df[df$view==m & df$group==g,]
      Y$sample <- factor(Y$sample, levels=samples)
      Y$feature <- factor(Y$feature, levels=features)
      if (nrow(Y)==0) {
        data_matrix[[m]][[g]] <- matrix(as.numeric(NA), ncol=length(samples), nrow=length(features))
        rownames(data_matrix[[m]][[g]]) <- features
        colnames(data_matrix[[m]][[g]]) <- samples
      } else {
        data_matrix[[m]][[g]] <- .df_to_matrix( reshape2::dcast(Y, feature~sample, value.var="value", fill=NA, drop=FALSE) )
      }
    }
  }
  
  # Create MOFA object
  object <- new("MOFA")
  object@status <- "untrained"
  object@data <- data_matrix
  
  # Set dimensionalities
  object@dimensions[["M"]] <- length(levels(df$view))
  object@dimensions[["D"]] <- sapply(levels(df$view), function(m) length(unique(df[df$view==m,]$feature)))
  object@dimensions[["G"]] <- length(levels(df$group))
  object@dimensions[["N"]] <- sapply(levels(df$group), function(g) length(unique(df[df$group==g,]$sample)))
  object@dimensions[["K"]] <- 0
  
  # Set view names
  views_names(object) <- levels(df$view)
  
  # Set group names
  groups_names(object) <- levels(df$group)
  
  return(object)
}


#' @title create a MOFA object from a SingleCellExperiment object
#' @name create_mofa_from_SingleCellExperiment
#' @description Method to create a \code{\link{MOFA}} object from a SingleCellExperiment object
#' @param sce SingleCellExperiment object
#' @param groups a string specifying column name of the colData to use it as a group variable. 
#' Alternatively, a character vector with group assignment for every sample.
#' Default is \code{NULL} (no group structure).
#' @param assay assay to use, default is \code{logcounts}.
#' @param save_metadata logical indicating whether to incorporate the metadata from the SingleCellExperiment object into the MOFA object
#' @return Returns an untrained \code{\link{MOFA}} object
#' @export
create_mofa_from_SingleCellExperiment <- function(sce, groups = NULL, assay = "logcounts", save_metadata = FALSE) {
  
  # Check is SingleCellExperiment is installed
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package \"SingleCellExperiment\" is required but is not installed.", call. = FALSE)
  }
  else if(!requireNamespace("SummarizedExperiment", quietly = TRUE)){
    stop("Package \"SummarizedExperiment\" is required but is not installed.", call. = FALSE)
  } else {
    stopifnot(assay%in%names(SummarizedExperiment::assays(sce)))
    
    # Define groups of cells
    if (is.null(groups)) {
      # message("No groups provided as argument... we assume that all samples are coming from the same group.\n")
      groups <- rep("group1", dim(sce)[2])
    } else {
      if (is(groups,'character')) {
        if (length(groups) == 1) {
          stopifnot(groups %in% colnames(colData(sce)))
          groups <- colData(sce)[,groups]
        } else {
          stopifnot(length(groups) == ncol(sce))
        }
      } else {
        stop("groups wrongly specified. Please see the documentation and the examples")
      }
    }
   
    # Extract data matrices
    data_matrices <- list( .split_sce_into_groups(sce, groups, assay) )
    names(data_matrices) <- assay
    
    # Create MOFA object
    object <- new("MOFA")
    object@status <- "untrained"
    object@data <- data_matrices
    
    # Define dimensions
    object@dimensions[["M"]] <- length(assay)
    object@dimensions[["D"]] <- vapply(data_matrices, function(m) nrow(m[[1]]), 1L)
    object@dimensions[["G"]] <- length(data_matrices[[1]])
    object@dimensions[["N"]] <- vapply(data_matrices[[1]], function(g) ncol(g), 1L)
    object@dimensions[["K"]] <- 0
    
    # Set views & groups names
    groups_names(object) <- as.character(names(data_matrices[[1]]))
    views_names(object)  <- assay
    
    # Set metadata
    if (isTRUE(save_metadata)) {
      object@samples_metadata <- as.data.frame(colData(sce))
      object@features_metadata <- as.data.frame(rowData(sce))
    }
    
    return(object)
  }
}

#' @title create a MOFA object from a Seurat object
#' @name create_mofa_from_Seurat
#' @description Method to create a \code{\link{MOFA}} object from a Seurat object
#' @param seurat Seurat object
#' @param groups a string specifying column name of the samples metadata to use it as a group variable. 
#' Alternatively, a character vector with group assignment for every sample.
#' Default is \code{NULL} (no group structure).
#' @param assays assays to use, default is \code{NULL}, it fetched all assays available
#' @param slot assay slot to be used such as scale.data or data (default).
#' @param features a list with vectors, which are used to subset features, with names corresponding to assays; a vector can be provided when only one assay is used
#' @param save_metadata logical indicating whether to incorporate the metadata from the Seurat object into the MOFA object
#' @return Returns an untrained \code{\link{MOFA}} object
#' @export
create_mofa_from_Seurat <- function(seurat, groups = NULL, assays = NULL, slot = "data", features = NULL, save_metadata = FALSE) {
  
  # Check is Seurat is installed
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package \"Seurat\" is required but is not installed.", call. = FALSE)
  } else {
  
    # Define assays
    if (is.null(assays)) {
      assays <- Seurat::Assays(seurat)
      message(paste0("No assays specified, using all assays by default: ", paste(assays,collapse=" ")))
    } else {
      stopifnot(assays%in%Seurat::Assays(seurat))
    }
    
    # Define groups of cells
    if (is(groups, 'character') && (length(groups) == 1)) {
      if (!(groups %in% colnames(seurat@meta.data)))
        stop(paste0(groups, " is not found in the Seurat@meta.data.\n",
                    "please ensure to provide a column name that exists. The columns of meta data are:\n",
                    paste0(colnames(seurat@meta.data), sep = ", ")))
      groups <- seurat@meta.data[,groups]
    }
    
    # If features to subset are provided,
    # make sure they are a list with respective views (assays) names.
    # A vector is accepted if there's one assay to be used
    if (is(features, "list")) {
      if (!is.null(features) && !all(names(features) %in% assays)) {
        stop("Please make sure all the names of the features list correspond to views (assays) names being used for the model")
      }
    } else {
      # By default select highly variable features if present in the Seurat object
      if (is.null(features)) {
        message("No features specified, using variable features from the Seurat object...")
        features <- lapply(assays, function(i) seurat@assays[[i]]@var.features)
        names(features) <- assays
        if (any(sapply(features,length)==0)) stop("No list of features provided and variable features are detected in the Seurat object")
      } else if (all(is(features, "character"))) {
        features <- list(features)
        names(features) <- assays
      } else {
        stop("Features not recognised. Please either provide a list of features (per assay) or calculate variable features in the Seurat object")
      }
    }
    
    # If no groups provided, treat all samples as coming from one group
    if (is.null(groups)) {
      # message("No groups provided as argument... we assume that all samples are coming from the same group.\n")
      groups <- rep("group1", dim(seurat)[2])
    }
    
    # Extract data matrices
    data_matrices <- lapply(assays, function(i) 
      .split_seurat_into_groups(seurat, groups = groups, assay = i, slot = slot, features = features[[i]]))
    names(data_matrices) <- assays
    
    # Create MOFA object
    object <- new("MOFA")
    object@status <- "untrained"
    object@data <- data_matrices
    
    # Define dimensions
    object@dimensions[["M"]] <- length(assays)
    object@dimensions[["D"]] <- vapply(data_matrices, function(m) nrow(m[[1]]), 1L)
    object@dimensions[["G"]] <- length(data_matrices[[1]])
    object@dimensions[["N"]] <- vapply(data_matrices[[1]], function(g) ncol(g), 1L)
    object@dimensions[["K"]] <- 0
    
    # Set views & groups names
    groups_names(object) <- as.character(names(data_matrices[[1]]))
    views_names(object)  <- assays
    
    # Set metadata
    if (isTRUE(save_metadata)) {
      object@samples_metadata <- seurat@meta.data
      object@features_metadata <- do.call(rbind, lapply(assays, function(a) seurat@assays[[a]]@meta.features))
    }
    
    return(object)
  }
}


#' @title create a MOFA object from a a list of matrices
#' @name create_mofa_from_matrix
#' @description Method to create a \code{\link{MOFA}} object from a list of matrices
#' @param data A list of matrices, where each entry corresponds to one view.
#'   Samples are stored in columns and features in rows.
#'   Missing values must be filled in prior to creating the MOFA object (see for example the CLL tutorial)
#' @param groups A character vector with group assignment for every sample. Default is \code{NULL}, no group structure.
#' @return Returns an untrained \code{\link{MOFA}} object
#' @export
#' @examples 
#' m <- make_example_data()
#' create_mofa_from_matrix(m$data)

create_mofa_from_matrix <- function(data, groups = NULL) {
  
  # Quality control: check that the matrices are all numeric
  stopifnot(all(sapply(data, function(g) all(is.numeric(g)))) || all(sapply(data, function(x) class(x) %in% c("dgTMatrix", "dgCMatrix"))))
  
  # Quality control: check that all matrices have the same samples
  tmp <- lapply(data, function(m) colnames(m))
  if(length(unique(sapply(tmp,length)))>1)
    stop("Views have different number of samples (columns)... please make sure that all views contain the same samples in the same order (see documentation)")
  if (length(unique(tmp))>1) 
    stop("Views have different sample names (columns)... please make sure that all views contain the same samples in the same order (see documentation)")
  
  # Make a dgCMatrix out of dgTMatrix
  if (all(sapply(data, function(x) is(x, "dgTMatrix")))) {
    data <- lapply(data, function(m) as(m, "dgCMatrix"))
  }
  
  # Set groups names
  if (is.null(groups)) {
    # message("No groups provided as argument... we assume that all samples are coming from the same group.\n")
    groups <- rep("group1", ncol(data[[1]]))
  }
  
  # Set views names
  if (is.null(names(data))) {
    default_views <- paste0("view_", seq_len(length(data)))
    message(paste0("View names are not specified in the data, using default: ", paste(default_views, collapse=", "), "\n"))
    names(data) <- default_views
  }
  views_names <- as.character(names(data))
  
  # Initialise MOFA object
  object <- new("MOFA")
  object@status <- "untrained"
  object@data <- .split_data_into_groups(data, groups)
  
  # groups_names <- as.character(unique(groups))
  groups_names <- names(object@data[[1]])
  
  # Set dimensionalities
  object@dimensions[["M"]] <- length(data)
  object@dimensions[["G"]] <- length(groups_names)
  object@dimensions[["D"]] <- sapply(data, nrow)
  object@dimensions[["N"]] <- sapply(groups_names, function(x) sum(groups == x))
  object@dimensions[["K"]] <- 0
  
  # Set features names
  for (m in seq_len(length(data))) {
    if (is.null(rownames(data[[m]]))) {
      warning(sprintf("Feature names are not specified for view %d, using default: feature1_v%d, feature2_v%d...", m, m, m))
      for (g in seq_len(length(object@data[[m]]))) {
        rownames(object@data[[m]][[g]]) <- paste0("feature_", seq_len(nrow(object@data[[m]][[g]])), "_v", m)
      }
    }
  }
  
  # Set samples names
  for (g in seq_len(object@dimensions[["G"]])) {
    if (is.null(colnames(object@data[[1]][[g]]))) {
      warning(sprintf("Sample names for group %d are not specified, using default: sample1_g%d, sample2_g%d,...", g, g, g))
      for (m in seq_len(object@dimensions[["M"]])) {
        colnames(object@data[[m]][[g]]) <- paste0("sample_", seq_len(ncol(object@data[[m]][[g]])), "_g", g)
      }
    }
  }
  
  # Set view names
  views_names(object) <- views_names
  
  # Set samples group names
  groups_names(object) <- groups_names
  
  return(object)
}

# (Hidden) function to split a list of matrices into a nested list of matrices
.split_data_into_groups <- function(data, groups) {
  group_indices <- split(seq_along(groups), factor(groups, exclude = character(0))) # factor call avoids dropping NA
  lapply(data, function(x) {
    lapply(group_indices, function(idx) {
      x[, idx, drop = FALSE]
    })
  })
}

# (Hidden) function to split data in Seurat object into a list of matrices
.split_seurat_into_groups <- function(seurat, groups, assay = "RNA", slot = "data", features = NULL) {
  data <- Seurat::GetAssayData(object = seurat, assay = assay, slot = slot)
  if (!is.null(features)) data <- data[features, , drop=FALSE]
  .split_data_into_groups(list(data), groups)[[1]]
}

# (Hidden) function to split data in a SingleCellExperiment object into a list of matrices
.split_sce_into_groups <- function(sce, groups, assay) {
  
  if(!requireNamespace("SummarizedExperiment", quietly = TRUE)){
    stop("Package \"SummarizedExperiment\" is required but is not installed.", call. = FALSE)
  } else {
  
    data <- SummarizedExperiment::assay(sce, assay = assay)
    .split_data_into_groups(list(data), groups)[[1]]
  }
}

# (Hidden) function to fill NAs for missing samples
.subset_augment<-function(mat, samp) {
  samp <- unique(samp)
  mat <- t(mat)
  aug_mat<-matrix(NA, ncol=ncol(mat), nrow=length(samp))
  aug_mat<-mat[match(samp,rownames(mat)),,drop=FALSE]
  rownames(aug_mat)<-samp
  colnames(aug_mat)<-colnames(mat)
  return(t(aug_mat))
}

.df_to_matrix <- function(x) {
  m <- as.matrix(x[,-1])
  rownames(m) <- x[[1]]
  m
}