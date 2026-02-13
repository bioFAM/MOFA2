
#' @title create a MOFA object
#' @name create_mofa
#' @description Method to create a \code{\link{MOFA}} object. Depending on the input data format, this method calls one of the following functions:
#' \itemize{
#'   \item{\strong{long data.frame}: \code{\link{create_mofa_from_df}}}
#'   \item{\strong{List of matrices}: \code{\link{create_mofa_from_matrix}}}
#'   \item{\strong{MultiAssayExperiment}: \code{\link{create_mofa_from_MultiAssayExperiment}}}
#'   \item{\strong{Seurat}: \code{\link{create_mofa_from_Seurat}}}
#'   \item{\strong{SingleCellExperiment}: \code{\link{create_mofa_from_SingleCellExperiment}}}
#'   }
#'  Please read the documentation of the corresponding function for more details on your specific data format.
#' @param data one of the formats above
#' @param groups group information, only relevant when using the multi-group framework. 
#' @param assays assay name(s). Relevant for MultiAssayExperiment, SingleCellExperiment and Seurat objects. For MAE and Seurat objects a vector or list of strings, for SCE objects a single string.
#' @param alt_experiments altExp slots to load, only relevant for SingleCellExperiment objects. Default is \code{NULL}.
#' @param extract_metadata logical indicating whether to incorporate the sample metadata from the input object into the MOFA object (
#' not relevant when the input is a list of matrices). Default is \code{TRUE}.
#' @param ... further arguments that can be passed to the function depending on the input data format.
#' See the documentation of above functions for details.
#' @return Returns an untrained \code{\link{MOFA}} object
#' @export
#' @examples
#' # Using an existing simulated data with two groups and two views
#' file <- system.file("extdata", "test_data.RData", package = "MOFA2")
#' 
#' # Load data (in long data.frame format)
#' load(file) 
#' MOFAmodel <- create_mofa(dt)
create_mofa <- function(data, assays = NULL, groups = NULL, alt_experiments = NULL, extract_metadata = TRUE, ...) {
  
  # Creating MOFA object from a Seurat object
  if (is(data, "Seurat")) {
    
    message("Creating MOFA object from a Seurat object...")
    object <- create_mofa_from_Seurat(data, groups, assays, extract_metadata = extract_metadata, ...)
    
    # Creating MOFA object from a SingleCellExperiment object
  } else if (is(data, "SingleCellExperiment")) {
    
    message("Creating MOFA object from a SingleCellExperiment object...")
    object <- create_mofa_from_SingleCellExperiment(data, assays = assays, groups=groups, alt_experiments = alt_experiments, extract_metadata = extract_metadata, ...)
    
    
    # Creating MOFA object from a data.frame object
  } else if (is(data, "data.frame")) {
    
    message("Creating MOFA object from a data.frame...")
    object <- create_mofa_from_df(data, extract_metadata = extract_metadata)
    
    # Creating MOFA object from a (sparse) matrix object
  } else if (is(data, "list") && (length(data) > 0) && 
             (all(sapply(data, function(x) is(x, "matrix"))) || 
              all(sapply(data, function(x) is(x, "dgCMatrix"))) || 
              all(sapply(data, function(x) is(x, "dgTMatrix"))))) {
    
    message("Creating MOFA object from a list of matrices (features as rows, sample as columns)...\n")
    object <- create_mofa_from_matrix(data, groups)
    
    # Creating MOFA object from MultiAssayExperiment object
  } else if(is(data, "MultiAssayExperiment")){
    
    message("Creating MOFA object from a MultiAssayExperiment...")
    object <- create_mofa_from_MultiAssayExperiment(data, assays, groups, extract_metadata = extract_metadata, ...)
    
  } else {
    stop("Error: input data has to be provided as a list of matrices, a data frame or a Seurat object. Please read the documentation for more details.")
  }
  
  return(object)
}

#' @title create a MOFA object from a MultiAssayExperiment object
#' @name create_mofa_from_MultiAssayExperiment
#' @description Method to create a \code{\link{MOFA}} object from a \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}} object
#' @param mae a \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}} object
#' @param assays assays to use from the \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}} object
#' @param groups a string specifying column name of the colData to use it as a group variable. 
#' Alternatively, a character vector with group assignment for every sample.
#' Default is \code{NULL} (no group structure).
#' @param extract_metadata logical indicating whether to incorporate the metadata from the MultiAssayExperiment object into the MOFA object
#' @return Returns an untrained \code{\link{MOFA}} object
#' @export
#' @examples
#' # Using the miniACC data transform RNASeq assays and select three experiments for MOFA
#' library(MultiAssayExperiment)
#' data(miniACC)
#' # apply log1p to the RNASeq and miRNAseq counts
#' for (nm in c("RNASeq2GeneNorm", "miRNASeqGene")) {
#'   assay(miniACC[[nm]], "log1p") <-
#'     log1p(assay(miniACC[[nm]], "exprs"))
#' }
#' # create model with the transformed assay, omitting the 2nd and 4th experiments
#' model <- create_mofa_from_MultiAssayExperiment(miniACC,
#'                     assays = list("log1p", NULL, "exprs", NULL,'log1p'),
#'                     extract_metadata = TRUE)
create_mofa_from_MultiAssayExperiment <- function(mae, assays = NULL, groups = NULL, extract_metadata = FALSE) {
  
  # Sanity check
  if(!requireNamespace("MultiAssayExperiment", quietly = TRUE)){
    stop("Package \"MultiAssayExperiment\" is required but is not installed.", call. = FALSE)
  } else {
    # Select assays of each experiment for MOFA
    mae <- .select_assays(mae, assays)
    
    # Re-arrange data for training in MOFA to matrices, fill in NAs
    data_list <- lapply(names(mae), function(m) {
      
      # Extract general sample names
      primary <- unique(MultiAssayExperiment::sampleMap(mae)[,"primary"])
      
      # Extract view
      subdata <- as.matrix(MultiAssayExperiment::assays(mae)[[m]])
      
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
    
    # Extract metadata
    if (extract_metadata) {
      if (ncol(MultiAssayExperiment::colData(mae)) > 0) {
        object@samples_metadata <- data.frame(MultiAssayExperiment::colData(mae))
      }
    }

    # Create sample metadata
    object <- .create_samples_metadata(object)

    # Create features metadata
    object <- .create_features_metadata(object)

    # Rename duplicated features
    object <- .rename_duplicated_features(object)

    # Do quality control
    object <- .quality_control(object)
    
    return(object)
  }
}


#' @title create a MOFA object from a data.frame object
#' @name create_mofa_from_df
#' @description Method to create a \code{\link{MOFA}} object from a data.frame object
#' @param df \code{data.frame} object with at most 5 columns: \code{sample}, \code{group}, \code{feature}, \code{view}, \code{value}. 
#'   The \code{group} column (optional) indicates the group of each sample when using the multi-group framework.
#'   The \code{view} column (optional) indicates the view of each feature when having multi-view data.
#' @param extract_metadata  logical indicating whether to incorporate the extra columns as sample metadata into the MOFA object
#' @return Returns an untrained \code{\link{MOFA}} object
#' @export
#' @examples
#' # Using an existing simulated data with two groups and two views
#' file <- system.file("extdata", "test_data.RData", package = "MOFA2")
#' 
#' # Load data (in long data.frame format)
#' load(file) 
#' MOFAmodel <- create_mofa_from_df(dt)
create_mofa_from_df <- function(df, extract_metadata = TRUE) {
  
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
  stopifnot(all(c("sample","feature","value") %in% colnames(df)))
  # stopifnot(all(colnames(df) %in% (c("sample","feature","value","group","view"))))
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
    features <- as.character( unique( df[df$view==m,"feature",drop=TRUE] ) )
    for (g in levels(df$group)) {
      samples <- as.character( unique( df[df$group==g,"sample",drop=TRUE] ) )
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
  
  # save other sample-level columns to samples metadata (e.g. covariates)
  if(extract_metadata && !all(colnames(df) %in% (c("sample","feature","value","group","view")))) {
    cols2keep <- df %>% group_by(sample) %>% select(-c("view", "feature", "value", "group", "value")) %>%
      summarise(across(!starts_with("sample"), function(x) length(unique(x)),
                       .names = "{col}")) 
    cols2keep <- colnames(cols2keep)[apply(cols2keep, 2, function(x) all(x  == 1))]
    if (length(cols2keep) > 0){
      df_meta <- df[, c("sample",cols2keep)] %>% distinct()
      object@samples_metadata <- df_meta %>% select(-sample)
      rownames(object@samples_metadata) <- df_meta$sample
    }
  }

    # Create sample metadata
    object <- .create_samples_metadata(object)

    # Create features metadata
    object <- .create_features_metadata(object)

    # Rename duplicated features
    object <- .rename_duplicated_features(object)

    # Do quality control
    object <- .quality_control(object)

  return(object)
}


#' @title create a MOFA object from a SingleCellExperiment object
#' @name create_mofa_from_SingleCellExperiment
#' @description Method to create a \code{\link{MOFA}} object from a SingleCellExperiment object
#' @param sce SingleCellExperiment object
#' @param assays assays of the sce object to use, default is \code{c('logcounts')}. The first entry specifies the assay of the main experiment, subsequent entries are used for the \code{alt_experiments}. Should match length of \code{alt_experiments} plus one.
#' @param groups a string specifying column name of the colData to use it as a group variable. 
#' Alternatively, a character vector with group assignment for every sample.
#' Default is \code{NULL} (no group structure).
#' @param alt_experiments list, that indicates altExp's of the SingleCellExperiment object to load. Default is \code{NULL} (no altExp's loaded).
#' @param extract_metadata logical indicating whether to incorporate the metadata from the SingleCellExperiment object into the MOFA object
#' @return Returns an untrained \code{\link{MOFA}} object
#' @export
#' @examples
#' # Minimal Example
#' library(SingleCellExperiment)
#' # Create SingleCellExperiment with sample data
#' data1 = matrix(rnorm(50 * 20), 
#'              nrow = 50, ncol = 20,
#'              dimnames = list(paste0("gene1_", 1:50),
#'              paste0("sample_", 1:20)))
#' 
#' data2 = matrix(rnorm(50 * 20), 
#'              nrow = 50, ncol = 20,
#'              dimnames = list(paste0("gene2_", 1:50),
#'              paste0("sample_", 1:20)))
#' alt_data = SummarizedExperiment(list(counts=data2))
#' 
#' sce <- SingleCellExperiment(
#'   assays = list(counts = data1),
#'   altExps = list(alt_data=alt_data)
#' )
#' # Create MOFA object
#' MOFAobject <- create_mofa_from_SingleCellExperiment(
#'   sce, assays = c("counts","counts"), alt_experiments = c("alt_data")
#' )
create_mofa_from_SingleCellExperiment <- function(sce, assays = c("logcounts"), groups = NULL, alt_experiments = NULL, extract_metadata = FALSE) {
  
  # Check is SingleCellExperiment is installed
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package \"SingleCellExperiment\" is required but is not installed.", call. = FALSE)
  }
  else if(!requireNamespace("SummarizedExperiment", quietly = TRUE)){
    stop("Package \"SummarizedExperiment\" is required but is not installed.", call. = FALSE)
  } else {
    # Check that and alt_experiments lists/vectors match in length
    if ((length(assays)-length(alt_experiments))!=1) {
      stop("length of passed `assays` must match `alt_experiments` plus one.")
    }
    assay <- assays[1]
    if (length(assays)==0) {
      alt_assays <- NULL
    } else {
      alt_assays <- assays[-1]
    }

    stopifnot(assay%in%names(SummarizedExperiment::assays(sce)))
    
    # check passed alt_experiments
    if (!is.null(alt_experiments )) {
      if (!(length(alt_experiments) == length(alt_assays))) {
        stop("length of passed `assays` must match `alt_experiments` plus one.")
      }
      stopifnot(all(alt_experiments %in% names(altExps(sce))))
      stopifnot(all(sapply(alt_experiments, function(x) is(altExps(sce)[[x]], "SummarizedExperiment"))))
    }

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
    data_matrices <- .split_sce_into_groups_alt(sce, groups, assay, alt_experiments, alt_assays)
    
    # Create MOFA object
    object <- new("MOFA")
    object@status <- "untrained"
    object@data <- data_matrices
    
    # Define dimensions
    object@dimensions[["M"]] <- 1 + length(alt_experiments)
    object@dimensions[["D"]] <- vapply(data_matrices, function(m) nrow(m[[1]]), 1L)
    object@dimensions[["G"]] <- length(data_matrices[[1]])
    object@dimensions[["N"]] <- vapply(data_matrices[[1]], function(g) ncol(g), 1L)
    object@dimensions[["K"]] <- 0

        
    # Set views & groups names
    groups_names(object) <- as.character(names(data_matrices[[1]]))
    views_names(object)  <- c("Main", alt_experiments)
    
    # Set metadata
    if (extract_metadata) {
      object@samples_metadata <- as.data.frame(colData(sce))
      # object@features_metadata <- as.data.frame(rowData(sce))
    }
    
    # Create sample metadata
    object <- .create_samples_metadata(object)

    # Create features metadata
    object <- .create_features_metadata(object)

    # Rename duplicated features
    object <- .rename_duplicated_features(object)

    # Do quality control
    object <- .quality_control(object)

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
#' @param layer layer to be used (default is data).
#' @param features a list with vectors, which are used to subset features, with names corresponding to assays; a vector can be provided when only one assay is used
#' @param extract_metadata logical indicating whether to incorporate the metadata from the Seurat object into the MOFA object
#' @return Returns an untrained \code{\link{MOFA}} object
#' @export
create_mofa_from_Seurat <- function(seurat, groups = NULL, assays = NULL, layer = "data", features = NULL, extract_metadata = FALSE) {
  
  # Check is Seurat is installed
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package \"Seurat\" is required but is not installed.", call. = FALSE)
  } else {
    
    # Check Seurat version
    if (SeuratObject::Version(seurat)$major != 5) stop("Please install Seurat v5")
    
    # Define assays
    if (is.null(assays)) {
      assays <- SeuratObject::Assays(seurat)
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
        # features <- lapply(assays, function(i) seurat@assays[[i]]@var.features)
        features <- lapply(assays, function(i) SeuratObject::VariableFeatures(seurat@assays[[i]]))
        names(features) <- assays
        if (any(sapply(features,length)==0)) stop("No list of features provided and variable features not detected in the Seurat object")
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
      .split_seurat_into_groups(seurat, groups = groups, assay = i, layer = layer, features = features[[i]]))
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
    if (extract_metadata) {
      object@samples_metadata <- seurat@meta.data
      # object@features_metadata <- do.call(rbind, lapply(assays, function(a) seurat@assays[[a]]@meta.features))
    }

    # Create sample metadata
    object <- .create_samples_metadata(object)

    # Create features metadata
    object <- .create_features_metadata(object)

    # Rename duplicated features
    object <- .rename_duplicated_features(object)

    # Do quality control
    object <- .quality_control(object)
    
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

  # Create sample metadata
  object <- .create_samples_metadata(object)

  # Create features metadata
  object <- .create_features_metadata(object)

  # Rename duplicated features
  object <- .rename_duplicated_features(object)

  # Do quality control
  object <- .quality_control(object)

  return(object)
}



###########################
###### HELP FUNCTIONS #####
###########################

# Select assays from MultiAssayExperiment
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SummarizedExperiment assay assays assayNames
.select_assays <- function(mae, assay_names) {
  # skip if assay_names is NULL
  if(!is.null(assay_names)) {
    # replace `NULL` entries with character(0)
    assay_names = as.list(assay_names)
    assay_names <- lapply(assay_names, function(el) {
      if (is.null(el) || length(el) == 0 || el == "") character(0) else el
    })
    # check that assay_names has right length
    if (length(assay_names) != length(experiments(mae))){
      stop(paste0("Length of assay_names list must match number of experiments (num=",length(experiments(mae)),")"))
    }
    # Give corresponding experiment names to assays
    names(assay_names) <- names(experiments(mae))
    # For every experiment in MAE
    for (exp in names(experiments(mae)) ){
      # If assay.name is `NULL` drop the experiment from MAE
      if ( length(assay_names[[exp]]) == 0 ) {
        # Remove experiment from ExperimentList
        message(paste0("removing experiment: ",exp))
        experiments(mae)[[exp]] <- NULL
      } else {
        # Select given assay if object is SE, otherwise skip
        if (is(mae[[exp]], "SummarizedExperiment")) {
          # check that assay.name is in experiment
          if (! assay_names[[exp]] %in%  assayNames(mae[[exp]])) {
            stop("Cannot find assay '", assay_names[[exp]], "' in experiment '",
            exp, "'.")
          }
          # Keep only specified assay.name from the given experiment
          assays(mae[[exp]]) <- list(assay(mae[[exp]], assay_names[[exp]]))
          # Update assay names (needed?)
          assayNames(mae[[exp]]) <- assay_names[[exp]]
        }
      }
    }
  }
  return(mae)
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
.split_seurat_into_groups <- function(seurat, groups, assay = "RNA", layer = "data", features = NULL) {
  data <- SeuratObject::GetAssayData(object = seurat, assay = assay, layer = layer)
  if(is.null(data) | any(dim(data) == 0)){
    stop(paste("No data present in the layer",layer, "of the assay",assay ,"in the Seurat object."))
  }
  if (!is.null(features)) data <- data[features, , drop=FALSE]
  .split_data_into_groups(list(data), groups)[[1]]
}

# (Hidden) function to select altExps from SingleCellExperiment object and split data into groups
.split_sce_into_groups_alt <-function(sce, groups, assay, alt_experiments, alt_assays) {
    # 1. get main matrix and put it into 1st entry of data_list
    main_data <- SummarizedExperiment::assay(sce, i = assay)
    data_list <- list(main_data)

    if (!is.null(alt_experiments)) {
      # 2. get subsequent altExps and put them into following entries
      for (i in seq_along(alt_experiments)) {
        alt_se <- altExp(sce, alt_experiments[i])
        if(!(alt_assays[i] %in% assayNames(alt_se))){
          stop(paste0("Error: assay `",alt_assays[i],"` not found in altExp `",alt_experiments[i],"`."))
        }
        data <- assays(alt_se)[[alt_assays[i]]]
        data_list[[i+1]] <- data
      }
    }
    # 3. set matrices names
    names(data_list) <- c("main", alt_experiments)
    # 2. spit matrix list into groups with .split_data_into_groups
    data_matrices <- .split_data_into_groups(data_list, groups)
    return(data_matrices)
}

# (Hidden) function to split data in a SingleCellExperiment object into a list of matrices
#.split_sce_into_groups <- function(sce, groups, assay) {
#  
#  if(!requireNamespace("SummarizedExperiment", quietly = TRUE)){
#    stop("Package \"SummarizedExperiment\" is required but is not installed.", call. = FALSE)
#  } else {
#    
#    data <- SummarizedExperiment::assay(sce, i = assay)
#    .split_data_into_groups(list(data), groups)[[1]]
#  }
#}

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
  if (ncol(m) == 1)
    colnames(m) <- colnames(x)[2:ncol(x)]
  m
}

.create_samples_metadata <- function(object) {
  # TO-DO: CHECK SAMPLE AND GROUP COLUMN IN PROVIDED METADATA
  foo <- lapply(object@data[[1]], colnames)
  tmp <- data.frame(
    sample = unname(unlist(foo)),
    group = unlist(lapply(names(foo), function(x) rep(x, length(foo[[x]])) )),
    stringsAsFactors = FALSE
  )
  if (.hasSlot(object, "samples_metadata") && (length(object@samples_metadata) > 0)) {
    object@samples_metadata <- cbind(tmp, object@samples_metadata[match(tmp$sample, rownames(object@samples_metadata)),, drop = FALSE])
  } else {
    object@samples_metadata <- tmp
  }
  return(object)
}

.create_features_metadata <- function(object) {
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
  return(object)
}

.rename_duplicated_features <- function(object) {
  feature_names <- unname(unlist(lapply(object@data, function(x) rownames(x[[1]]))))
  duplicated_names <- unique(feature_names[duplicated(feature_names)])
  if (length(duplicated_names)>0) 
    warning("There are duplicated features names across different views. We will add the suffix *_view* only for those features 
            Example: if you have both TP53 in mRNA and mutation data it will be renamed to TP53_mRNA, TP53_mutation")
  # Rename data matrices
  for (m in names(object@data)) {
    for (g in names(object@data[[m]])) {
      tmp <- which(rownames(object@data[[m]][[g]]) %in% duplicated_names)
      if (length(tmp)>0) {
        rownames(object@data[[m]][[g]])[tmp] <- paste(rownames(object@data[[m]][[g]])[tmp], m, sep="_")
      }
    }
  }
  
  # Rename features metadata
  tmp <- object@features_metadata[["feature"]] %in% duplicated_names
  object@features_metadata[tmp,"feature"] <- paste(object@features_metadata[tmp,"feature"], object@features_metadata[tmp,"view"], sep="_")
  return(object)
}
