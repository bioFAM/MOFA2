
#' @title create a BioFAM object
#' @name create_biofam
#' @description Method to create a \code{\link{BioFAModel}} object
#' @param data TO SPECIFY
#' @return Returns an untrained \code{\link{BioFAModel}} object
#' @import BiocGenerics
#' @export
create_biofam <- function(data, samples_groups = NULL) {
  
  if (is(data, "MultiAssayExperiment")) {
    stop("Not functional")

  } else if (is(data, "seurat")) {
    message("Creating BioFAM object from a Seurat object...")
    
    object <- .create_biofam_from_seurat(data, samples_groups)
    
  } else if (is(data, "data.frame")) {
    message("Creating BioFAM object from a dataframe...")
    
    object <- .create_biofam_from_df(data)
    
    # Set dimensionalities
    object@dimensions[["M"]] <- length(unique(data$feature_group))
    object@dimensions[["D"]] <- sapply(unique(data$feature_group), function(m) length(unique(data[data$feature_group==m,]$feature)))
    object@dimensions[["P"]] <- length(unique(data$sample_group))
    object@dimensions[["N"]] <- sapply(unique(data$sample_group), function(p) length(unique(data[data$sample_group==p,]$sample)))
    object@dimensions[["K"]] <- 0
    
    # Set view names
    views_names(object) <- as.character(unique(data$feature_group))
    
    # Set sample group names
    groups_names(object) <- as.character(unique(data$sample_group))
    
  } else if (is(data, "list") && (length(data) > 0) && 
             (all(sapply(data, function(x) is(x, "matrix"))) || 
              all(sapply(data, function(x) is(x, "dgCMatrix"))) || 
              all(sapply(data, function(x) is(x, "dgTMatrix"))))) {
    
    message("Creating BioFAM object from a list of matrices... make sure that features are stored in columns and samples in rows\n")
    
    # Quality controls
    stopifnot(all(sapply(data, function(g) all(is.numeric(g)))) || all(sapply(data, function(x) class(x) %in% c("dgTMatrix", "dgCMatrix"))))

    # Make a dgCMatrix out of dgTMatrix
    if (all(sapply(data, function(x) is(x, "dgTMatrix")))) {
      data <- lapply(data, function(m) as(m, "dgCMatrix"))
    }
    
    # Set samples groups
    if (is.null(samples_groups)) {
      message("No samples_groups provided as argument... we assume that all samples are coming from the same group.\n")
      samples_groups <- rep("group1", nrow(data[[1]]))
    }
    groups_names <- as.character(unique(samples_groups))
    
    # Set views names
    if (is.null(names(data))) {
      default_views_names <- paste0("view_", 1:length(data))
      message(paste0("View names are not specified in the data, using default: ", paste(default_views_names, collapse=", "), "\n"))
      names(data) <- default_views_names
    }
    
    # Initialise BioFAM object
    object <- new("BioFAModel")
    object@status <- "untrained" # define status as untrained
    object@input_data <- .split_data_into_groups(data, samples_groups) # pass input data
    
    # Set dimensionalities
    object@dimensions[["M"]] <- length(data)
    object@dimensions[["P"]] <- length(groups_names)
    object@dimensions[["D"]] <- sapply(data, function(m) ncol(m))
    object@dimensions[["N"]] <- sapply(groups_names, function(x) sum(samples_groups == x))
    object@dimensions[["K"]] <- 0

    # Set features names
    for (m in 1:length(data)) {
      if (is.null(colnames(data[[m]]))) {
        warning(sprintf("Feature names are not specified for view %d, using default: feature1_v%d, feature2_v%d...", m, m, m))
        for (g in 1:length(object@input_data[[m]])) {
          colnames(object@input_data[[m]][[g]]) <- paste0("feature_", 1:ncol(object@input_data[[m]][[g]]), "_v", m)
        }
      }
    }
    
    # Set samples names
    for (g in 1:object@dimensions[["P"]]) {
      if (is.null(rownames(object@input_data[[1]][[g]]))) {
        warning(sprintf("Sample names for group %d are not specified, using default: sample1_g%d, sample2_g%d,...", g, g, g))
        for (m in 1:object@dimensions[["M"]]) {
          rownames(object@input_data[[m]][[g]]) <- paste0("sample_", 1:nrow(object@input_data[[m]][[g]]), "_g", g)
        }
      }
    }
    
    # Set view names
    views_names(object) <- names(object@input_data)

    # Set samples group names
    groups_names(object) <- names(object@input_data[[1]])
    
  } else {
    stop("Error: input data has to be provided as a list of matrices, a data frame (long format), or a suerat object.")
  }
  
  print(object)
  
  return(object)
}



# (Hidden) function to initialise a BioFAModel object using a Dataframe
.create_biofam_from_df <- function(df) {
  
  # Quality controls
  stopifnot(all(colnames(df) %in% (c("sample","feature","value","sample_group","feature_group"))))
  stopifnot(all(is.numeric(df$value)))
  
  # Convert 'sample' and'feature' columns to factors
  if (!is.factor(df$sample))
    df$sample <- as.factor(df$sample)
  if (!is.factor(df$feature))
    df$feature <- as.factor(df$feature)
  
  # Convert 'sample_group' columns to factors
  if (!"sample_group" %in% colnames(df)) {
    df$sample_group <- factor("sample_group1")
  } else {
    df$sample_group <- factor(df$sample_group)
  }
  
  # Convert 'feature_group' columns to factors
  if (!"feature_group" %in% colnames(df)) {
    df$feature_group <- factor("feature_group1")
  } else {
    df$feature_group <- factor(df$feature_group)
  }
  
  # TO-FINISH.....
  
  # samples <- as.character(unique(df$sample))
  # data_matrix <- lapply(split(df,df$feature_group),
  #   function(x) lapply(split(x,x$sample_group),
  #     function(y) {
  #       y <- .df_to_matrix( reshape2::dcast(y, sample~feature, value.var="value", fill=NA, drop=TRUE))
  #       y <- .subset_augment(y,samples)
  #       y
  #     }
  # ))
  
  # data_matrix <- lapply(split(df,df$feature_group), function(x) 
  #   lapply(split(droplevels.data.frame(x),x$sample_group), function(y) {
  #       y <- droplevels.data.frame(y)
  #       if (nrow(y)==0) {
  #         # print(length(levels(y$sample)))
  #         matrix(NA, nrow=length(levels(y$sample)), ncol=length(levels(y$feature)))
  #       } else {
  #         tmp <- .df_to_matrix( reshape2::dcast(y, sample~feature, value.var="value", fill=NA, drop=TRUE) )
  #       }
  #     }
  #   )
  # )
  
  # data_matrix <- reshape2::dcast(df, sample~feature, value.var="value", fill=NA, drop=FALSE)
  
  levels(df$feature_group)
  object <- new("BioFAModel")
  object@status <- "untrained" # define status as untrained
  object@input_data <- data_matrix
  
  return(object)
}

.create_biofam_from_seurat <- function(srt, samples_groups) {
  if (is(samples_groups, 'character') && (length(samples_groups) == 1)) {
    if (!(samples_groups %in% colnames(srt@meta.data)))
      stop(paste0(samples_groups, " is not found in the Seurat@meta.data.\n",
                  "If you want to use samples_groups information from Seurat@meta.data,\n",
                  "please ensure to provide a column name that exists. The columns of meta data are:\n",
                  paste0(colnames(srt@meta.data), sep = ", ")))
    samples_groups <- srt@meta.data[,samples_groups]
  }

  if (is.null(samples_groups)) {
    message("No samples_groups provided as argument... we assume that all samples are coming from the same group.\n")
    samples_groups <- rep("group1", ncol(srt@data))
  }

  data_matrices <- list("rna" = .split_seurat_into_groups(srt, samples_groups))

  object <- new("BioFAModel")
  object@status <- "untrained"
  object@input_data <- data_matrices
  
  # Define dimensions
  object@dimensions[["M"]] <- 1
  object@dimensions[["D"]] <- ncol(data_matrices[[1]][[1]])
  object@dimensions[["P"]] <- length(data_matrices[[1]])
  object@dimensions[["N"]] <- sapply(data_matrices[[1]], function(m) nrow(m))
  object@dimensions[["K"]] <- 0

  # Set views & groups names
  groups_names(object) <- as.character(names(data_matrices[[1]]))
  views_names(object)  <- c("rna")

  return(object)
}


.split_data_into_groups <- function(data, samples_groups) {
  groups_names <- unique(samples_groups)
  tmp <- lapply(data, function(view) {
    tmp_view <- lapply(groups_names, function(p) {
      view[samples_groups == p,]
    })
    names(tmp_view) <- groups_names
    tmp_view
  })
  names(tmp) <- names(data)
  tmp
}

.split_seurat_into_groups <- function(srt, samples_groups) {
  groups_names <- unique(samples_groups)
  tmp <- lapply(groups_names, function(g) {
    # If group name is NA, it has to be treated separately
    # due to the way R handles NAs and equal signs
    if (is.na(g)) {
      t(srt@data[,is.na(samples_groups)])
    } else {
      BiocGenerics::t(srt@data[,which(samples_groups == g)])
    }
  })
  names(tmp) <- groups_names
  tmp
}

.df_to_matrix <- function(x) {
  m <- as.matrix(x[,-1])
  rownames(m) <- x[[1]]
  m
}

# (Hidden) function to fill NAs for missing samples
.subset_augment<-function(mat, samples) {
  aug_mat <- matrix(NA, ncol=ncol(mat), nrow=length(samples))
  aug_mat <- mat[match(samples,rownames(mat)),,drop=FALSE]
  rownames(aug_mat) <- samples
  colnames(aug_mat) <- colnames(mat)
  return(aug_mat)
}


