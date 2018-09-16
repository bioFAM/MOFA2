
#' @title create a BioFAM object
#' @name create_biofam
#' @description Method to create a \code{\link{BioFAModel}} object
#' @param data TO SPECIFY
#' @return Returns an untrained \code{\link{BioFAModel}} object
#' @export
create_biofam <- function(data, samples_groups = NULL) {
  
  if (is(data, "MultiAssayExperiment")) {stop("Not functional")
    # message("Creating BioFAM object from a MultiAssayExperiment object...")
    # object <- .createBioFAMobjectFromMAE(data)
    
  } else if (is(data, "SummarizedExperiment")) {
    stop("Not functional")
    
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
      
  } else if (is(data, "list") && 
             (length(data) > 0) && 
             (all(sapply(data, function(x) is(x, "matrix"))) || 
              all(sapply(data, function(x) is(x, "dgCMatrix"))) || 
              all(sapply(data, function(x) is(x, "dgTMatrix"))))) {

      message("Creating BioFAM object from a list of matrices...")

      # Quality controls
      stopifnot(all(sapply(data, function(p) all(is.numeric(p)))))

      if (is.null(samples_groups)) {
        warning("When providing a list of matrices, one matrix per view, samples annotation (groups) can be provided via samples_groups argument.")
        warning("No samples_groups provided. Considering samples are coming from one group.")
        samples_groups <- rep("group1", nrow(data[0]))
      }
    
      object <- .create_biofam(.split_into_groups(data, samples_groups))
      groups_names <- as.character(unique(samples_groups))
      
      # Set dimensionalities
      object@dimensions[["M"]] <- length(data)
      object@dimensions[["P"]] <- length(groups_names)
      object@dimensions[["D"]] <- sapply(data, function(m) ncol(m))
      object@dimensions[["N"]] <- sapply(groups_names, function(x) sum(samples_groups == x))
      object@dimensions[["K"]] <- 0
      
      # Set view names
      views_names(object) <- as.character(names(data))
      
      # Set sample group names
      groups_names(object) <- groups_names
      
  } else {
    stop("Error: input data has to be provided either as .....")
  }
  
  print(object)
  
  return(object)
}



# (Hidden) function to initialise a BioFAModel object
.create_biofam <- function(data) {
  
  # Initialise BioFAM object
  object <- new("BioFAModel")
  object@status <- "untrained"
  object@input_data <- data
  
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
    df$sample_group <- as.factor("sample_group1")
  } else {
    df$sample_group <- as.factor(df$sample_group)
  }
  
  # Convert 'feature_group' columns to factors
  if (!"feature_group" %in% colnames(df)) {
    df$feature_group <- as.factor("feature_group1")
  } else {
    df$feature_group <- as.factor(df$feature_group)
  }
  
  data_matrix <- lapply(split(df,df$feature_group), 
    function(x) lapply(split(x,x$sample_group),
      function(y) .matrix.please( reshape2::dcast(y, sample~feature, value.var="value", fill=NA)
  )))
  object <- .create_biofam(data_matrix)
  
  return(object)
}


.split_into_groups <- function(data, samples_groups) {
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

.matrix.please <- function(x) {
  m <- as.matrix(x[,-1])
  rownames(m) <- x[[1]]
  m
}
