set_covariates <- function(object, covariates = NULL) {

  # Sanity checks
  if (!is(object, "MOFA")) 
    stop("'object' has to be an instance of MOFA")
  if (object@status=="trained") 
    stop("The model is already trained! Covariates must be added before training")
  
  # get sample names
  samples_data <- lapply(object@data[[1]], colnames)
  # samples <- unlist(samples_data)
  samples_data_vec <- unlist(samples_names(object))
  
  # covariates passed as characters: extract from the metadata
  if (is(covariates, "character")) {
    if (!all(covariates %in% colnames(samples_metadata(object)))) {
      stop("Columns specified in covariates do not exist in the MOFA object metadata slot.")
    }
    covariates <- samples_metadata(object)[,c("sample",covariates),drop=F]
    # TO-DO: Check that they are numeric and continuous
    # Add covariates to the MOFA object
    # object <- .add_covariate(object, covariate)
  
  # covariates passed in data.frame format
  } else if (all(class(covariates) %in% c("data.frame", "tibble", "Data.Frame"))) { # TO-DO: USE is()
      if (!all(c("sample", "covariates", "value") %in% colnames(covariates)))
        stop("If covariates is provided as data.frame it needs to contain the columns: sample, covariate, value")
      if (!is.numeric(covariates$value)) {
        stop("Values in covariates need to be numeric")
      }
      samples <- covariates$sample
      # covariates <- covariates[!duplicated(covariates), ]
      covariates <- reshape2::acast(covariates, covariate ~ sample)
      
  # covariates passed in matrix format
  # TO-DO: CHECK THIS
  } else if (all(is.numeric(covariates)) || class(covariates) %in% c("dgTMatrix", "dgCMatrix")) {
    samples <- colnames(covariates)
    if (!is.null(samples)) {
      if(!(all(samples %in% samples_data_vec) & all(samples_data_vec %in% samples)))
        stop("Sample names of the data and the sample covariates do not match.")
      covariates <- covariates[ , samples_data_vec, drop = FALSE]
    } else {
      # warnings and checks if no matching sample names
      if(sum(object@dimensions[['N']]) != ncol(covariates))
        stop("Number of columns in sample covariates does not match the number of samples")
      if(!is.null(samples_data) & length(samples_data_vec) > 0) {
        warning("No sample names in covariates - we will use the sample names in data. Please ensure that the order matches.")
        colnames(covariates) <- samples
      } else {
        stop("No sample names found!")
      }
    }
    
  # covariates format not recognised
  } else {
    stop("covariates needs to be a character vector, a dataframe, a matrix or NULL.")
  }
    
  # Set covariate dimensionality
  object@dimensions[["C"]] <- nrow(covariates)
    
  # Set covariate names
  if (is.null(rownames(covariates))) {
    message("No covariates names provided - using generic: covariate1, covariate2, ...")
    rownames(covariates) <- paste0("covariate", seq_len(nrow(covariates)))
  }
  
  # split covariates by groups
  covariates <- lapply(samples_names(object), function(i)   covariates[, i, drop = FALSE])
  names(covariates) <- groups_names(object)
  
  # Sanity checks
  stopifnot(all(sapply(object@covariates, ncol) == object@dimensions[["N"]]))
  
  # add covariates to the MOFA object
  object@covariates <- covariates
    
  return(object)
}

