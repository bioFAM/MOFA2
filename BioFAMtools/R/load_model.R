
############################################
## Functions to load a trained BioFAModel ##
############################################

#' @title loading a trained BioFAModel
#' @name load_model
#' @description Method to load a trained BioFAModel \cr
#' The training of biofam is done using a Python framework, and the model output is saved as an .hdf5 file, which has to be loaded in the R package.
#' @param file an hdf5 file saved by the biofam Python framework
#' @param object either NULL (default) or an an existing untrained biofam object. If NULL, the \code{\link{BioFAModel}} object is created from the scratch.
#' @param sort_factors boolean indicating whether factors should be sorted by variance explained (default is TRUE)
#' @param on_disk boolean indicating whether to work from memory (FALSE) or disk (TRUE). \cr
#' This should be set to TRUE when the training data is so big that cannot fit into memory. \cr
#' On-disk operations are performed using the \code{\link{HDF5Array}} and \code{\link{DelayedArray}} framework.
#' @param load_training_data boolean indicating whether to load the training data (default is TRUE, it can be memory expensive)
#' @return a \code{\link{BioFAModel}} model
#' @importFrom rhdf5 h5read h5ls
#' @importFrom HDF5Array HDF5ArraySeed
#' @importFrom DelayedArray DelayedArray
#' @export

load_model <- function(file, object = NULL, sort_factors = TRUE, on_disk = FALSE, load_training_data = TRUE) {

  # Create new bioFAModel object
  if (is.null(object)) object <- new("BioFAModel")
  object@status <- "trained"
  
  # Set on_disk option
  if (.hasSlot(object, "on_disk")) {
    if (on_disk) { object@on_disk <- TRUE } else { object@on_disk <- FALSE }
  }
  
  # Get groups and data set names from the hdf5 file object
  foo <- h5ls(file, datasetinfo = F)

  ########################
  ## Load training data ##
  ########################

  # Load identity of features and samples
  feature_names <- h5read(file, "features")
  sample_names  <- h5read(file, "samples")
  view_names <- foo[foo$group=="/data","name"]
  group_names <- foo[foo$group==paste0("/data/",view_names[1]),"name"]

  # Load training data as matrices
  training_data <- list()
  if (load_training_data) {
    for (m in view_names) {
      training_data[[m]] <- list()
      for (p in group_names) {
        if (on_disk) {
          # as DelayedArrays
          training_data[[m]][[p]] <- DelayedArray( HDF5ArraySeed(file, name = sprintf("data/%s/%s", m, p) ) )
        } else {
          # as matrices
          training_data[[m]][[p]] <- h5read(file, sprintf("data/%s/%s", m, p) )
        }
      }
    }
  } else {
    n_features <- lapply(feature_names, length)
    n_samples  <- lapply(sample_names, length)
    for (m in view_names) {
      training_data[[m]] <- list()
      for (p in group_names) {
        training_data[[m]][[p]] <- .create_matrix_placeholder(rownames = feature_names[[m]], colnames = sample_names[[p]])
      }
    }
  }

  # Replace NaN by NA
  # RICARD: IF USING ON_DISK, I THINK THIS REALISES EVERYTHING INTO MEMORY, TO CHECK
  for (m in view_names) {
    for (p in group_names) {
      training_data[[m]][[p]][is.nan(training_data[[m]][[p]])] <- NA
    }
  }

  # Give corresponding names for rows (features) and columns (samples)
  # RICARD: I THINK THIS REALISES EVERYTHING INTO MEMORY, TO CHECK
  # tryCatch( {
  #   for (m in 1:length(training_data)) {
  #     for (p in 1:length(training_data[[m]])) {
  #       rownames(training_data[[m]][[p]]) <- feature_names[[m]]
  #       colnames(training_data[[m]][[p]]) <- sample_names[[p]]
  #     }
  #   }
  #   object@training_data <- training_data
  # }, error = function(x) { cat("Error defining feature and sample names\n") })

  object@training_data <- training_data

  #######################
  ## Load expectations ##
  #######################

  expectations <- list()
  node_names <- foo[foo$group=="/expectations","name"]

  if ("AlphaW" %in% node_names)
    expectations[["AlphaW"]] <- h5read(file, "expectations/AlphaW")
  if ("AlphaZ" %in% node_names)
    expectations[["AlphaZ"]] <- h5read(file, "expectations/AlphaZ")
  if ("Z" %in% node_names)
    expectations[["Z"]] <- h5read(file, "expectations/Z")
  if ("W" %in% node_names)
    expectations[["W"]] <- h5read(file, "expectations/W")
  if ("ThetaW" %in% node_names)
    expectations[["ThetaW"]] <- h5read(file, "expectations/ThetaW")
  if ("ThetaZ" %in% node_names)
    expectations[["ThetaZ"]] <- h5read(file, "expectations/ThetaZ")
  if ("Tau" %in% node_names)
    expectations[["Tau"]] <- h5read(file, "expectations/Tau") # TO-DO: DELAYEDARRAY
  
  if ("Y" %in% node_names) {
    expectations[["Y"]] <- list()
    for (m in view_names) {
      expectations[["Y"]][[m]] <- list()
      for (p in group_names) {
        if (on_disk) {
          expectations[["Y"]][[m]][[p]] <- DelayedArray( HDF5ArraySeed(file, name=sprintf("expectations/Y/%s/%s", m, p)) )
        } else {
          expectations[["Y"]][[m]][[p]] <- h5read(file, sprintf("expectations/Y/%s/%s", m, p))
        }
      }
    }
  }

  object@expectations <- expectations

  ##############################
  ## Specify dimensionalities ##
  ##############################

  # Specify dimensionality of the data
  object@dimensions[["M"]] <- length(training_data)                           # number of views (groups of features)
  object@dimensions[["P"]] <- length(training_data[[1]])                      # number of groups (groups of samples)
  object@dimensions[["N"]] <- sapply(training_data[[1]], ncol)                # number of samples per sample_group
  object@dimensions[["D"]] <- sapply(training_data, function(e) nrow(e[[1]])) # number of features per feature_group (view)
  object@dimensions[["K"]] <- ncol(object@expectations$Z[[1]])                # number of factors


  # Create default samples names if they are null
  if (is.null(sample_names)) {
    print("Samples names not found, generating default: sample1, ..., sampleN")
    sample_names <- lapply(object@dimensions[["N"]], function(n) paste0("sample", as.character(1:n)))
  }
  samples_names(object)  <- sample_names
  
  # Create default features names if they are null
  if (is.null(feature_names)) {
    print("Features names not found, generating default: feature1_view1, ..., featureD_viewM")
    feature_names <- lapply(1:object@dimensions[["M"]],
      function(m) sprintf("feature%d_view_&d", as.character(1:object@dimensions[["D"]][m]), m))
  }
  features_names(object) <- feature_names

  # Set views names
  if (is.null(names(object@training_data))) {
    print("Views names not found, generating default: view1, ..., viewM")
    view_names <- paste0("view", as.character(1:object@dimensions[["M"]]))
  }
  views_names(object) <- view_names

  # Set groups names
  if (is.null(names(object@training_data[[1]]))) {
    print("Groups names not found, generating default: group1, ..., groupG")
    group_names <- paste0("group", as.character(1:object@dimensions[["P"]]))
  }
  groups_names(object) <- group_names

  # Set factors names
  factors_names(object)  <- paste0("Factor", as.character(1:object@dimensions[["K"]]))
  
  ########################
  ## Load model options ##
  ########################

  tryCatch( {
    object@model_options <- as.list(h5read(file, 'model_options', read.attributes=T))
  }, error = function(x) { print("Model options not found, not loading it...") })

  # Convert True/FalsesStrings to logical values
  for (opt in names(object@model_options)) {
    if (object@model_options[opt] == "False" | object@model_options[opt] == "True") {
      object@model_options[opt] <- as.logical(object@model_options[opt])
    } else {
      object@model_options[opt] <- object@model_options[opt]
    }
  }

  # Add views names to likelihood vector
  if (is.null(names(object@model_options$likelihood))) {
    names(object@model_options$likelihood) <- views_names(object)
  }

  ##########################################
  ## Load training options and statistics ##
  ##########################################

  # Load training options
  if (length(object@training_options) == 0) {
    tryCatch( {
      object@training_options <- as.list(h5read(file, 'training_opts', read.attributes=T))
    }, error = function(x) { print("Training opts not found, not loading it...") })
  }

  # Load training statistics
  tryCatch( {
    object@training_stats <- h5read(file, 'training_stats', read.attributes=T)
    colnames(object@training_stats$elbo_terms) <- attr(h5read(file,"training_stats/elbo_terms", read.attributes=T),"colnames")
  }, error = function(x) { print("Training stats not found, not loading it...") })

  ###################
  ## Parse factors ##
  ###################

  # Order factors in order of variance explained
  if (sort_factors) {
    object@cache[["variance_explained"]] <- calculate_variance_explained(object)
    r2 <- rowSums(sapply(object@cache[["variance_explained"]]$r2_per_factor, function(e) rowSums(e)))
    order_factors <- c(names(r2)[order(r2, decreasing = T)])
    object <- subset_factors(object, order_factors)
  }

  # Mask passenger samples
  object <- detect_passengers(object)

  ############################
  ## Update previous models ##
  ############################

  # object <- .update_old_model(object)

  ######################
  ## Quality controls ##
  ######################

  # qualityControl(object)

  return(object)
}
