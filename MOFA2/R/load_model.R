
############################################
## Functions to load a trained MOFA ##
############################################

#' @title Load a trained MOFA
#' @name load_model
#' @description Method to load a trained MOFA \cr
#' The training of mofa is done using a Python framework, and the model output is saved as an .hdf5 file, which has to be loaded in the R package.
#' @param file an hdf5 file saved by the mofa Python framework
#' @param sort_factors logical indicating whether factors should be sorted by variance explained (default is TRUE)
#' @param on_disk logical indicating whether to work from memory (FALSE) or disk (TRUE). \cr
#' This should be set to TRUE when the training data is so big that cannot fit into memory. \cr
#' On-disk operations are performed using the \code{\link{HDF5Array}} and \code{\link{DelayedArray}} framework.
#' @param load_data logical indicating whether to load the training data (default is TRUE, it can be memory expensive)
#' @param remove_outliers logical indicating whether to mask outlier values.
#' @return a \code{\link{MOFA}} model
#' @importFrom rhdf5 h5read h5ls
#' @importFrom HDF5Array HDF5ArraySeed
#' @importFrom DelayedArray DelayedArray
#' @export

load_model <- function(file, sort_factors = TRUE, 
                       on_disk = FALSE, load_data = TRUE, load_imputed_data = FALSE, 
                       remove_outliers = FALSE, verbose = FALSE) {

  # Create new MOFAodel object
  object <- new("MOFA")
  object@status <- "trained"
  
  # Set on_disk option
  if (.hasSlot(object, "on_disk")) {
    if (on_disk) { object@on_disk <- TRUE } else { object@on_disk <- FALSE }
  }
  
  # Get groups and data set names from the hdf5 file object
  h5ls.out <- h5ls(file, datasetinfo = FALSE)

  ########################
  ## Load training data ##
  ########################

  # Load identity of features and samples
  feature_names <- h5read(file, "features")
  sample_names  <- h5read(file, "samples")
  view_names <- names(feature_names)
  group_names <- names(sample_names)

  # Load training data (as nested list of matrices)
  if (isTRUE(verbose)) message("Loading data...")
  data <- list()
  intercepts <- list()
  if (load_data) {
    for (m in view_names) {
      data[[m]] <- list()
      intercepts[[m]] <- list()
      for (g in group_names) {
        if (on_disk) {
          # as DelayedArrays
          data[[m]][[g]] <- DelayedArray( HDF5ArraySeed(file, name = sprintf("data/%s/%s", m, g) ) )
        } else {
          # as matrices
          data[[m]][[g]] <- h5read(file, sprintf("data/%s/%s", m, g) )
          tryCatch(intercepts[[m]][[g]] <- as.numeric( h5read(file, sprintf("intercepts/%s/%s", m, g) ) ), error = function(e) { NULL })
        }
        # Replace NaN by NA
        data[[m]][[g]][is.nan(data[[m]][[g]])] <- NA # this realised into memory, TO FIX
      }
    }
    
  # Create empty training data (as nested list of empty matrices, with the correct dimensions)
  } else {
    for (m in view_names) {
      data[[m]] <- list()
      for (g in group_names) {
        data[[m]][[g]] <- .create_matrix_placeholder(rownames = feature_names[[m]], colnames = sample_names[[g]])
      }
    }
  }

  # Give corresponding names for rows (features) and columns (samples)
  # RICARD: I THINK THIS REALISES EVERYTHING INTO MEMORY, TO CHECK
  # tryCatch( {
  #   for (m in 1:length(data)) {
  #     for (p in 1:length(data[[m]])) {
  #       rownames(data[[m]][[p]]) <- feature_names[[m]]
  #       colnames(data[[m]][[p]]) <- sample_names[[p]]
  #     }
  #   }
  #   object@data <- data
  # }, error = function(x) { cat("Error defining feature and sample names\n") })
  
  object@data <- data
  object@intercepts <- intercepts
  
  #######################
  ## Load imputed data ##
  #######################
  
  imputed_data <- list()
  if (load_imputed_data) {
    if (isTRUE(verbose)) message("Loading imputed data...")
    for (m in view_names) {
      imputed_data[[m]] <- list()
      for (g in group_names) {
        imputed_data[[m]][[g]] <- list()
        if (on_disk) {
          # as DelayedArrays
          # imputed_data[[m]][[g]] <- DelayedArray( HDF5ArraySeed(file, name = sprintf("imputed_data/%s/%s", m, g) ) )
        } else {
          # as matrices
          imputed_data[[m]][[g]][["mean"]] <- h5read(file, sprintf("imputed_data/%s/%s/mean", m, g) )
          imputed_data[[m]][[g]][["variance"]] <- h5read(file, sprintf("imputed_data/%s/%s/variance", m, g) )
        }
        # Replace NaN by NA
        # imputed_data[[m]][[g]][is.nan(imputed_data[[m]][[g]])] <- NA # this realises into memory, TO FIX
      }
    }
  }
  object@imputed_data <- imputed_data
  
  #######################
  ## Load expectations ##
  #######################

  expectations <- list()
  node_names <- h5ls.out[h5ls.out$group=="/expectations","name"]

  if (isTRUE(verbose)) message(paste0("Loading expectations for ", length(node_names), " nodes..."))

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
    expectations[["Tau"]] <- h5read(file, "expectations/Tau")
  
  tmp <- as.list(h5read(file, 'model_options', read.attributes = TRUE))[["likelihoods"]]
  names(tmp) <- names(data)
  
  # Load expectations of Y
  # expectations[["Y"]] <- list()
  # for (m in view_names) {
  #   expectations[["Y"]][[m]] <- list()
  #   for (p in group_names) {
  #     if (on_disk) {
  #       if (tmp[[m]]=="gaussian") {
  #         expectations[["Y"]][[m]][[p]] <- data[[m]][[p]]
  #       } else {
  #         expectations[["Y"]][[m]][[p]] <- DelayedArray( HDF5ArraySeed(file, name=sprintf("expectations/Y/%s/%s", m, p)) )
  #       }
  #     } else {
  #       if (tmp[[m]]=="gaussian") {
  #         expectations[["Y"]][[m]][[p]] <- data[[m]][[p]]
  #       } else {
  #         expectations[["Y"]][[m]][[p]] <- h5read(file, sprintf("expectations/Y/%s/%s", m, p))
  #       }
  #     }
  #   }
  # }
  object@expectations <- expectations

  
  ########################
  ## Load model options ##
  ########################

  if (isTRUE(verbose)) message("Loading model options...")

  tryCatch( {
    object@model_options <- as.list(h5read(file, 'model_options', read.attributes = TRUE))
  }, error = function(x) { print("Model options not found, not loading it...") })

  # Convert True/False strings to logical values
  for (opt in names(object@model_options)) {
    if (object@model_options[opt] == "False" | object@model_options[opt] == "True") {
      object@model_options[opt] <- as.logical(object@model_options[opt])
    } else {
      object@model_options[opt] <- object@model_options[opt]
    }
  }

  # Add views names to likelihood vector
  # if (is.null(names(object@model_options$likelihoods))) {
  #   names(object@model_options$likelihoods) <- views(object)
  # }

  ##########################################
  ## Load training options and statistics ##
  ##########################################

  if (isTRUE(verbose)) message("Loading training options and statistics...")

  # Load training options
  if (length(object@training_options) == 0) {
    tryCatch( {
      object@training_options <- as.list(h5read(file, 'training_opts', read.attributes = TRUE))
    }, error = function(x) { print("Training opts not found, not loading it...") })
  }

  # Load training statistics
  tryCatch( {
    object@training_stats <- h5read(file, 'training_stats', read.attributes = TRUE)
  }, error = function(x) { print("Training stats not found, not loading it...") })

  #######################################
  ## Load variance explained estimates ##
  #######################################
  
  if ("variance_explained"%in%h5ls.out$name) {
    r2_list <- list(
      r2_total = h5read(file, "variance_explained/r2_total"),
      r2_per_factor = h5read(file, "variance_explained/r2_per_factor")
    )
    object@cache[["variance_explained"]] <- r2_list
  }
  
  ##############################
  ## Specify dimensionalities ##
  ##############################
  
  # Specify dimensionality of the data
  object@dimensions[["M"]] <- length(data)                            # number of views
  object@dimensions[["G"]] <- length(data[[1]])                       # number of groups
  object@dimensions[["N"]] <- sapply(data[[1]], ncol)                 # number of samples (per group)
  object@dimensions[["D"]] <- sapply(data, function(e) nrow(e[[1]]))  # number of features (per view)
  object@dimensions[["K"]] <- ncol(object@expectations$Z[[1]])        # number of factors
  
  # Assign sample and feature names (slow for large matrices)
  
  # Create default features names if they are null
  if (is.null(feature_names)) {
    print("Features names not found, generating default: feature1_view1, ..., featureD_viewM")
    feature_names <- lapply(seq_len(object@dimensions[["M"]]),
                            function(m) sprintf("feature%d_view_&d", as.character(seq_len(object@dimensions[["D"]][m])), m))
  } else {
    # Check duplicated features names
    all_names <- unname(unlist(feature_names))
    duplicated_names <- unique(all_names[duplicated(all_names)])
    if (length(duplicated_names)>0) 
      warning("There are duplicated features names across different views. We will add the suffix *_view* only for those features 
            Example: if you have both TP53 in mRNA and mutation data it will be renamed to TP53_mRNA, TP53_mutation")
    for (m in names(feature_names)) {
      tmp <- which(feature_names[[m]] %in% duplicated_names)
      if (length(tmp)>0) feature_names[[m]][tmp] <- paste(feature_names[[m]][tmp], m, sep="_")
    }
  }
  features(object) <- feature_names
  
  # Create default samples names if they are null
  if (is.null(sample_names)) {
    print("Samples names not found, generating default: sample1, ..., sampleN")
    sample_names <- lapply(object@dimensions[["N"]], function(n) paste0("sample", as.character(seq_len(n))))
  }
  samples(object) <- sample_names
  
  # Set views names
  if (is.null(names(object@data))) {
    print("Views names not found, generating default: view1, ..., viewM")
    view_names <- paste0("view", as.character(seq_len(object@dimensions[["M"]])))
  }
  views(object) <- view_names
  
  # Set groups names
  if (is.null(names(object@data[[1]]))) {
    print("Groups names not found, generating default: group1, ..., groupG")
    group_names <- paste0("group", as.character(seq_len(object@dimensions[["G"]])))
  }
  groups(object) <- group_names
  
  # Set factors names
  factors(object)  <- paste0("Factor", as.character(seq_len(object@dimensions[["K"]])))
  
  ###################
  ## Parse factors ##
  ###################

  # Order factors in order of variance explained
  if (sort_factors) {
    if (isTRUE(verbose)) message("Re-ordering factors by their variance explained...")
    if (is.null(object@cache[["variance_explained"]])) 
      object@cache[["variance_explained"]] <- calculate_variance_explained(object)
    r2 <- rowSums(sapply(object@cache[["variance_explained"]]$r2_per_factor, function(e) rowSums(e, na.rm = TRUE)))
    order_factors <- c(names(r2)[order(r2, decreasing = TRUE)])
    object <- subset_factors(object, order_factors)
  }

  # Mask outliers
  if (remove_outliers) {
    object <- .detect_outliers(object)
  }

  ######################
  ## Quality controls ##
  ######################

  object <- quality_control(object, verbose = verbose)

  return(object)
}
