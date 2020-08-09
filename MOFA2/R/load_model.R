
############################################
## Functions to load a trained MOFA model ##
############################################

#' @title Load a trained MOFA
#' @name load_model
#' @description Method to load a trained MOFA \cr
#' The training of mofa is done using a Python framework, and the model output is saved as an .hdf5 file, which has to be loaded in the R package.
#' @param file an hdf5 file saved by the mofa Python framework
# #' @param sort_factors logical indicating whether factors should be sorted by variance explained (default is TRUE)
#' @param on_disk logical indicating whether to work from memory (FALSE) or disk (TRUE). \cr
#' This should be set to TRUE when the training data is so big that cannot fit into memory. \cr
#' On-disk operations are performed using the \code{\link{HDF5Array}} and \code{\link{DelayedArray}} framework.
#' @param load_data logical indicating whether to load the training data (default is TRUE, it can be memory expensive)
#' @param load_imputed_data logical indicating whether to load the imputed data (default is FALSE)
#' @param remove_outliers logical indicating whether to mask outlier values.
#' @param remove_inactive_factors logical indicating whether to remove inactive factors from the model.
# #' @param remove_intercept_factors logical indicating whether to remove intercept factors for non-Gaussian views.
#' @param verbose logical indicating whether to print verbose output (default is FALSE)
#' @return a \code{\link{MOFA}} model
#' @importFrom rhdf5 h5read h5ls
#' @importFrom HDF5Array HDF5ArraySeed
# #' @importFrom DelayedArray DelayedArray
#' @export
load_model <- function(file, sort_factors = TRUE, on_disk = FALSE, load_data = TRUE, load_imputed_data = FALSE, 
                       remove_outliers = FALSE, remove_inactive_factors = TRUE, verbose = FALSE) {

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

  # Load names
  if ("views" %in% h5ls.out$name) {
    view_names <- as.character( h5read(file, "views")[[1]] )
    group_names <- as.character( h5read(file, "groups")[[1]] )
    feature_names <- h5read(file, "features")[view_names]
    sample_names  <- h5read(file, "samples")[group_names] 
  } else {  # for old models
    feature_names <- h5read(file, "features")
    sample_names  <- h5read(file, "samples")
    view_names <- names(feature_names)
    group_names <- names(sample_names)
    h5ls.out <- h5ls.out[grep("variance_explained", h5ls.out$name, invert = TRUE),]
  }

  # Load training data (as nested list of matrices)
  data <- list(); intercepts <- list()
  if (isTRUE(load_data) & ("data"%in%h5ls.out$name) ) {
    
    object@data_options[["loaded"]] <- TRUE
    if (isTRUE(verbose)) message("Loading data...")
    
    for (m in view_names) {
      data[[m]] <- list()
      intercepts[[m]] <- list()
      for (g in group_names) {
        if (on_disk) {
          # as DelayedArrays
          data[[m]][[g]] <- DelayedArray::DelayedArray( HDF5ArraySeed(file, name = sprintf("data/%s/%s", m, g) ) )
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
    
    object@data_options[["loaded"]] <- FALSE
    
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


  # Load metadata if any
  if ("samples_metadata" %in% h5ls.out$name) {
    object@samples_metadata <- do.call(rbind, lapply(group_names, function(g) as.data.frame(h5read(file, sprintf("samples_metadata/%s", g)))))
  }
  if ("features_metadata" %in% h5ls.out$name) {
    object@features_metadata <- do.call(rbind, lapply(view_names, function(m) as.data.frame(h5read(file, sprintf("features_metadata/%s", m)))))
  }
  
  #######################
  ## Load imputed data ##
  #######################
  
  imputed_data <- list()
  if (isTRUE(load_imputed_data)) {
    
    if (isTRUE(verbose)) message("Loading imputed data...")
    
    for (m in view_names) {
      imputed_data[[m]] <- list()
      for (g in group_names) {
        imputed_data[[m]][[g]] <- list()
        if (on_disk) {
          # as DelayedArrays
          # imputed_data[[m]][[g]] <- DelayedArray::DelayedArray( HDF5ArraySeed(file, name = sprintf("imputed_data/%s/%s", m, g) ) )
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
    expectations[["AlphaW"]] <- h5read(file, "expectations/AlphaW")[view_names]
  if ("AlphaZ" %in% node_names)
    expectations[["AlphaZ"]] <- h5read(file, "expectations/AlphaZ")[group_names]
  if ("Z" %in% node_names)
    expectations[["Z"]] <- h5read(file, "expectations/Z")[group_names]
  if ("W" %in% node_names)
    expectations[["W"]] <- h5read(file, "expectations/W")[view_names]
  if ("ThetaW" %in% node_names)
    expectations[["ThetaW"]] <- h5read(file, "expectations/ThetaW")[view_names]
  if ("ThetaZ" %in% node_names)
    expectations[["ThetaZ"]] <- h5read(file, "expectations/ThetaZ")[group_names]
  # if ("Tau" %in% node_names)
  #   expectations[["Tau"]] <- h5read(file, "expectations/Tau")
  
  # Load expectations of Y
  # tmp <- as.list(h5read(file, 'model_options', read.attributes = TRUE))[["likelihoods"]]
  # names(tmp) <- view_names
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
  for (i in names(object@model_options)) {
    if (object@model_options[i] == "False" | object@model_options[i] == "True") {
      object@model_options[i] <- as.logical(object@model_options[i])
    } else {
      object@model_options[i] <- object@model_options[i]
    }
  }

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
    object@training_stats <- h5read(file, 'training_stats', read.attributes = TRUE)
  }, error = function(x) { print("Training stats not found, not loading it...") })

  #######################################
  ## Load variance explained estimates ##
  #######################################
  
  if ("variance_explained" %in% h5ls.out$name) {
    r2_list <- list(
      r2_total = h5read(file, "variance_explained/r2_total"),
      r2_per_factor = h5read(file, "variance_explained/r2_per_factor")
    )
    object@cache[["variance_explained"]] <- r2_list
  }
  
  # Hack to fix the problems where variance explained values range from 0 to 1 (%)
  if (max(sapply(object@cache$variance_explained$r2_total,max))<1) {
    for (m in 1:length(view_names)) {
      for (g in 1:length(group_names)) {
        object@cache$variance_explained$r2_total[[g]][[m]] <- 100 * object@cache$variance_explained$r2_total[[g]][[m]]
        object@cache$variance_explained$r2_per_factor[[g]][,m] <- 100 * object@cache$variance_explained$r2_per_factor[[g]][,m]
      }
    }
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
  if (isTRUE(verbose)) message("Assigning names to the different dimensions...")

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
  features_names(object) <- feature_names
  
  # Create default samples names if they are null
  if (is.null(sample_names)) {
    print("Samples names not found, generating default: sample1, ..., sampleN")
    sample_names <- lapply(object@dimensions[["N"]], function(n) paste0("sample", as.character(seq_len(n))))
  }
  samples_names(object) <- sample_names
  
  # Set views names
  if (is.null(names(object@data))) {
    print("Views names not found, generating default: view1, ..., viewM")
    view_names <- paste0("view", as.character(seq_len(object@dimensions[["M"]])))
  }
  views_names(object) <- view_names
  
  # Set groups names
  if (is.null(names(object@data[[1]]))) {
    print("Groups names not found, generating default: group1, ..., groupG")
    group_names <- paste0("group", as.character(seq_len(object@dimensions[["G"]])))
  }
  groups_names(object) <- group_names
  
  # Set factors names
  factors_names(object)  <- paste0("Factor", as.character(seq_len(object@dimensions[["K"]])))
  
  ###################
  ## Parse factors ##
  ###################
  
  # Calculate variance explained estimates per factor
  if (is.null(object@cache[["variance_explained"]])) {
    object@cache[["variance_explained"]] <- calculate_variance_explained(object)
  } 
  
  # Remove inactive factors
  if (isTRUE(remove_inactive_factors)) {
    r2 <- rowSums(do.call('cbind', lapply(object@cache[["variance_explained"]]$r2_per_factor, rowSums, na.rm=TRUE)))
    var.threshold <- 0.0001
    if (all(r2 < var.threshold)) {
      warning(sprintf("All %s factors were found to explain little or no variance so remove_inactive_factors option has been disabled.", length(r2)))
    } else if (any(r2 < var.threshold)) {
      object <- subset_factors(object, which(r2>=var.threshold))
      message(sprintf("%s factors were found to explain little or no variance and they were removed for downstream analysis. You can disable this option by setting load_model(..., remove_inactive_factors = F)", sum(r2 < var.threshold)))
    }
  }
  
  # Remove inactive factors
  # if (isTRUE(remove_intercept_factors) & any(object@model_options$likelihoods!="gaussian")) {
  #   non_gaussian_views <- names(which(object@model_options$likelihoods!="gaussian"))
  #   for (m in non_gaussian_views) {
  #     W <- get_weights(object, views=m)[[1]]
  #     intercept_factor <- which.max(abs(colMeans(W)))
  #     object <- subset_factors(object, factors_names(object)[-intercept_factor])
  #     message(sprintf("%s is determined to be an intercept Factor for the %s non-gaussian view and it has been removed from the model. To disable this behaviour set remove_intercept_factors=FALSE",names(intercept_factor),m))
  #   }
  # }
  
  # [Done in mofapy2] Sort factors by total variance explained
  if (isTRUE(sort_factors) && object@dimensions$K>1) {

    # Sanity checks
    if (isTRUE(verbose)) message("Re-ordering factors by their variance explained...")

    # Calculate variance explained per factor across all views
    r2 <- rowSums(sapply(object@cache[["variance_explained"]]$r2_per_factor, function(e) rowSums(e, na.rm = TRUE)))
    order_factors <- c(names(r2)[order(r2, decreasing = TRUE)])

    # re-order factors
    object <- subset_factors(object, order_factors)
  }

  # Mask outliers
  if (isTRUE(remove_outliers)) {
    if (isTRUE(verbose)) message("Removing outliers...")
    object <- .detect_outliers(object)
  }
  
  # Mask intercepts for non-Gaussian data
  if (any(object@model_options$likelihoods!="gaussian")) {
    for (m in names(which(object@model_options$likelihoods!="gaussian"))) {
      for (g in names(object@intercepts[[m]])) {
        object@intercepts[[m]][[g]] <- NA
      }
    }
  }

  ######################
  ## Quality controls ##
  ######################

  if (isTRUE(verbose)) message("Doing quality control...")
  object <- quality_control(object, verbose = verbose)
  
  return(object)
}
