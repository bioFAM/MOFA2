# (Hidden) General function to set names
# Entity is features, samples, or factors
.set_expectations_names <- function(object, entity, values, views="all", groups="all") {

  stopifnot(entity %in% c("features", "samples", "factors"))

  # Define what entities should be updated for which nodes
  # Notation for axes: 2 is for columns, 1 is for rows, 0 is for vectors
  node_lists_options <- list(features = list(nodes = c("Y", "Tau", "W"), axes = c(1, 1, 1, 1)),
                             samples  = list(nodes = c("Y", "Tau", "Z"), axes = c(2, 2, 1, 1)),
                             factors  = list(nodes = c("Z", "W", "AlphaZ", "AlphaW", "ThetaZ", "ThetaW"), axes = c(2, 2, 0, 0, 0, 0)))

  if (paste0(views, collapse = "") == "all") { 
    views <- names(object@dimensions$D)
  } else {
    stopifnot(all(views %in% names(object@dimensions$D)))
  }

  if (paste0(groups, collapse = "") == "all") {
    groups <- names(object@dimensions$N)
  } else {
    stopifnot(all(groups %in% names(object@dimensions$N)))
  }

  # Iterate over node list depending on the entity
  nodes <- node_lists_options[[entity]]$nodes
  axes  <- node_lists_options[[entity]]$axes
  for (i in 1:length(nodes)) {
    node <- nodes[i]
    axis <- axes[i]
    # Update the nodes for which expectations do exist
    if (node %in% names(object@expectations)) {
      # Update nodes with one level of nestedness (views or groups)
      if (any(node %in% object@model_options$nodes$multiview_node, 
              node %in% object@model_options$nodes$multigroup_nodes)) {
        sub_dim <- length(object@expectations[[node]])
        for (ind in 1:sub_dim) {
          if (all(length(object@expectations[[node]][[ind]]) == 1, 
                  names(object@expectations[[node]][[ind]]) == c("E"))) {
            object@expectations[[node]][[ind]] <- object@expectations[[node]][[ind]]$E
          }
          # No nestedness in values if factors
          vals <- if (entity == "factors") values else values[[ind]]
          dim  <- length(vals)
          # Set names for rows
          if (axis == 1) {
            stopifnot(nrow(object@expectations[[node]][[ind]]) == dim)
            rownames(object@expectations[[node]][[ind]]) <- vals
          # ... or set names for columns
          } else if (axis == 2) {
            stopifnot(ncol(object@expectations[[node]][[ind]]) == dim)
            colnames(object@expectations[[node]][[ind]]) <- vals
          # ... or set vector names
          } else {
            stopifnot(length(object@expectations[[node]][[ind]]) == dim)
            names(object@expectations[[node]][[ind]]) <- vals
          }
        }
      # Update nodes with two levels of nestedness (e.g. Y or Tau)
      } else if (node %in% object@model_options$nodes$twodim_nodes) {
        sub_dim <- length(object@expectations[[node]])
        for (ind in 1:sub_dim) {
          sub_dim2 <- length(object@expectations[[node]][[ind]])
          for (ind2 in 1:sub_dim2) {
            if (all(length(object@expectations[[node]][[ind]][[ind2]]) == 1, 
                    names(object@expectations[[node]][[ind]][[ind2]]) == c("E"))) {
              object@expectations[[node]][[ind]][[ind2]] <- object@expectations[[node]][[ind]][[ind2]]$E
            }
            # Infer which index to use to iterate over a provided list of values
            deduced_ind <- if (entity == "features") ind else ind2  # since ind corresponds to views (groups of features)
            dim <- length(values[[deduced_ind]])
            # Set names for rows
            if (axis == 1) {
              stopifnot(nrow(object@expectations[[node]][[ind]][[ind2]]) == dim)
              rownames(object@expectations[[node]][[ind]][[ind2]]) <- values[[deduced_ind]]
            # ... or set names for columns
            } else if (axis == 2) {
              stopifnot(ncol(object@expectations[[node]][[ind]][[ind2]]) == dim)
              colnames(object@expectations[[node]][[ind]][[ind2]]) <- values[[deduced_ind]]
            # ... or set vector names
            } else {
              stopifnot(length(object@expectations[[node]][[ind]]) == dim)
              names(object@expectations[[node]][[ind]]) <- vals
            }
          }
        }
      } else {
        print(paste0("DEV :: NOTE: There are no expectations for the node ", node))
      }
    }
  }

  object
}

.set_data_names <- function(object, entity, values, views = "all", groups = "all") {
  
  stopifnot(entity %in% c("features", "samples"))

  axes_options <- list(features = 1, samples = 2)        # features are stored in rows

  if (paste0(views, collapse = "") == "all") { 
    views <- names(object@dimensions$D)
  } else {
    stopifnot(all(views %in% names(object@dimensions$D)))
  } 

  if (paste0(groups, collapse = "") == "all") {
    groups <- names(object@dimensions$N)
  } else {
    stopifnot(all(groups %in% names(object@dimensions$N)))
  }
  for (m in views) {
    for (p in groups) {
      deduced_ind <- if (entity == "features") m else p  # since ind corresponds to views (groups of features)
      if (axes_options[[entity]] == 1) {
        rownames(object@training_data[[m]][[p]]) <- values[[deduced_ind]]
      } else {
        colnames(object@training_data[[m]][[p]]) <- values[[deduced_ind]]
      }
    }
  }

  object
}


# (Hidden) General function to set names
# Entity is features, samples, or factors
.set_parameters_names <- function(object, entity, values, views="all", groups="all") {

  tryCatch({
    
    stopifnot(entity %in% c("features", "samples", "factors"))

    # Define what entities should be updated for which nodes
    # Notation for axes: 2 is for columns, 1 is for rows, 0 is for vectors
    node_lists_options <- list(features = list(nodes = c("W"), axes = c(1)),
                               samples  = list(nodes = c("Z"), axes = c(1)),
                               factors  = list(nodes = c("Z", "W", "AlphaZ", "AlphaW", "ThetaZ", "ThetaW"), axes = c(2, 2, 0, 0, 0, 0)))

    views  <- .check_and_get_views(object, views)
    groups <- .check_and_get_views(object, groups)

    # Define the dimensionality of noise
    # If not found in model_options, the default setting is 
    # having a noise parameter per feature
    if ('noise_on' %in% names(object@model_options)) {
      if (object@model_options$noise_on == "samples") {
        node_lists_options$samples$nodes <- c(node_lists_options$samples$nodes, "Tau")
        node_lists_options$samples$axes  <- c(node_lists_options$samples$axes, 0)
      } else {
        node_lists_options$features$nodes <- c(node_lists_options$features$nodes, "Tau")
        node_lists_options$features$axes  <- c(node_lists_options$features$axes, 0)
      }
    } else {
      node_lists_options$features$nodes <- c(node_lists_options$features$nodes, "Tau")
      node_lists_options$features$axes  <- c(node_lists_options$features$axes, 0)
    }

    # Iterate over node list depending on the entity
    nodes <- node_lists_options[[entity]]$nodes
    axes  <- node_lists_options[[entity]]$axes
    for (i in 1:length(nodes)) {
      node <- nodes[i]
      axis <- axes[i]
      # Update the nodes for which parameters do exist
      if (node %in% names(object@parameters)) {
        # Update nodes with one level of nestedness (views or groups)
        if (any(node %in% object@model_options$nodes$multiview_node, 
                node %in% object@model_options$nodes$multigroup_nodes)) {
          sub_dim <- length(object@parameters[[node]])
          for (ind in 1:sub_dim) {
            # No nestedness in values if factors
            vals <- if (entity == "factors") values else values[[ind]]
            dim  <- length(vals)
            for (par in names(object@parameters[[node]][[ind]])) {
              # Set names for rows
              if (axis == 1) {
                stopifnot(nrow(object@parameters[[node]][[ind]][[par]]) == dim)
                rownames(object@parameters[[node]][[ind]][[par]]) <- vals
              # ... or set names for columns
              } else if (axis == 2) {
                stopifnot(ncol(object@parameters[[node]][[ind]][[par]]) == dim)
                colnames(object@parameters[[node]][[ind]][[par]]) <- vals
              # ... or set vector names
              } else {
                stopifnot(length(object@parameters[[node]][[ind]][[par]]) == dim)
                names(object@parameters[[node]][[ind]][[par]]) <- vals
              }
            }
          }
        # Update nodes with two levels of nestedness (e.g. Y or Tau)
        } else if (node %in% object@model_options$nodes$twodim_nodes) {
          sub_dim <- length(object@parameters[[node]])
          for (ind in 1:sub_dim) {
            sub_dim2 <- length(object@parameters[[node]][[ind]])
            for (ind2 in 1:sub_dim2) {
              # Infer which index to use to iterate over a provided list of values
              deduced_ind <- if (entity == "features") ind else ind2  # since ind corresponds to views (groups of features)
              dim <- length(values[[deduced_ind]])
              for (par in names(object@parameters[[node]][[ind]][[ind2]])) {
                # Set names for rows
                if (axis == 1) {
                  stopifnot(nrow(object@parameters[[node]][[ind]][[ind2]][[par]]) == dim)
                  rownames(object@parameters[[node]][[ind]][[ind2]][[par]]) <- values[[deduced_ind]]
                # ... or set names for columns
                } else if (axis == 2) {
                  stopifnot(ncol(object@parameters[[node]][[ind]][[ind2]][[par]]) == dim)
                  colnames(object@parameters[[node]][[ind]][[ind2]][[par]]) <- values[[deduced_ind]]
                # ... or set vector names
                } else {
                  stopifnot(length(object@parameters[[node]][[ind]][[par]]) == dim)
                  names(object@parameters[[node]][[ind]][[par]]) <- vals
                }
              }
            }
          }
        } else {
          print(paste0("DEV :: NOTE: There are no parameters for the node ", node))
        }
      }
    }
  
  }, error = function(e) { print("Oh, something is wrong with setting parameters names.") })

  object

  
}


###################################
## Set and retrieve factor names ##
###################################

#' @rdname factors_names
#' @param object a \code{\link{BioFAM}} object.
#' @aliases factors_names,BioFAModel-method
#' @return character vector with the features names
#' @export
setMethod("factors_names", signature(object="BioFAModel"), 
  function(object) {
    colnames(object@expectations$Z[[1]]) 
  }
)

#' @rdname factors_names
#' @param value a character vector of factor names
#' @import methods
#' @export
setReplaceMethod("factors_names", signature(object="BioFAModel", value="vector"), 
  function(object, value) {
    if (!methods::.hasSlot(object, "expectations") | length(object@expectations) == 0)
      stop("Before assigning factor names you have to assign expectations")
    if (methods::.hasSlot(object, "dimensions") & length(object@dimensions) != 0)
      if (length(value) != object@dimensions["K"])
        stop("Length of factor names does not match the dimensionality of the latent variable matrix")
    # DEPRECATED: If it's an old model, it should be object@expectations$Z$E
    # if(length(value) != ncol(object@expectations$Z[[1]])) 
    #   stop("Factor names do not match the number of columns in the latent variable matrix")
 
    object <- .set_expectations_names(object, entity = 'factors', value)
    object <- .set_parameters_names(object, entity = 'factors', value)
 
    object
  })



###################################
## Set and retrieve sample names ##
###################################

#' @rdname samples_names
#' @param object a \code{\link{BioFAModel}} object.
#' @aliases samples_names,BioFAModel-method
#' @return list of character vectors with the sample names for each group
#' @export
setMethod("samples_names", signature(object="BioFAModel"), 
  function(object) { 
    object@data_options$samples_names
  })

#' @rdname samples_names
#' @param value list of character vectors with the sample names for every group
#' @import methods
#' @export
setReplaceMethod("samples_names", signature(object="BioFAModel", value="list"), 
  function(object, value) {
   if (!methods::.hasSlot(object, "training_data") | length(object@training_data) == 0 | length(object@training_data[[1]]) == 0)
     stop("Before assigning sample names you have to assign the training data")
   if (!methods::.hasSlot(object, "expectations") | length(object@expectations) == 0)
     stop("Before assigning sample names you have to assign the expectations")
   if (methods::.hasSlot(object, "dimensions") & length(object@dimensions) != 0)
    if (!all(sapply(value, length) == object@dimensions[["N"]]))
        stop("Length of sample names does not match the dimensionality of the model")
    if (!all(sapply(value, length) == sapply(object@training_data[[1]], ncol)))
      stop("sample names do not match the dimensionality of the data (columns)")

    object@data_options$samples_names <- value
    object <- .set_expectations_names(object, entity = 'samples', value)
    object <- .set_parameters_names(object, entity = 'samples', value)
    object <- .set_data_names(object, entity = 'samples', value)

    object
  })

####################################
## Set and retrieve feature names ##
####################################

#' @rdname features_names
#' @param object a \code{\link{BioFAModel}} object.
#' @aliases features_names,BioFAModel-method
#' @return list of character vectors with the feature names for each view
#' @export
setMethod("features_names", signature(object="BioFAModel"), 
  function(object) { 
    object@data_options$features_names
  })

#' @rdname features_names
#' @param value list of character vectors with the feature names for every view
#' @import methods
#' @export
setReplaceMethod("features_names", signature(object="BioFAModel", value="list"),
  function(object, value) {
    if (!methods::.hasSlot(object, "training_data")  | length(object@training_data) == 0)
      stop("Before assigning feature names you have to assign the training data")
    if (!methods::.hasSlot(object, "expectations") | length(object@expectations) == 0)
      stop("Before assigning feature names you have to assign the expectations")
    if (methods::.hasSlot(object, "dimensions")  | length(object@dimensions) == 0)
      if (!all(sapply(value, length) == object@dimensions[["D"]]))
        stop("Length of feature names does not match the dimensionality of the model")
    if (!all(sapply(value, length) == sapply(object@training_data, function(e) nrow(e[[1]]))))
      stop("Feature names do not match the dimensionality of the data (rows)")
    
    object@data_options$features_names <- value
    object <- .set_expectations_names(object, entity = 'features', value)
    object <- .set_parameters_names(object, entity = 'features', value)
    object <- .set_data_names(object, entity = 'features', value)
    
    return(object)
  })

#################################
## Set and retrieve view names ##
#################################

#' @rdname views_names
#' @param object a \code{\link{BioFAModel}} object.
#' @return character vector with the names for each view
#' @rdname views_names
#' @export
setMethod("views_names", signature(object="BioFAModel"), 
  function(object) {
    object@data_options$views_names
  })


#' @rdname views_names
#' @param value character vector with the names for each view
#' @import methods
#' @export
setMethod("views_names<-", signature(object="BioFAModel", value="character"), 
  function(object, value) {
    # if (!methods::.hasSlot(object, "training_data") | length(object@training_data) == 0)
    #   stop("Before assigning view names you have to assign the training data")
    if (methods::.hasSlot(object, "dimensions") & length(object@dimensions) != 0)
      if (length(value) != object@dimensions["M"])
      stop("Length of view names does not match the dimensionality of the model")
    # if (length(value) != length(object@training_data))
    #   stop("View names do not match the number of views in the training data")

    # Set view names in data options
    object@data_options$views_names <- value
  
    # Set view names in expectations
    for (node in names(object@expectations)) {
      if (node %in% object@model_options$nodes$multiview_nodes |
          node %in% object@model_options$nodes$twodim_nodes) {
        if (class(object@expectations[[node]]) == "list" & length(object@expectations[[node]]) == object@dimensions["M"]) {
          names(object@expectations[[node]]) <- value 
        }
      }
    }
    
    # Set view names in the raining data
    if (length(object@training_data)>0) {
      names(object@training_data) <- value
    }
    
    # Set view names in the dimensionalities
    names(object@dimensions$D) <- value
    
    return(object)
  })


#################################
## Set and retrieve group names ##
#################################

#' @rdname groups_names
#' @param object a \code{\link{BioFAModel}} object.
#' @return character vector with the names for each view
#' @rdname groups_names
#' @export
setMethod("groups_names", signature(object="BioFAModel"), 
  function(object) {
    object@data_options$samples_groups
  })


#' @rdname groups_names
#' @param value character vector with the names for each view
#' @import methods
#' @export
setMethod("groups_names<-", signature(object="BioFAModel", value="character"), 
  function(object, value) {
    # if (!methods::.hasSlot(object, "training_data") | length(object@training_data) == 0)
    #   stop("Before assigning group names you have to assign the training data")
    if (methods::.hasSlot(object,"dimensions") & length(object@dimensions) != 0)
      if(length(value) != object@dimensions["P"])
        stop("Length of group names does not match the dimensionality of the model")
    # if (length(value) != length(object@training_data[[1]]))
    #   stop("Group names do not match the number of groups in the training data")
    
    # Set sample group names in data options
    object@data_options$samples_groups <- value
    
    # Set sample group names in expectations
    for (node in object@model_options$nodes$multigroup_nodes) {
      if (node %in% names(object@expectations)) {
        if (class(object@expectations[[node]]) == "list" & 
            length(object@expectations[[node]]) == object@dimensions["P"]) {
          names(object@expectations[[node]]) <- value 
        }
      }
    }
    for (node in object@model_options$nodes$twodim_nodes) {
      if (node %in% names(object@expectations)) {
        for (m in 1:length(object@expectations$Y)) {
            if (class(object@expectations[[node]][[m]]) == "list" & 
                length(object@expectations[[node]][[m]]) == object@dimensions["P"]) {
              names(object@expectations[[node]][[m]]) <- value 
            }
        }
      }
    }
    
    # Set sample group names in training data
    if (length(object@training_data)>0) {
      for (m in names(object@training_data))
        names(object@training_data[[m]]) <- value
    }
    
    # Set sample group names in dimensionalities
    names(object@dimensions$N) <- value
    
    return(object)
  })

#################################
## Set and retrieve input data ##
#################################

#' @title Set and retrieve input data
#' @name input_data
#' @param object a \code{\link{BioFAModel}} object.
#' @rdname input_data
#' @export
setMethod("input_data", signature(object="BioFAModel"), 
  function(object) {
    object@input_data
  })

#' @title Set and retrieve input data
#' @docType methods
#' @name input_data
#' @param object a \code{\link{BioFAModel}} object.
#' @rdname input_data
#' @aliases input_data<-
#' @export
# setMethod(".input_data<-", signature(object="BioFAModel", value="MultiAssayExperiment"),
setMethod(".input_data<-", signature(object="BioFAModel", value="data.frame"),
  function(object,value) {
    object@input_data <- value
    object
  })

####################################
## Set and retrieve training data ##
####################################

#' @rdname training_data
#' @param object a \code{\link{BioFAModel}} object.
#' @return list of numeric matrices that contain the training data
#' @rdname training_data
#' @export
setMethod("training_data", signature(object="BioFAModel"),
  function(object) {
    object@training_data 
  })

#' @import methods
setMethod(".training_data<-", signature(object="BioFAModel", value="list"),
  function(object,value) {
    # N <- unique(sapply(value,nrow))
    # if (length(N) > 1) 
    #   stop("Views do not have the same number of samples (rows)")
    # if (methods::.hasSlot(object,"dimensions")) {
    #   if (object@dimensions["M"] != length(value))
    #     if (object@dimensions["N"] != N)
    #       stop("Number of samples in the data do not match the specified dimensionality of the model")
    #   if (all(object@dimensions["D"] != sapply(value,ncol)))
    #     stop("Number of features in the data do not match the specified dimensionality of the model")
    # }
    object@training_data <- value
    object
  })

####################################
## Set and retrieve imputed data ##
####################################

#' @rdname imputed_data
#' @param object a \code{\link{BioFAModel}} object.
#' @return list of numeric matrices that contain the training data
#' @rdname imputed_data
#' @export
setMethod("imputed_data", signature(object="BioFAModel"), function(object) { object@imputed_data } )

#' @import methods
setMethod(".imputed_data<-", signature(object="BioFAModel", value="list"),
  function(object,value) {
    # to do sanity checks
    object@imputed_data <- value
    object
  })

#######################################
## Set and retrieve training options ##
#######################################

#' @rdname training_options
#' @param object a \code{\link{BioFAModel}} object.
#' @rdname training_options
#' @return list of training options
#' @export
setMethod("training_options", "BioFAModel", function(object) { object@training_options } )
setMethod(".training_options<-", signature(object="BioFAModel", value="list"),
   function(object,value) {
     object@training_options <- value
     object
   })

#######################################
## Set and retrieve model options ##
#######################################

#' @rdname model_options
#' @param object a \code{\link{BioFAModel}} object.
#' @rdname model_options
#' @return list of model options
#' @export
setMethod("model_options", "BioFAModel", function(object) { object@model_options } )
setMethod(".model_options<-", signature(object="BioFAModel", value="list"),
    function(object,value) {
      object@model_options <- value
      object
    })

##########################################
## Set and retrieve training statistics ##
##########################################

#' @rdname training_stats
#' @param object a \code{\link{BioFAModel}} object.
#' @return list of training statistics
#' @export
setMethod("training_stats", "BioFAModel", function(object) { object@training_stats } )
setMethod(".training_stats<-", signature(object="BioFAModel", value="list"),
  function(object,value) {
    object@training_stats <- value
    object
  })

###################################
## Set and retrieve expectations ##
###################################

#' @rdname expectations
#' @param object a \code{\link{BioFAModel}} object.
#' @rdname expectations
#' @return list of expectations
#' @export
setMethod("expectations", "BioFAModel", function(object) { object@expectations } )
setMethod(".expectations<-", signature(object="BioFAModel", value="list"),
  function(object,value) {
    object@expectations <- value
    object
  })
