

####################################
## Set and retrieve factors names ##
####################################

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
                   
                   object <- .set_expectations_names(object, entity = 'factors', value)
                   
                   object
                 })



####################################
## Set and retrieve samples names ##
####################################

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
                   object <- .set_data_names(object, entity = 'samples', value)
                   
                   object
                 })


#############################
## Retrieve samples groups ##
#############################

#' @rdname samples_groups
#' @param object a \code{\link{BioFAModel}} object.
#' @aliases samples_groups, BioFAModel-method
#' @return data.frame with the sample names and a group for each sample
#' @export
setMethod("samples_groups", signature(object="BioFAModel"), 
          function(object, format = "default") {
            samples <- object@data_options$samples_names
            tmp <- data.frame(sample = unlist(samples),
                              group = rep(names(samples), sapply(samples, length)),
                              row.names = c())
            if (format == "pheatmap") {
              rownames(tmp) <- unlist(samples)
              tmp <- tmp[,"group", drop = FALSE]
              colnames(tmp) <- "ID"
            }
            tmp
          })


#####################################
## Set and retrieve features names ##
#####################################

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
                   object <- .set_data_names(object, entity = 'features', value)
                   
                   return(object)
                 })

##################################
## Set and retrieve views names ##
##################################

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
            
            # Define types of nodes
            nodes_types <- list(
              multiview_nodes  = c("W", "AlphaW", "ThetaW"),
              multigroup_nodes = c("Z", "AlphaZ", "ThetaZ"),
              twodim_nodes     = c("Y", "Tau")
            )
            
            # Set view names in data options
            object@data_options$views_names <- value
            if (!is.null(object@data_options$features_names)) {
              names(object@data_options$features_names) <- value 
            }
            
            # Set view names in expectations
            for (node in names(object@expectations)) {
              if (node %in% nodes_types$multiview_nodes | node %in% nodes_types$twodim_nodes) {
                if (class(object@expectations[[node]]) == "list" & length(object@expectations[[node]]) == object@dimensions["M"]) {
                  names(object@expectations[[node]]) <- value 
                }
              }
            }
            
            # Set view names in the training data
            if (length(object@training_data)>0)
              names(object@training_data) <- value
            
            # Set view names in the input data
            if (length(object@input_data)>0)
              names(object@input_data) <- value
            
            # Set view names in the dimensionalities
            names(object@dimensions$D) <- value
            
            return(object)
          })


###################################
## Set and retrieve groups names ##
###################################

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
            
            # Define types of nodes
            nodes_types <- .get_nodes_types()

            # Set sample group names in data options
            object@data_options$samples_groups <- value
            if (!is.null(object@data_options$samples_names)) {
              names(object@data_options$samples_names) <- value
            }
              
            
            # Set sample group names in expectations
            for (node in nodes_types$multigroup_nodes) {
              if (node %in% names(object@expectations)) {
                if (class(object@expectations[[node]])=="list" & length(object@expectations[[node]])==object@dimensions["P"]) {
                  names(object@expectations[[node]]) <- value 
                }
              }
            }
            
            for (node in nodes_types$twodim_nodes) {
              if (node %in% names(object@expectations)) {
                for (m in 1:length(object@expectations[[node]])) {
                  if (class(object@expectations[[node]][[m]])=="list" & length(object@expectations[[node]][[m]])==object@dimensions["P"]) {
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
            stopifnot(length(object@dimensions$N)==length(value))
            names(object@dimensions$N) <- value
            
            return(object)
          })


# (Hidden) General function to set dimension names for the expectations
# Entity is features, samples, or factors
.set_expectations_names <- function(object, entity, values, views="all", groups="all") {
  
  # Define types of nodes
  nodes_types <- .get_nodes_types()
  
  # Define what entities should be updated for which nodes
  #   Notation for axes: 2 is for columns, 1 is for rows, 0 is for vectors
  stopifnot(entity %in% c("features", "samples", "factors"))
  node_lists_options <- list(
    features = list(nodes = c("Y", "Tau", "W"), axes = c(1, 1, 1, 1)),
    samples  = list(nodes = c("Y", "Tau", "Z"), axes = c(2, 2, 1, 1)),
    factors  = list(nodes = c("Z", "W", "AlphaZ", "AlphaW", "ThetaZ", "ThetaW"), axes = c(2, 2, 0, 0, 0, 0))
  )
  
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
      
      # Update nodes with one level of nestedness (e.g. W or Z)
      if (any(node %in% nodes_types$multiview_node, node %in% nodes_types$multigroup_nodes)) {
        sub_dim <- length(object@expectations[[node]])
        for (ind in 1:sub_dim) {
          
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
          } else if (axis == 0) {
            stopifnot(length(object@expectations[[node]][[ind]]) == dim)
            names(object@expectations[[node]][[ind]]) <- vals
          }
        }
        
      # Update nodes with two levels of nestedness (e.g. Y or Tau)
      } else if (node %in% nodes_types$twodim_nodes) {
        sub_dim <- length(object@expectations[[node]])
        for (ind in 1:sub_dim) {
          sub_dim2 <- length(object@expectations[[node]][[ind]])
          for (ind2 in 1:sub_dim2) {
            
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


# (Hidden) General function to set dimensions names for the data
.set_data_names <- function(object, entity, values) {
  
  stopifnot(entity %in% c("features", "samples"))
  
  axes_options <- list(features = 1, samples = 2)
  
  for (m in 1:length(object@training_data)) {
    for (p in 1:length(object@training_data[[m]])) {
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

