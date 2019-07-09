

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
                   
                   # Modify expectations
                   object <- .set_expectations_names(object, entity = 'factors', value)
                   
                   # Modify cache
                   if ((methods::.hasSlot(object, "cache")) & ("variance_explained" %in% names(object@cache))) {
                     for (i in 1:length(object@cache$variance_explained$r2_per_factor)) {
                       rownames(object@cache$variance_explained$r2_per_factor[[i]]) <- value
                     }
                   }
                     
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
            
            # When the model is not trained, the samples slot is not initialized yet
            if (!("samples_metadata" %in% slotNames(object)) || (length(object@samples_metadata) == 0)) {
              return(list())
            }
            
            # The default case when samples are initialized (trained model)
            samples_list <- lapply(object@data_options$groups, function(g) {
              with(object@samples_metadata, object@samples_metadata[group_name == g, "sample_name"])
            })
            
            names(samples_list) <- object@data_options$groups
            return(samples_list)
          })

#' @rdname samples_names
#' @param value list of character vectors with the sample names for every group
#' @import methods
#' @export
setReplaceMethod("samples_names", signature(object="BioFAModel", value="list"), 
                 function(object, value) {
                   if (!methods::.hasSlot(object, "data") | length(object@data) == 0 | length(object@data[[1]]) == 0)
                     stop("Before assigning sample names you have to assign the training data")
                   if (!methods::.hasSlot(object, "expectations") | length(object@expectations) == 0)
                     stop("Before assigning sample names you have to assign the expectations")
                   if (methods::.hasSlot(object, "dimensions") & length(object@dimensions) != 0)
                     if (!all(sapply(value, length) == object@dimensions[["N"]]))
                       stop("Length of sample names does not match the dimensionality of the model")
                   if (!all(sapply(value, length) == sapply(object@data[[1]], ncol)))
                     stop("sample names do not match the dimensionality of the data (columns)")
                   
                   value_groups <- rep(names(value), lengths(value))

                   # Modify sample names in the sample metadata
                   object@samples_metadata$sample_name <- unlist(value, use.names = FALSE)
                   object@samples_metadata$group_name  <- as.factor( value_groups )
                   if (class(object@samples_metadata) == "list") {
                    object@samples_metadata <- data.frame(object@samples_metadata, stringsAsFactors=FALSE)
                   }
                   
                   # Add samples names to the expectations matrices
                   object <- .set_expectations_names(object, entity = 'samples', value)
                   
                   # Add samples names to the data matrices
                   object <- .set_data_names(object, entity = 'samples', value)
                   
                   object
                 })


#############################
## Retrieve samples groups ##
#############################

#' @rdname groups
#' @param object a \code{\link{BioFAModel}} object.
#' @aliases groups, BioFAModel-method
#' @return data.frame with the sample names and a group for each sample
#' @export
setMethod("groups", signature(object="BioFAModel"), 
          function(object, format = "default") {
            tmp <- data.frame(sample_name = object@samples_metadata$sample_name,
                              group_name  = object@samples_metadata$group_name,
                              row.names   = c())
            if (format == "pheatmap") {
              rownames(tmp) <- tmp$sample_name
              tmp <- tmp[,"group_name", drop = FALSE]
              colnames(tmp) <- "ID"
            }
            tmp
          })

#######################################
## Set and retriveve sample metadata ##
#######################################

#' @rdname samples_metadata
#' @param object a \code{\link{BioFAModel}} object.
#' @return a data frame with sample metadata
#' @export
setMethod("samples_metadata", signature(object="BioFAModel"), 
          function(object) { 
            object@samples_metadata
          })

#' @rdname samples_metadata
#' @param value data frame with sample information, has to contain columns sample_name and group_name
#' @import methods
#' @export
setReplaceMethod("samples_metadata", signature(object="BioFAModel", value="data.frame"), 
                 function(object, value) {
                   if (!methods::.hasSlot(object, "data") | length(object@data) == 0 | length(object@data[[1]]) == 0)
                     stop("Before assigning samples metadata you have to assign the input data")
                   # if (!methods::.hasSlot(object, "expectations") | length(object@expectations) == 0)
                   #   stop("Before assigning samples metadata you have to assign the expectations")
                   if (methods::.hasSlot(object, "dimensions") & length(object@dimensions) != 0)
                     if (nrow(value) != sum(object@dimensions[["N"]]))
                       stop("Number of rows in samples metadata does not match the dimensionality of the model")
                   if (nrow(value) != sum(sapply(object@data[[1]], ncol)))
                     stop("sample names do not match the dimensionality of the data (columns)")
                   if (!("sample_name" %in% colnames(value)))
                     stop("Metadata has to contain the column sample_name")
                   if (!("group_name" %in% colnames(value)))
                     stop("Metadata has to contain the column group_name")
                   
                   object@samples_metadata <- as.data.frame(value)
                   
                   object
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
            # When the model is not trained, the features slot is not initialized yet
            if (!("features_metadata" %in% slotNames(object)) || (length(object@features_metadata) == 0)) {
              return(list())
            }
            # The default case when features are initialized (trained model)
            features_list <- lapply(object@data_options$views, function(g) {
              with(object@features_metadata, object@features_metadata[view_name == g, "feature_name"])
            })
            names(features_list) <- object@data_options$views
            return(features_list)
          })

#' @rdname features_names
#' @param value list of character vectors with the feature names for every view
#' @import methods
#' @export
setReplaceMethod("features_names", signature(object="BioFAModel", value="list"),
                 function(object, value) {
                   if (!methods::.hasSlot(object, "data")  | length(object@data) == 0)
                     stop("Before assigning feature names you have to assign the training data")
                   if (!methods::.hasSlot(object, "expectations") | length(object@expectations) == 0)
                     stop("Before assigning feature names you have to assign the expectations")
                   if (methods::.hasSlot(object, "dimensions")  | length(object@dimensions) == 0)
                     if (!all(sapply(value, length) == object@dimensions[["D"]]))
                       stop("Length of feature names does not match the dimensionality of the model")
                   if (!all(sapply(value, length) == sapply(object@data, function(e) nrow(e[[1]]))))
                     stop("Feature names do not match the dimensionality of the data (rows)")
                   
                   value_groups <- rep(names(value), lengths(value))

                   object@features_metadata$feature_name <- unlist(value, use.names = FALSE)
                   object@features_metadata$view_name <- value_groups

                   if (class(object@features_metadata) == "list") {
                    object@features_metadata <- data.frame(object@features_metadata, stringsAsFactors = FALSE)
                   }
                   
                   # Add features names to the expectations matrices
                   object <- .set_expectations_names(object, entity = 'features', value)
                   
                   # Add features names to the data matrices
                   # if (!is.null(dim(object@data[[1]][[1]]))) {
                   # }
                   object <- .set_data_names(object, entity = 'features', value)
                   
                   object
                 })

########################################
## Set and retriveve feature metadata ##
########################################

#' @rdname features
#' @param object a \code{\link{BioFAModel}} object.
#' @return a data frame with sample metadata
#' @export
setMethod("features_metadata", signature(object="BioFAModel"), 
          function(object) { 
            object@features_metadata
          })

#' @rdname features_metadata
#' @param value data frame with feature information, has to contain columns feature_name and view_name
#' @import methods
#' @export
setReplaceMethod("features_metadata", signature(object="BioFAModel", value="data.frame"), 
                 function(object, value) {
                   if (!methods::.hasSlot(object, "data") | length(object@data) == 0 | length(object@data[[1]]) == 0)
                     stop("Before assigning features metadata you have to assign the training data")
                   if (!methods::.hasSlot(object, "expectations") | length(object@expectations) == 0)
                     stop("Before assigning features metadata you have to assign the expectations")
                   if (methods::.hasSlot(object, "dimensions") & length(object@dimensions) != 0)
                     if (nrow(value) != sum(object@dimensions[["G"]]))
                       stop("Number of rows in features metadata does not match the dimensionality of the model")
                   if (nrow(value) != sum(sapply(object@data[[1]], nrow)))
                     stop("Features names do not match the dimensionality of the data (rows)")
                   if (!("feature_name" %in% colnames(value)))
                     stop("Metadata has to contain the column feature_name")
                   if (!("view_name" %in% colnames(value)))
                     stop("Metadata has to contain the column view_name")
                   if (colnames(value)[1] != "feature_name")
                     message("Note that feature_name is currently not the first column of the features metadata.")
                   
                   object@features_metadata <- value
                   
                   object
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
            object@data_options$views
          })


#' @rdname views_names
#' @param value character vector with the names for each view
#' @import methods
#' @export
setMethod("views_names<-", signature(object="BioFAModel", value="character"), 
          function(object, value) {
            # if (!methods::.hasSlot(object, "data") | length(object@data) == 0)
            #   stop("Before assigning view names you have to assign the training data")
            if (methods::.hasSlot(object, "dimensions") & length(object@dimensions) != 0)
              if (length(value) != object@dimensions[["M"]])
                stop("Length of view names does not match the dimensionality of the model")
            # if (length(value) != length(object@data))
            #   stop("View names do not match the number of views in the training data")
            
            # Define types of nodes
            nodes_types <- .get_nodes_types()
            
            # Set view names in data options
            old_views <- object@data_options$views
            object@data_options$views <- value
            
            # Set view names in model options
            if (length(object@model_options$likelihoods)>0)
              object@model_options$likelihoods <- value
          
            
            # Set view names in features_metadata 
            if (!is.null(object@features_metadata) && (length(object@features_metadata) != 0)) {
              # object@features_metadata$view_name <- as.character(object@features_metadata$view_name)
              for (i in 1:object@dimensions[["M"]]) {
                old_name <- old_views[i]
                new_name <- value[i]
                object@features_metadata[object@features_metadata$view_name == old_name, "view_name"] <- new_name
              }
            }
            
            # Set view names in cache
            if (!is.null(object@cache$variance_explained)) {
              for (i in names(object@cache$variance_explained$r2_total)) {
                names(object@cache$variance_explained$r2_total[[i]]) <- value
                colnames(object@cache$variance_explained$r2_per_factor[[i]]) <- value
              }
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
            if (length(object@data)>0)
              names(object@data) <- value
            
            # Set view names in the dimensionalities
            names(object@dimensions$D) <- value
            
            return(object)
          })


###################################
## Set and retrieve groups names ##
###################################

#' @rdname groups_names
#' @param object a \code{\link{BioFAModel}} object.
#' @return character vector with the names for each sample group
#' @rdname groups_names
#' @export
setMethod("groups_names", signature(object="BioFAModel"), 
          function(object) {
            object@data_options$groups
          })


#' @rdname groups_names
#' @param value character vector with the names for each sample group
#' @import methods
#' @export
setMethod("groups_names<-", signature(object="BioFAModel", value="character"), 
          function(object, value) {
            # if (!methods::.hasSlot(object, "data") | length(object@data) == 0)
            #   stop("Before assigning group names you have to assign the training data")
            if (methods::.hasSlot(object, "dimensions") & length(object@dimensions) != 0)
              if(length(value) != object@dimensions[["G"]])
                stop("Length of group names does not match the dimensionality of the model")
            # if (length(value) != length(object@data[[1]]))
            #   stop("Group names do not match the number of groups in the training data")
            
            # Define types of nodes
            nodes_types <- .get_nodes_types()

            # Set sample group names in data options
            old_groups <- object@data_options$groups
            object@data_options$groups <- value
            
            # Set sample group names in samples_metadata
            if (!is.null(object@samples_metadata) && (length(object@samples_metadata) != 0)) {
              object@samples_metadata$group_name <- as.character(object@samples_metadata$group_name)
              for (i in 1:object@dimensions[["G"]]) {
                old_name <- old_groups[i]
                new_name <- value[i]
                object@samples_metadata[object@samples_metadata$group_name == old_name, "group_name"] <- new_name
              }
              object@samples_metadata$group_name <- factor(object@samples_metadata$group_name, levels=value)
            }
              
            # Set sample group names in cache
            if (!is.null(object@cache$variance_explained)) {
              names(object@cache$variance_explained$r2_total) <- value
              names(object@cache$variance_explained$r2_per_factor) <- value
            }
          
            # Set sample group names in expectations
            for (node in nodes_types$multigroup_nodes) {
              if (node %in% names(object@expectations)) {
                if (class(object@expectations[[node]])=="list" & length(object@expectations[[node]])==object@dimensions["G"]) {
                  names(object@expectations[[node]]) <- value 
                }
              }
            }
            for (node in nodes_types$twodim_nodes) {
              if (node %in% names(object@expectations)) {
                for (m in 1:length(object@expectations[[node]])) {
                  if (class(object@expectations[[node]][[m]])=="list" & length(object@expectations[[node]][[m]])==object@dimensions["G"]) {
                    names(object@expectations[[node]][[m]]) <- value 
                  }
                }
              }
            }
            
            # Set sample group names in training data
            if (length(object@data)>0) {
              for (m in names(object@data))
                names(object@data[[m]]) <- value
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


# (Hidden) General function to set dimensions names for the data and intercepts
.set_data_names <- function(object, entity, values) {
  
  stopifnot(entity %in% c("features", "samples"))
  
  axes_options <- list(features = 1, samples = 2)
  
  for (m in 1:length(object@data)) {
    for (p in 1:length(object@data[[m]])) {
      deduced_ind <- if (entity == "features") m else p  # since ind corresponds to views (groups of features)
      if (axes_options[[entity]] == 1) {
        rownames(object@data[[m]][[p]]) <- values[[deduced_ind]]
      } else {
        colnames(object@data[[m]][[p]]) <- values[[deduced_ind]]
      }
      
      if (entity=="features")
        tryCatch(names(object@intercepts[[m]][[p]]) <- values[[deduced_ind]], error = function(e) { NULL })
    }
  }
  
  object
}

