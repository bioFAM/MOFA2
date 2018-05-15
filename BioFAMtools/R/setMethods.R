
# (Hidden) General function to set names
# TODO: re-write for expectations as data frames
.set_expectations_names <- function(object, nodes, axes, values) {
  nodes_with_factors <- c("Z", "W", "AlphaZ", "AlphaW", "ThetaZ", "ThetaW")
  axes_for_factors   <- c(2, 2, 1, 1, 1, 1)  # 2 for columns, 1 for rows or vectors
  object <- .set_names(object, nodes, axes, values)

  stopifnot(len(nodes) == len(axes))
  stopifnot(axes[i] == 1 | axes[i] == 2)
  dim <- length(values)

  if (object@status == "trained") {
    for (i in 1:length(nodes)) {
      node <- nodes[i]
      if (node %in% names(object@expectations)) {
        if (class(object@expectations[[node]]) == "list") {
          for (sub in 1:length(object@expectations[[node]])) {
            if (class(object@expectations[[node]][[sub]]) == "list") {  # Y node is a list of lists
              for (sub2 in 1:length(object@expectations[[node]][[sub]])) {
                if (axes[i] == 1) {
                  stopifnot(nrow(object@expectations[[node]][[sub]][[sub2]]) == dim)
                  rownames(object@expectations[[node]][[sub]][[sub2]]) <- values
                } else {
                  stopifnot(ncol(object@expectations[[node]][[sub]][[sub2]]) == dim)
                  colnames(object@expectations[[node]][[sub]][[sub2]]) <- values
                }
              }
            } else {
              if (length(dim(object@expectations[[node]][[sub]])) == 1) {  # set names for a vector
                stopifnot(length(object@expectations[[node]][[sub]]) == dim)
                stopifnot(axes[i] == 1)
                names(object@expectations[[node]][[sub]]) <- values
              } else {  # set names for a matrix
                if (axes[i] == 1) {
                  stopifnot(nrow(object@expectations[[node]][[sub]]) == dim)
                  rownames(object@expectations[[node]][[sub]]) <- values
                } else {
                  stopifnot(ncol(object@expectations[[node]][[sub]]) == dim)
                  colnames(object@expectations[[node]][[sub]]) <- values
                }
              }
            }
          }
        } else {
          if (length(dim(object@expectations[[node]])) == 1) {  # set names for a vector
            stopifnot(length(object@expectations[[node]]) == dim)
            stopifnot(axes[i] == 1)
            names(object@expectations[[node]]) <- values
          } else {  # set names for a matrix
            if (axes[i] == 1) {
              stopifnot(nrow(object@expectations[[node]]) == dim)
              rownames(object@expectations[[node]]) <- values
            } else {
              stopifnot(ncol(object@expectations[[node]]) == dim)
              colnames(object@expectations[[node]]) <- values
            }
          }
        }
      }
    }
  }

  object
}

.set_data_names <- function(object, axes, values) {

}

# .set_names <- function(object, values, dimensionality, views="all", groups="all") {
#   nodes <- names(object@expectations)

#   if (paste0(views, collapse = "") == "all") { 
#     views <- names(object@Dimensions$D)
#   } else {
#     stopifnot(all(views  %in% names(object@Dimensions$D)))
#   } 

#   if (paste0(groups, collapse = "") == "all") {
#     groups <- names(object@Dimensions$N)
#   } else {
#     stopifnot(all(groups %in% names(object@Dimensions$N)))
#   }
  
#   # Loop over training data
#   for (m in views) {
#     for (h in groups) {
#       if (nrow(object@training_data[[m]][[h]]) == dimensionality)
#         rownames(object@training_data[[m]][[h]]) <- values
#       if (ncol(object@training_data[[m]][[h]]) == dimensionality)
#         colnames(object@training_data[[m]][[h]]) <- values
#     }
#   }
  
#   # Loop over nodes
#   for (node in nodes) {
    

#     if (node %in% c("Y")) {

#       for (m in views) {
#         for (h in groups) {

#           # Loop over expectations
#           if (class(object@expectations[[node]][[m]][[h]]) == "matrix") {
#             if (nrow(object@expectations[[node]][[m]][[h]]) == dimensionality)
#               rownames(object@expectations[[node]][[m]][[h]]) <- values
#             if (ncol(object@expectations[[node]][[m]][[h]]) == dimensionality)
#               colnames(object@expectations[[node]][[m]][[h]]) <- values
#           } else if (class(object@expectations[[node]][[m]][[h]]) == "array") {
#             if (length(object@expectations[[node]][[m]][[h]]) == dimensionality)
#               names(object@expectations[[node]][[m]][[h]]) <- values
#           }
          
#           # Loop over parameters
#           for (Parameter in attributes(object@Parameters[[node]][[m]][[h]])$names){
            
#             if (class(object@Parameters[[node]][[m]][[h]][[Parameter]]) == "matrix") {
#               if (nrow(object@Parameters[[node]][[m]][[h]][[Parameter]]) == dimensionality)
#                 rownames(object@Parameters[[node]][[m]][[h]][[Parameter]]) <- values
#               if (ncol(object@Parameters[[node]][[m]][[h]][[Parameter]]) == dimensionality)
#                 colnames(object@Parameters[[node]][[m]][[h]][[Parameter]]) <- values
#             } else if (class(object@Parameters[[node]][[m]][[h]][[Parameter]]) %in% c("array",'list')) {
#               if (length(object@Parameters[[node]][[m]][[h]][[Parameter]]) == dimensionality)
#                 names(object@Parameters[[node]][[m]][[h]][[Parameter]]) <- values
#             }
            
#           }

#         }
#       }
    
#     # Multi-view nodes
#     } else if (!(node %in% c("Z", "AlphaZ", "SigmaZ", "ThetaZ"))) {
#     #if (node != "Z") {
      
#     # if (setequal(names(object@expectations[[node]]),views)) {
      
#       # Loop over views
#       for (m in views) {
        
#         # Loop over expectations
#         if (class(object@expectations[[node]][[m]]) == "matrix") {
#           if (nrow(object@expectations[[node]][[m]]) == dimensionality)
#             rownames(object@expectations[[node]][[m]]) <- values
#           if (ncol(object@expectations[[node]][[m]]) == dimensionality)
#             colnames(object@expectations[[node]][[m]]) <- values
#         } else if (class(object@expectations[[node]][[m]]) == "array") {
#           if (length(object@expectations[[node]][[m]]) == dimensionality)
#             names(object@expectations[[node]][[m]]) <- values
#         }
        
#         # Loop over parameters
#         for (Parameter in attributes(object@Parameters[[node]][[m]])$names){
          
#           if (class(object@Parameters[[node]][[m]][[Parameter]]) == "matrix") {
#             if (nrow(object@Parameters[[node]][[m]][[Parameter]]) == dimensionality)
#               rownames(object@Parameters[[node]][[m]][[Parameter]]) <- values
#             if (ncol(object@Parameters[[node]][[m]][[Parameter]]) == dimensionality)
#               colnames(object@Parameters[[node]][[m]][[Parameter]]) <- values
#           } else if (class(object@Parameters[[node]][[m]][[Parameter]]) %in% c("array",'list')) {
#             if (length(object@Parameters[[node]][[m]][[Parameter]]) == dimensionality)
#               names(object@Parameters[[node]][[m]][[Parameter]]) <- values
#           }
          
#         }
        
#       }
      
#     # Single-view nodes
#     } else {

#       # Loop over groups
#       for (h in groups) {
        
#         # Loop over expectations
#         if (class(object@expectations[[node]][[h]]) == "matrix") {
#           if (nrow(object@expectations[[node]][[h]]) == dimensionality)
#             rownames(object@expectations[[node]][[h]]) <- values
#           if (ncol(object@expectations[[node]][[h]]) == dimensionality)
#             colnames(object@expectations[[node]][[h]]) <- values
#         } else if (class(object@expectations[[node]][[h]]) == "array") {
#           if (length(object@expectations[[node]][[h]]) == dimensionality)
#             names(object@expectations[[node]][[h]]) <- values
#         }
        
#         # Loop over parameters
        
#         for (Parameter in attributes(object@Parameters[[node]][[h]])$names){
#           if (class(object@Parameters[[node]][[h]][[Parameter]]) == "matrix") {
#             if (nrow(object@Parameters[[node]][[h]][[Parameter]]) == dimensionality)
#               rownames(object@Parameters[[node]][[h]][[Parameter]]) <- values
#             if (ncol(object@Parameters[[node]][[h]][[Parameter]]) == dimensionality)
#               colnames(object@Parameters[[node]][[h]][[Parameter]]) <- values
#           } else if (class(object@Parameters[[node]][[h]][[Parameter]]) %in% c("array", 'list')) {
#             if (length(object@Parameters[[node]][[h]][[Parameter]]) == dimensionality){
#               names(object@Parameters[[node]][[h]][[Parameter]]) <- values
#             }
#           }
#         }

#       }
      
#     }
    
#   }
  
  
#   return(object)
# }

###################################
## Set and retrieve factor names ##
###################################

#' @rdname factors_names
#' @param object a \code{\link{BioFAM}} object.
#' @aliases factors_names,BioFAModel-method
#' @return character vector with the features names
#' @export
setMethod("factors_names", signature(object="BioFAModel"), function(object) { colnames(object@expectations$Z[[1]]) })

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
    if(length(value) != ncol(object@expectations$Z[[1]])) 
      stop("Factor names do not match the number of columns in the latent variable matrix")
 
    nodes_with_factors <- c("Z", "W", "AlphaZ", "AlphaW", "ThetaZ", "ThetaW")
    axes_for_factors   <- c(2, 2, 1, 1, 1, 1)
    object <- .set_expectations_names(object, nodes_with_factors, axes_for_factors, value)
 
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
setMethod("samples_names", signature(object="BioFAModel"), function(object) { tmp <- lapply(object@training_data[[1]], colnames); names(tmp) <- groups_names(object); return(tmp) } )

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

    nodes_with_samples <- c("Z", "Tau", "Y")
    axes_for_samples   <- c(1, 1, 1)
    object <- .set_expectations_names(object, nodes_with_samples, axes_for_samples, value)


    for (h in 1:length(object@training_data[[1]])) {
      object <- .set_names(object, value[[h]], object@dimensions[["N"]][h], groups = names(object@dimensions[["N"]][h]))
    }
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
setMethod("features_names", signature(object="BioFAModel"), function(object) { tmp <- lapply(object@training_data, function(e) rownames(e[[1]])); names(tmp) <- views_names(object); return(tmp) } )

#' @rdname features_names
#' @param value list of character vectors with the feature names for every view
#' @import methods
#' @export
setReplaceMethod("features_names", signature(object="BioFAModel", value="list"),
  function(object, value) {
    if (!methods::.hasSlot(object,"training_data")  | length(object@training_data) == 0)
      stop("Before assigning feature names you have to assign the training data")
    if (!methods::.hasSlot(object,"expectations") | length(object@expectations) == 0)
      stop("Before assigning feature names you have to assign the expectations")
    if (methods::.hasSlot(object,"Dimensions")  | length(object@Dimensions) == 0)
      if (!all(sapply(value, length) == object@Dimensions[["D"]]))
        stop("Length of feature names does not match the dimensionality of the model")
    if (!all(sapply(value, length) == sapply(object@training_data, function(e) nrow(e[[1]]))))
      stop("Feature names do not match the dimensionality of the data (rows)")
    
    for (m in 1:length(object@training_data)) {
      object <- .setNames(object, value[[m]], object@dimensions[["D"]][m], names(object@dimensions[["D"]][m]))
    }
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
setMethod("views_names", signature(object="BioFAModel"), function(object) { names(object@training_data) } )


#' @rdname views_names
#' @param value character vector with the names for each view
#' @import methods
#' @export
setMethod("views_names<-", signature(object="BioFAModel", value="character"), 
  function(object, value) {
    if (!methods::.hasSlot(object, "training_data") | length(object@training_data) == 0)
      stop("Before assigning view names you have to assign the training data")
    if (methods::.hasSlot(object, "dimensions") & length(object@dimensions) != 0)
      if (length(value) != object@dimensions["M"])
      stop("Length of view names does not match the dimensionality of the model")
    if (length(value) != length(object@training_data))
      stop("View names do not match the number of views in the training data")
    
    if (object@status == "trained") {
      multiview_nodes <- c("W", "AlphaW", "Tau", "Theta", "ThetaW", "Y")
      for (node in multiview_nodes) {
        if (node %in% names(object@expectations)) {
          if (class(object@expectations[[node]]) == "list" & 
              length(object@expectations[[node]]) == object@dimensions["M"]) {
            names(object@expectations[[node]]) <- value 
          }
        }
      }
    }
    
    names(object@training_data) <- value
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
setMethod("groups_names", signature(object="BioFAModel"), function(object) { names(object@training_data[[1]]) } )


#' @rdname groups_names
#' @param value character vector with the names for each view
#' @import methods
#' @export
setMethod("groups_names<-", signature(object="BioFAModel", value="character"), 
  function(object, value) {
    if (!methods::.hasSlot(object,"training_data") | length(object@training_data) == 0)
      stop("Before assigning group names you have to assign the training data")
    if (methods::.hasSlot(object,"Dimensions") & length(object@Dimensions) != 0)
      if(length(value) != object@Dimensions["P"])
        stop("Length of group names does not match the dimensionality of the model")
    if (length(value) != length(object@training_data[[1]]))
      stop("Group names do not match the number of groups in the training data")
    
    if (object@Status == "trained"){
      multigroup_nodes <- c("Z", "Y", "AlphaZ")
      for (node in multigroup_nodes) {
        if (node %in% names(object@expectations)) {
          if (node == "Y") {  # the only multi-view and multi-group node
            for (m in 1:length(object@expectations$Y)) {
              names(object@expectations$Y[[m]]) <- value 
            }
          } else if (class(object@expectations[[node]]) == "list" & 
              length(object@expectations[[node]]) == object@Dimensions["P"]) {
            names(object@expectations[[node]]) <- value 
          }
        }
      }
    }
    
    for (m in names(object@training_data)) {
      names(object@training_data[[m]]) <- value
    }
    names(object@Dimensions$N) <- value
    
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
setMethod("input_data", signature(object="BioFAModel"), function(object) { object@input_data } )

#' @title Set and retrieve input data
#' @docType methods
#' @name input_data
#' @param object a \code{\link{BioFAModel}} object.
#' @rdname input_data
#' @aliases inputData<-
#' @export
setMethod(".input_data<-", signature(object="BioFAModel", value="MultiAssayExperiment"),
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
setMethod("training_data", signature(object="BioFAModel"), function(object) { object@training_data } )

#' @import methods
setMethod(".training_data<-", signature(object="BioFAModel", value="list"),
  function(object,value) {
    # N <- unique(sapply(value,nrow))
    # if (length(N) > 1) 
    #   stop("Views do not have the same number of samples (rows)")
    # if (methods::.hasSlot(object,"Dimensions")) {
    #   if (object@Dimensions["M"] != length(value))
    #     if (object@Dimensions["N"] != N)
    #       stop("Number of samples in the data do not match the specified dimensionality of the model")
    #   if (all(object@Dimensions["D"] != sapply(value,ncol)))
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
