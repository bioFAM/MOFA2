
# (Hidden) General function to set names
.setNames <- function(object, values, dimensionality, views="all", batches="all") {
  nodes <- names(object@Expectations)

  if (paste0(views, collapse = "") == "all") { 
    views <- names(object@Dimensions$D)
  } else {
    stopifnot(all(views  %in% names(object@Dimensions$D)))
  } 

  if (paste0(batches, collapse = "") == "all") {
    batches <- names(object@Dimensions$N)
  } else {
    stopifnot(all(batches %in% names(object@Dimensions$N)))
  }
  
  # Loop over training data
  for (m in views) {
    for (h in batches) {
      if (nrow(object@TrainData[[m]][[h]]) == dimensionality)
        rownames(object@TrainData[[m]][[h]]) <- values
      if (ncol(object@TrainData[[m]][[h]]) == dimensionality)
        colnames(object@TrainData[[m]][[h]]) <- values
    }
  }
  
  # Loop over nodes
  for (node in nodes) {
    

    if (node %in% c("Y")) {

      for (m in views) {
        for (h in batches) {

          # Loop over expectations
          if (class(object@Expectations[[node]][[m]][[h]]) == "matrix") {
            if (nrow(object@Expectations[[node]][[m]][[h]]) == dimensionality)
              rownames(object@Expectations[[node]][[m]][[h]]) <- values
            if (ncol(object@Expectations[[node]][[m]][[h]]) == dimensionality)
              colnames(object@Expectations[[node]][[m]][[h]]) <- values
          } else if (class(object@Expectations[[node]][[m]][[h]]) == "array") {
            if (length(object@Expectations[[node]][[m]][[h]]) == dimensionality)
              names(object@Expectations[[node]][[m]][[h]]) <- values
          }
          
          # Loop over parameters
          for (Parameter in attributes(object@Parameters[[node]][[m]][[h]])$names){
            
            if (class(object@Parameters[[node]][[m]][[h]][[Parameter]]) == "matrix") {
              if (nrow(object@Parameters[[node]][[m]][[h]][[Parameter]]) == dimensionality)
                rownames(object@Parameters[[node]][[m]][[h]][[Parameter]]) <- values
              if (ncol(object@Parameters[[node]][[m]][[h]][[Parameter]]) == dimensionality)
                colnames(object@Parameters[[node]][[m]][[h]][[Parameter]]) <- values
            } else if (class(object@Parameters[[node]][[m]][[h]][[Parameter]]) %in% c("array",'list')) {
              if (length(object@Parameters[[node]][[m]][[h]][[Parameter]]) == dimensionality)
                names(object@Parameters[[node]][[m]][[h]][[Parameter]]) <- values
            }
            
          }

        }
      }
    
    # Multi-view nodes
    } else if (!(node %in% c("Z", "AlphaZ", "SigmaZ", "ThetaZ"))) {
    #if (node != "Z") {
      
    # if (setequal(names(object@Expectations[[node]]),views)) {
      
      # Loop over views
      for (m in views) {
        
        # Loop over expectations
        if (class(object@Expectations[[node]][[m]]) == "matrix") {
          if (nrow(object@Expectations[[node]][[m]]) == dimensionality)
            rownames(object@Expectations[[node]][[m]]) <- values
          if (ncol(object@Expectations[[node]][[m]]) == dimensionality)
            colnames(object@Expectations[[node]][[m]]) <- values
        } else if (class(object@Expectations[[node]][[m]]) == "array") {
          if (length(object@Expectations[[node]][[m]]) == dimensionality)
            names(object@Expectations[[node]][[m]]) <- values
        }
        
        # Loop over parameters
        for (Parameter in attributes(object@Parameters[[node]][[m]])$names){
          
          if (class(object@Parameters[[node]][[m]][[Parameter]]) == "matrix") {
            if (nrow(object@Parameters[[node]][[m]][[Parameter]]) == dimensionality)
              rownames(object@Parameters[[node]][[m]][[Parameter]]) <- values
            if (ncol(object@Parameters[[node]][[m]][[Parameter]]) == dimensionality)
              colnames(object@Parameters[[node]][[m]][[Parameter]]) <- values
          } else if (class(object@Parameters[[node]][[m]][[Parameter]]) %in% c("array",'list')) {
            if (length(object@Parameters[[node]][[m]][[Parameter]]) == dimensionality)
              names(object@Parameters[[node]][[m]][[Parameter]]) <- values
          }
          
        }
        
      }
      
    # Single-view nodes
    } else {

      # Loop over batches
      for (h in batches) {
        
        # Loop over expectations
        if (class(object@Expectations[[node]][[h]]) == "matrix") {
          if (nrow(object@Expectations[[node]][[h]]) == dimensionality)
            rownames(object@Expectations[[node]][[h]]) <- values
          if (ncol(object@Expectations[[node]][[h]]) == dimensionality)
            colnames(object@Expectations[[node]][[h]]) <- values
        } else if (class(object@Expectations[[node]][[h]]) == "array") {
          if (length(object@Expectations[[node]][[h]]) == dimensionality)
            names(object@Expectations[[node]][[h]]) <- values
        }
        
        # Loop over parameters
        
        for (Parameter in attributes(object@Parameters[[node]][[h]])$names){
          if (class(object@Parameters[[node]][[h]][[Parameter]]) == "matrix") {
            if (nrow(object@Parameters[[node]][[h]][[Parameter]]) == dimensionality)
              rownames(object@Parameters[[node]][[h]][[Parameter]]) <- values
            if (ncol(object@Parameters[[node]][[h]][[Parameter]]) == dimensionality)
              colnames(object@Parameters[[node]][[h]][[Parameter]]) <- values
          } else if (class(object@Parameters[[node]][[h]][[Parameter]]) %in% c("array", 'list')) {
            if (length(object@Parameters[[node]][[h]][[Parameter]]) == dimensionality){
              names(object@Parameters[[node]][[h]][[Parameter]]) <- values
            }
          }
        }

      }
      
    }
    
  }
  
  
  return(object)
}

###################################
## Set and retrieve factor names ##
###################################

#' @rdname factorNames
#' @param object a \code{\link{BioFAM}} object.
#' @aliases factorNames,BioFAModel-method
#' @return character vector with the features names
#' @export
setMethod("factorNames", signature(object="BioFAModel"), function(object) { colnames(object@Expectations$Z[[1]]) })

#' @rdname factorNames
#' @param value a character vector of factor names
#' @import methods
#' @export
setReplaceMethod("factorNames", signature(object="BioFAModel", value="vector"), 
  function(object,value) {
   if (!methods::.hasSlot(object,"Expectations") | length(object@Expectations) == 0)
     stop("Before assigning factor names you have to assign expectations")
   if (methods::.hasSlot(object,"Dimensions") & length(object@Dimensions) != 0)
     if (length(value) != object@Dimensions["K"])
       stop("Length of factor names does not match the dimensionality of the latent variable matrix")
   if(length(value) != ncol(object@Expectations$Z[[1]])) 
     stop("Factor names do not match the number of columns in the latent variable matrix")
   
   object <- .setNames(object, value, object@Dimensions[["K"]])
   object
})



###################################
## Set and retrieve sample names ##
###################################

#' @rdname sampleNames
#' @param object a \code{\link{BioFAModel}} object.
#' @aliases sampleNames,BioFAModel-method
#' @return list of character vectors with the sample names for each batch
#' @export
setMethod("sampleNames", signature(object="BioFAModel"), function(object) { tmp <- lapply(object@TrainData[[1]], colnames); names(tmp) <- batchNames(object); return(tmp) } )

#' @rdname sampleNames
#' @param value list of character vectors with the sample names for every batch
#' @import methods
#' @export
setReplaceMethod("sampleNames", signature(object="BioFAModel", value="list"), 
  function(object,value) {
   if (!methods::.hasSlot(object,"TrainData") | length(object@TrainData) == 0 | length(object@TrainData[[1]]) == 0)
     stop("Before assigning sample names you have to assign the training data")
   if (!methods::.hasSlot(object,"Expectations") | length(object@Expectations) == 0)
     stop("Before assigning sample names you have to assign the expectations")
   if (methods::.hasSlot(object,"Dimensions") & length(object@Dimensions) != 0)
    if (!all(sapply(value,length) == object@Dimensions[["N"]]))
        stop("Length of sample names does not match the dimensionality of the model")
    if (!all(sapply(value,length) == sapply(object@TrainData[[1]], ncol)))
      stop("sample names do not match the dimensionality of the data (columns)")

    for (h in 1:length(object@TrainData[[1]])) {
      object <- .setNames(object, value[[h]], object@Dimensions[["N"]][h], batches = names(object@Dimensions[["N"]][h]))
    }
    object
})

####################################
## Set and retrieve feature names ##
####################################

#' @rdname featureNames
#' @param object a \code{\link{BioFAModel}} object.
#' @aliases featureNames,BioFAModel-method
#' @return list of character vectors with the feature names for each view
#' @export
setMethod("featureNames", signature(object="BioFAModel"), function(object) { tmp <- lapply(object@TrainData, function(e) rownames(e[[1]])); names(tmp) <- viewNames(object); return(tmp) } )

#' @rdname featureNames
#' @param value list of character vectors with the feature names for every view
#' @import methods
#' @export
setReplaceMethod("featureNames", signature(object="BioFAModel", value="list"),
  function(object, value) {
    if (!methods::.hasSlot(object,"TrainData")  | length(object@TrainData) == 0)
      stop("Before assigning feature names you have to assign the training data")
    if (!methods::.hasSlot(object,"Expectations") | length(object@Expectations) == 0)
      stop("Before assigning feature names you have to assign the expectations")
    if (methods::.hasSlot(object,"Dimensions")  | length(object@Dimensions) == 0)
      if (!all(sapply(value, length) == object@Dimensions[["D"]]))
        stop("Length of feature names does not match the dimensionality of the model")
    if (!all(sapply(value, length) == sapply(object@TrainData, function(e) nrow(e[[1]]))))
      stop("Feature names do not match the dimensionality of the data (rows)")
    
    for (m in 1:length(object@TrainData)) {
      object <- .setNames(object, value[[m]], object@Dimensions[["D"]][m], names(object@Dimensions[["D"]][m]))
    }
    return(object)
})

#################################
## Set and retrieve view names ##
#################################

#' @rdname viewNames
#' @param object a \code{\link{BioFAModel}} object.
#' @return character vector with the names for each view
#' @rdname viewNames
#' @export
setMethod("viewNames", signature(object="BioFAModel"), function(object) { names(object@TrainData) } )


#' @rdname viewNames
#' @param value character vector with the names for each view
#' @import methods
#' @export
setMethod("viewNames<-", signature(object="BioFAModel", value="character"), 
  function(object, value) {
    if (!methods::.hasSlot(object,"TrainData") | length(object@TrainData) == 0)
      stop("Before assigning view names you have to assign the training data")
    if (methods::.hasSlot(object,"Dimensions") & length(object@Dimensions) != 0)
      if (length(value) != object@Dimensions["M"])
      stop("Length of view names does not match the dimensionality of the model")
    if (length(value) != length(object@TrainData))
      stop("View names do not match the number of views in the training data")
    
    if (object@Status == "trained") {
      multiview_nodes <- c("W", "AlphaW", "Tau", "Theta", "ThetaW", "Y")
      for (node in multiview_nodes) {
        if (node %in% names(object@Expectations)) {
          if (class(object@Expectations[[node]]) == "list" & 
              length(object@Expectations[[node]]) == object@Dimensions["M"]) {
            names(object@Expectations[[node]]) <- value 
          }
        }
      }
    }
    
    names(object@TrainData) <- value
    names(object@Dimensions$D) <- value
    
    return(object)
})


#################################
## Set and retrieve batch names ##
#################################

#' @rdname batchNames
#' @param object a \code{\link{BioFAModel}} object.
#' @return character vector with the names for each view
#' @rdname batchNames
#' @export
setMethod("batchNames", signature(object="BioFAModel"), function(object) { names(object@TrainData[[1]]) } )


#' @rdname batchNames
#' @param value character vector with the names for each view
#' @import methods
#' @export
setMethod("batchNames<-", signature(object="BioFAModel", value="character"), 
  function(object, value) {
    if (!methods::.hasSlot(object,"TrainData") | length(object@TrainData) == 0)
      stop("Before assigning batch names you have to assign the training data")
    if (methods::.hasSlot(object,"Dimensions") & length(object@Dimensions) != 0)
      if(length(value) != object@Dimensions["H"])
        stop("Length of batch names does not match the dimensionality of the model")
    if (length(value) != length(object@TrainData[[1]]))
      stop("Batch names do not match the number of views in the training data")
    
    if (object@Status == "trained"){
      multibatch_nodes <- c("Z", "Z", "AlphaZ")
      for (node in multibatch_nodes) {
        if (node %in% names(object@Expectations)) {
          if (node == "Y") {  # the only multi-view and multi-batch node
            for (m in 1:length(object@Expectations$Y)) {
              names(object@Expectations$Y[[m]]) <- value 
            }
          } else if (class(object@Expectations[[node]]) == "list" & 
              length(object@Expectations[[node]]) == object@Dimensions["H"]) {
            names(object@Expectations[[node]]) <- value 
          }
        }
      }
    }
    
    for (m in names(object@TrainData)) {
      names(object@TrainData[[m]]) <- value
    }
    names(object@Dimensions$N) <- value
    
    return(object)
})

#################################
## Set and retrieve input data ##
#################################

#' @title Set and retrieve input data
#' @name InputData
#' @param object a \code{\link{BioFAModel}} object.
#' @rdname InputData
#' @export
setMethod("InputData", signature(object="BioFAModel"), function(object) { object@InputData } )

#' @title Set and retrieve input data
#' @docType methods
#' @name InputData
#' @param object a \code{\link{BioFAModel}} object.
#' @rdname InputData
#' @aliases inputData<-
#' @export
setMethod(".InputData<-", signature(object="BioFAModel", value="MultiAssayExperiment"),
          function(object,value) {
            object@InputData <- value
            object
          })

####################################
## Set and retrieve training data ##
####################################

#' @rdname TrainData
#' @param object a \code{\link{BioFAModel}} object.
#' @return list of numeric matrices that contain the training data
#' @rdname TrainData
#' @export
setMethod("TrainData", signature(object="BioFAModel"), function(object) { object@TrainData } )

#' @import methods
setMethod(".TrainData<-", signature(object="BioFAModel", value="list"),
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
    object@TrainData <- value
    object
})

####################################
## Set and retrieve imputed data ##
####################################

#' @rdname ImputedData
#' @param object a \code{\link{BioFAModel}} object.
#' @return list of numeric matrices that contain the training data
#' @rdname ImputedData
#' @export
setMethod("ImputedData", signature(object="BioFAModel"), function(object) { object@ImputedData } )

#' @import methods
setMethod(".ImputedData<-", signature(object="BioFAModel", value="list"),
          function(object,value) {
            # to do sanity checks
            object@ImputedData <- value
            object
          })

#######################################
## Set and retrieve training options ##
#######################################

#' @rdname TrainOpts
#' @param object a \code{\link{BioFAModel}} object.
#' @rdname TrainOpts
#' @return list of training options
#' @export
setMethod("TrainOpts", "BioFAModel", function(object) { object@TrainOpts } )
setMethod(".TrainOpts<-", signature(object="BioFAModel", value="list"),
          function(object,value) {
            object@TrainOpts <- value
            object
          })

#######################################
## Set and retrieve model options ##
#######################################

#' @rdname ModelOpts
#' @param object a \code{\link{BioFAModel}} object.
#' @rdname ModelOpts
#' @return list of model options
#' @export
setMethod("ModelOpts", "BioFAModel", function(object) { object@ModelOpts } )
setMethod(".ModelOpts<-", signature(object="BioFAModel", value="list"),
          function(object,value) {
            object@ModelOpts <- value
            object
          })

##########################################
## Set and retrieve training statistics ##
##########################################

#' @rdname TrainStats
#' @param object a \code{\link{BioFAModel}} object.
#' @return list of training statistics
#' @export
setMethod("TrainStats", "BioFAModel", function(object) { object@TrainStats } )
setMethod(".TrainStats<-", signature(object="BioFAModel", value="list"),
  function(object,value) {
    object@TrainStats <- value
    object
})

###################################
## Set and retrieve expectations ##
###################################

#' @rdname Expectations
#' @param object a \code{\link{BioFAModel}} object.
#' @rdname Expectations
#' @return list of expectations
#' @export
setMethod("Expectations", "BioFAModel", function(object) { object@Expectations } )
setMethod(".Expectations<-", signature(object="BioFAModel", value="list"),
  function(object,value) {
    object@Expectations <- value
    object
})
