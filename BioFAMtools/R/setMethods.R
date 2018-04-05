
# (Hidden) General function to set names
.setNames <- function(object, values, dimensionality, views="all") {
  nodes <- names(object@Expectations)
  if (paste0(views,collapse="") == "all") { 
    views <- names(object@Dimensions$D) 
  } else {
    stopifnot(all(views%in%names(object@Dimensions$D) ))
  } 
  
  # Loop over training data
  for (m in views) {
    if (nrow(object@TrainData[[m]]) == dimensionality)
      rownames(object@TrainData[[m]]) <- values
    if (ncol(object@TrainData[[m]]) == dimensionality)
      colnames(object@TrainData[[m]]) <- values
  }
  
  
  # Loop over nodes
  for (node in nodes) {
    
    # Multi-view nodes
    if (node != "Z") {
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
        
      }
      
    # Single-view nodes
    } else {
      
      # Loop over expectations
        if (class(object@Expectations[[node]]) == "matrix") {
          if (nrow(object@Expectations[[node]]) == dimensionality)
            rownames(object@Expectations[[node]]) <- values
          if (ncol(object@Expectations[[node]]) == dimensionality)
            colnames(object@Expectations[[node]]) <- values
        } else if (class(object@Expectations[[node]]) == "array") {
          if (length(object@Expectations[[node]]) == dimensionality)
            names(object@Expectations[[node]]) <- values
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
setMethod("factorNames", signature(object="BioFAModel"), function(object) { colnames(object@Expectations$Z) } )

#' @rdname factorNames
#' @param value a character vector of factor names
#' @import methods
#' @export
setReplaceMethod("factorNames", signature(object="BioFAModel", value="vector"), 
  function(object,value) {
   if (!methods::.hasSlot(object,"Expectations") | length(object@Expectations) == 0)
     stop("Before assigning factor names you have to assign expectations")
   if (methods::.hasSlot(object,"Dimensions") | length(object@Dimensions) == 0)
     if (!length(value)==object@Dimensions["K"])
       stop("Length of factor names does not match the dimensionality of the latent variable matrix")
   if(!length(value)==ncol(object@Expectations$Z)) 
     stop("factor names do not match the number of columns in the latent variable matrix")
    
   object <- .setNames(object, value, object@Dimensions[["K"]])
   object
})



###################################
## Set and retrieve sample names ##
###################################

#' @rdname sampleNames
#' @param object a \code{\link{BioFAModel}} object.
#' @aliases sampleNames,BioFAModel-method
#' @return character vector with the sample names
#' @export
setMethod("sampleNames", signature(object="BioFAModel"), function(object) { colnames(object@TrainData[[1]]) } )

#' @rdname sampleNames
#' @param value a character vector of sample names
#' @import methods
#' @export
setReplaceMethod("sampleNames", signature(object="BioFAModel", value="vector"), 
  function(object,value) {
   if (!methods::.hasSlot(object,"TrainData") | length(object@TrainData) == 0)
     stop("Before assigning sample names you have to assign the training data")
   if (!methods::.hasSlot(object,"Expectations") | length(object@Expectations) == 0)
     stop("Before assigning sample names you have to assign the expectations")
   if (methods::.hasSlot(object,"Dimensions") | length(object@Dimensions) == 0)
     if (!length(value)==object@Dimensions["N"])
       stop("Length of sample names does not match the dimensionality of the model")
   if(!length(value)==ncol(object@TrainData[[1]])) 
     stop("sample names do not match the dimensionality of the data (cols) ")
   
    object <- .setNames(object, value, object@Dimensions[["N"]])
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
setMethod("featureNames", signature(object="BioFAModel"), function(object) { tmp <- lapply(object@TrainData,rownames); names(tmp) <- viewNames(object); return(tmp) } )

#' @rdname featureNames
#' @param value list of character vectors with the feature names for each view
#' @import methods
#' @export
setReplaceMethod("featureNames", signature(object="BioFAModel", value="list"), 
  function(object,value) {
    if (!methods::.hasSlot(object,"TrainData")  | length(object@TrainData) == 0)
      stop("Before assigning feature names you have to assign the training data")
    if (!methods::.hasSlot(object,"Expectations")  | length(object@Expectations) == 0)
      stop("Before assigning feature names you have to assign the expectations")
    if (methods::.hasSlot(object,"Dimensions")  | length(object@Dimensions) == 0)
      if (!all(sapply(value,length) == object@Dimensions[["D"]]))
        stop("Length of feature names does not match the dimensionality of the model")
    if (!all(sapply(value,length)==sapply(object@TrainData,nrow)))
      stop("feature names do not match the dimensionality of the data (columns)")
    
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
  function(object,value) {
    if (!methods::.hasSlot(object,"TrainData") | length(object@TrainData) == 0)
      stop("Before assigning view names you have to assign the training data")
    if (methods::.hasSlot(object,"Dimensions")| length(object@Dimensions) == 0)
      if (!length(value) == object@Dimensions["M"])
        stop("Length of view names does not match the dimensionality of the model")
    if (!length(value)==length(object@TrainData))
      stop("view names do not match the number of views in the training data")
    
    # We have to modify this
    if (object@Status == "trained"){
      multiview_nodes <- c("AlphaW","W","Tau","Theta","Y")
      for (node in multiview_nodes) { 
        names(object@Expectations[[node]]) <- value 
      }
    }
    
    names(object@TrainData) <- value
    names(object@Dimensions$D) <- value
    
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
