
##################
## Factor Names ##
##################

#' @title factorNames: set and retrieve factor names
#' @name factorNames
#' @rdname factorNames
#' @export
setGeneric("factorNames", function(object) {standardGeneric("factorNames")})

#' @name factorNames
#' @rdname factorNames
#' @aliases factorNames<-
#' @export
setGeneric("factorNames<-", function(object, value) {standardGeneric("factorNames<-")})


##################
## Sample Names ##
##################

#' @title sampleNames: set and retrieve sample names
#' @name sampleNames
#' @rdname sampleNames
#' @export
setGeneric("sampleNames", function(object) {standardGeneric("sampleNames")})

#' @name sampleNames
#' @rdname sampleNames
#' @aliases sampleNames<-
#' @export
setGeneric("sampleNames<-", function(object, value) {standardGeneric("sampleNames<-")})

###################
## Feature Names ##
###################

#' @title featureNames: set and retrieve feature names
#' @name featureNames
#' @rdname featureNames
#' @export
setGeneric("featureNames", function(object) {standardGeneric("featureNames")})

#' @name featureNames
#' @rdname featureNames
#' @aliases featureNames<-
#' @export
setGeneric("featureNames<-", function(object, value) {standardGeneric("featureNames<-")})

################
## View Names ##
################

#' @title viewNames: set and retrieve view names
#' @name viewNames
#' @rdname viewNames
#' @export
setGeneric("viewNames", function(object) {standardGeneric("viewNames")})

#' @name viewNames
#' @rdname viewNames
#' @aliases viewNames<-
#' @export
setGeneric("viewNames<-", function(object, value) {standardGeneric("viewNames<-")})

################
## Batch Names ##
################

#' @title batchNames: set and retrieve view names
#' @name batchNames
#' @rdname batchNames
#' @export
setGeneric("batchNames", function(object) {standardGeneric("batchNames")})

#' @name batchNames
#' @rdname batchNames
#' @aliases batchNames<-
#' @export
setGeneric("batchNames<-", function(object, value) {standardGeneric("batchNames<-")})

################
## Input Data ##
################

#' @title Set and retrieve input data
#' @name InputData
#' @export
setGeneric("InputData", function(object) {standardGeneric("InputData")})

#' @name InputData
#' @aliases inputData<-
#' @export
setGeneric(".InputData<-", function(object, value) {standardGeneric(".InputData<-")})


##################
## Imputed Data ##
##################

#' @title ImputedData: set and retrieve imputed data
#' @name ImputedData
#' @export
setGeneric("ImputedData", function(object) {standardGeneric("ImputedData")})

#' @name ImputedData
#' @aliases ImputedData<-
#' @export
setGeneric(".ImputedData<-", function(object, value) {standardGeneric(".ImputedData<-")})

################
## Train Data ##
################

#' @title TrainData: set and retrieve training data
#' @name TrainData
#' @rdname TrainData
#' @export
setGeneric("TrainData", function(object) {standardGeneric("TrainData")})

#' @name TrainData
#' @aliases TrainData<-
#' @export
setGeneric(".TrainData<-", function(object, value) {standardGeneric(".TrainData<-")})


###################
## Train Options ##
###################

#' @title TrainOptions: set and retrieve training opts
#' @name TrainOptions
#' @rdname TrainOptions
#' @export
setGeneric("TrainOptions", function(object) {standardGeneric("TrainOptions")})

#' @name TrainOptions
#' @rdname TrainOptions
#' @aliases TrainOptions<-
#' @export
setGeneric(".TrainOptions<-", function(object, value) {standardGeneric(".TrainOptions<-")})


###################
## Model Options ##
###################

#' @title ModelOptions: set and retrieve Model options
#' @name ModelOptions
#' @export
setGeneric("ModelOptions", function(object) {standardGeneric("ModelOptions")})

#' @name ModelOptions
#' @aliases ModelOptions<-
#' @export
setGeneric(".ModelOptions<-", function(object, value) {standardGeneric(".ModelOptions<-")})


######################
## Train Statistics ##
######################

#' @title TrainStats: set and retrieve training statistics
#' @name TrainStats
#' @export
setGeneric("TrainStats", function(object) {standardGeneric("TrainStats")})

#' @name TrainStats
#' @aliases TrainStats<-
#' @export
setGeneric(".TrainStats<-", function(object, value) {standardGeneric(".TrainStats<-")})

##################
## Expectations ##
##################

#' @title Expectations: set and retrieve expectations
#' @name Expectations
#' @rdname Expectations
#' @export
setGeneric("Expectations", function(object) {standardGeneric("Expectations")})

#' @name Expectations
#' @rdname Expectations
#' @aliases Expectations<-
#' @export
setGeneric(".Expectations<-", function(object, value) {standardGeneric(".Expectations<-")})

