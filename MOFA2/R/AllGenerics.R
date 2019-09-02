
##################
## Factor Names ##
##################

#' @title factors: set and retrieve factor names
#' @name factors
#' @rdname factors
#' @export
setGeneric("factors", function(object) { standardGeneric("factors") })

#' @name factors
#' @rdname factors
#' @aliases factors<-
#' @export
setGeneric("factors<-", function(object, value) { standardGeneric("factors<-") })


##################
## Sample Names ##
##################

#' @title samples: set and retrieve sample names
#' @name samples
#' @rdname samples
#' @export
setGeneric("samples", function(object) { standardGeneric("samples") })

#' @name samples
#' @rdname samples
#' @aliases samples<-
#' @export
setGeneric("samples<-", function(object, value) { standardGeneric("samples<-") })

#####################
## Sample Metadata ##
#####################

#' @title samples_metadata: retrieve sample metadata
#' @name samples_metadata
#' @rdname samples_metadata
#' @export
setGeneric("samples_metadata", function(object, format = "default") { standardGeneric("samples_metadata") })

#' @name samples_metadata
#' @rdname samples_metadata
#' @aliases samples_metadata<-
#' @export
setGeneric("samples_metadata<-", function(object, value) { standardGeneric("samples_metadata<-") })

###################
## Feature Names ##
###################

#' @title features: set and retrieve feature names
#' @name features
#' @rdname features
#' @export
setGeneric("features", function(object) { standardGeneric("features") })

#' @name features
#' @rdname features
#' @aliases features<-
#' @export
setGeneric("features<-", function(object, value) { standardGeneric("features<-") })

######################
## Feature Metadata ##
######################

#' @title features_metadata: set and retrieve feature metadata
#' @name features_metadata
#' @rdname features_metadata
#' @export
setGeneric("features_metadata", function(object, format = "default") { standardGeneric("features_metadata") })

#' @name features_metadata
#' @rdname features_metadata
#' @aliases features_metadata<-
#' @export
setGeneric("features_metadata<-", function(object, value) { standardGeneric("features_metadata<-") })

################
## View Names ##
################

#' @title views: set and retrieve view names
#' @name views
#' @rdname views
#' @export
setGeneric("views", function(object) { standardGeneric("views") })

#' @name views
#' @rdname views
#' @aliases views<-
#' @export
setGeneric("views<-", function(object, value) { standardGeneric("views<-") })

################
## group Names ##
################

#' @title groups: set and retrieve view names
#' @name groups
#' @rdname groups
#' @export
setGeneric("groups", function(object) { standardGeneric("groups") })

#' @name groups
#' @rdname groups
#' @aliases groups<-
#' @export
setGeneric("groups<-", function(object, value) { standardGeneric("groups<-") })

################
## Input Data ##
################

#' @title Retrieve input data
#' @name data
#' @export
setGeneric("data", function(object) { standardGeneric("data") })


##################
## Imputed Data ##
##################

#' @title imputed_data: set and retrieve imputed data
#' @name imputed_data
#' @export
setGeneric("imputed_data", function(object) { standardGeneric("imputed_data") })

###################
## Train Options ##
###################

#' @title training_options: set and retrieve training opts
#' @name training_options
#' @rdname training_options
#' @export
setGeneric("training_options", function(object) { standardGeneric("training_options") })


###################
## Model Options ##
###################

#' @title model_options: set and retrieve model options
#' @name model_options
#' @export
setGeneric("model_options", function(object) { standardGeneric("model_options") })


######################
## Train Statistics ##
######################

#' @title training_stats: set and retrieve training statistics
#' @name training_stats
#' @export
setGeneric("training_stats", function(object) { standardGeneric("training_stats") })

##################
## Expectations ##
##################

#' @title expectations: set and retrieve expectations
#' @name expectations
#' @rdname expectations
#' @export
setGeneric("expectations", function(object) { standardGeneric("expectations") })
