
##################
## Factor Names ##
##################

#' @title factors_names: set and retrieve factor names
#' @name factors_names
#' @rdname factors_names
#' @export
setGeneric("factors_names", function(object) { standardGeneric("factors_names") })

#' @name factors_names
#' @rdname factors_names
#' @aliases factors_names<-
#' @export
setGeneric("factors_names<-", function(object, value) { standardGeneric("factors_names<-") })


##################
## Sample Names ##
##################

#' @title samples_names: set and retrieve sample names
#' @name samples_names
#' @rdname samples_names
#' @export
setGeneric("samples_names", function(object) { standardGeneric("samples_names") })

#' @name samples_names
#' @rdname samples_names
#' @aliases samples_names<-
#' @export
setGeneric("samples_names<-", function(object, value) { standardGeneric("samples_names<-") })

##################
## Sample Groups ##
##################

#' @title samples_groups: retrieve sample names and groups
#' @name samples_groups
#' @rdname samples_groups
#' @export
setGeneric("samples_groups", function(object, format = "default") { standardGeneric("samples_groups") })

###################
## Feature Names ##
###################

#' @title features_names: set and retrieve feature names
#' @name features_names
#' @rdname features_names
#' @export
setGeneric("features_names", function(object) { standardGeneric("features_names") })

#' @name features_names
#' @rdname features_names
#' @aliases features_names<-
#' @export
setGeneric("features_names<-", function(object, value) { standardGeneric("features_names<-") })

################
## View Names ##
################

#' @title views_names: set and retrieve view names
#' @name views_names
#' @rdname views_names
#' @export
setGeneric("views_names", function(object) { standardGeneric("views_names") })

#' @name views_names
#' @rdname views_names
#' @aliases views_names<-
#' @export
setGeneric("views_names<-", function(object, value) { standardGeneric("views_names<-") })

################
## group Names ##
################

#' @title groups_names: set and retrieve view names
#' @name groups_names
#' @rdname groups_names
#' @export
setGeneric("groups_names", function(object) { standardGeneric("groups_names") })

#' @name groups_names
#' @rdname groups_names
#' @aliases groups_names<-
#' @export
setGeneric("groups_names<-", function(object, value) { standardGeneric("groups_names<-") })

################
## Input Data ##
################

#' @title Set and retrieve input data
#' @name input_data
#' @export
setGeneric("input_data", function(object) { standardGeneric("input_data") })

#' @name input_data
#' @aliases inputData<-
#' @export
setGeneric(".input_data<-", function(object, value) { standardGeneric(".input_data<-") })


##################
## Imputed Data ##
##################

#' @title imputed_data: set and retrieve imputed data
#' @name imputed_data
#' @export
setGeneric("imputed_data", function(object) { standardGeneric("imputed_data") })

#' @name imputed_data
#' @aliases imputed_data<-
#' @export
setGeneric(".imputed_data<-", function(object, value) { standardGeneric(".imputed_data<-") })

################
## Train Data ##
################

#' @title training_data: set and retrieve training data
#' @name training_data
#' @rdname training_data
#' @export
setGeneric("training_data", function(object) { standardGeneric("training_data") })

#' @name training_data
#' @aliases training_data<-
#' @export
setGeneric(".training_data<-", function(object, value) { standardGeneric(".training_data<-") })


###################
## Train Options ##
###################

#' @title training_options: set and retrieve training opts
#' @name training_options
#' @rdname training_options
#' @export
setGeneric("training_options", function(object) { standardGeneric("training_options") })

#' @name training_options
#' @rdname training_options
#' @aliases training_options<-
#' @export
setGeneric(".training_options<-", function(object, value) { standardGeneric(".training_options<-") })


###################
## Model Options ##
###################

#' @title model_options: set and retrieve Model options
#' @name model_options
#' @export
setGeneric("model_options", function(object) { standardGeneric("model_options") })

#' @name model_options
#' @aliases model_options<-
#' @export
setGeneric(".model_options<-", function(object, value) { standardGeneric(".model_options<-") })


######################
## Train Statistics ##
######################

#' @title training_stats: set and retrieve training statistics
#' @name training_stats
#' @export
setGeneric("training_stats", function(object) { standardGeneric("training_stats") })

#' @name training_stats
#' @aliases training_stats<-
#' @export
setGeneric(".training_stats<-", function(object, value) { standardGeneric(".training_stats<-") })

##################
## Expectations ##
##################

#' @title expectations: set and retrieve expectations
#' @name expectations
#' @rdname expectations
#' @export
setGeneric("expectations", function(object) { standardGeneric("expectations") })

#' @name expectations
#' @rdname expectations
#' @aliases expectations<-
#' @export
setGeneric(".expectations<-", function(object, value) { standardGeneric(".expectations<-") })

# Misc.
setGeneric(".cache_variance_explained<-", function(object, value) { standardGeneric(".cache_variance_explained<-") })
