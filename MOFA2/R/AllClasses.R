
##########################################################
## Define a general class to store a MOFA trained model ##
##########################################################

#' @title Class to store a mofa model
#' @description
#' The \code{MOFA} is an S4 class used to store all relevant data to analyse a MOFA model
#' @section Slots:
#'  \itemize{
#'    \item{\code{data}:}{ The input data. }
#'    \item{\code{intercepts}:}{ Feature intercepts. }
#'    \item{\code{samples_metadata}:}{ Samples metadata. }
#'    \item{\code{features_metadata}:}{ Features metadata. }
#'    \item{\code{imputed_data}:}{ The imputed data. }
#'    \item{\code{expectations}:}{ expected values of the factors and the loadings. }
#'    \item{\code{dim_red}:}{ non-linear dimensionality reduction manifolds. }
#'    \item{\code{training_stats}:}{ model training statistics. }
#'    \item{\code{data_options}:}{ Data processing options. }
#'    \item{\code{training_options}:}{ Model training options. }
#'    \item{\code{stochastic_options}:}{ Stochastic variational inference options. }
#'    \item{\code{model_options}:}{ Model options. }
#'    \item{\code{dimensions}:}{ Dimensionalities of the model: 
#'    M for the number of views, 
#'    G for the number of groups, 
#'    N for the number of samples (per group), 
#'    D for the number of features (per view),
#'    K for the number of factors.}
#'    \item{\code{on_disk}:}{ Logical indicating whether data is loaded from disk. }
#'    \item{\code{cache}:}{ Cache.}
#'    \item{\code{status}:}{ Auxiliary variable indicating whether the model has been trained.}
#'}
#' @name MOFA
#' @rdname MOFA
#' @aliases MOFA-class
#' @exportClass MOFA
setClass("MOFA", 
        slots=c(
            data                = "list",
            intercepts          = "list",
            imputed_data        = "list",
            samples_metadata    = "list",
            features_metadata   = "list",
            expectations        = "list", 
            training_stats      = "list",
            training_options    = "list",
            stochastic_options  = "list",
            data_options        = "list",
            model_options       = "list",
            dimensions          = "list",
            on_disk             = "logical",
            dim_red             = "list",
            cache               = "list",
            status              = "character"
        )
)

# Printing method
setMethod("show", "MOFA", function(object) {
  
  if (!.hasSlot(object, "dimensions") | length(object@dimensions) == 0)
    stop("Error: dimensions not defined")
  if (!.hasSlot(object, "status") | length(object@status) == 0)
    stop("Error: status not defined")
  
  if (object@status == "trained") {
    nfactors <- object@dimensions[["K"]]
    cat(sprintf("Trained MOFA with the following characteristics: \n Number of views: %d \n Views names: %s \n Number of features (per view): %s \n Number of groups: %d \n Groups names: %s \n Number of samples (per group): %s \n Number of factors: %d \n",
                object@dimensions[["M"]], paste(views_names(object),  collapse=" "), paste(as.character(object@dimensions[["D"]]), collapse=" "), 
                object@dimensions[["G"]], paste(groups_names(object), collapse=" "), paste(as.character(object@dimensions[["N"]]), collapse=" "), 
                nfactors))
  } else {
    cat(sprintf("Untrained MOFA model with the following characteristics: \n Number of views: %d \n Views names: %s \n Number of features (per view): %s \n Number of groups: %d \n Groups names: %s \n Number of samples (per group): %s \n",
                object@dimensions[["M"]], paste(views_names(object),  collapse=" "), paste(as.character(object@dimensions[["D"]]), collapse=" "),
                object@dimensions[["G"]], paste(groups_names(object), collapse=" "), paste(as.character(object@dimensions[["N"]]), collapse=" ")))
  }
  cat("\n")
})


