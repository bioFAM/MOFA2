
##########################################################
## Define a general class to store a MOFA trained model ##
##########################################################

#' @title Class to store a model with biofam
#' @description
#' The \code{BioFAModel} is an S4 class used to store all relevant data to analyse a BioFAModel
#' @section Slots:
#'  \itemize{
#'    \item{\code{input_data}:}{ the input data before being parsed to Training Data. 
#'    Either a MultiAssayExperiment object or a list of matrices, one per view.}
#'    \item{\code{training_data}:}{ the parsed data used to fit the BioFAModel 
#'    A list with one matrix per view.}
#'    \item{\code{imputed_data}:}{ the parsed data with the missing values imputed using the BioFAModel 
#'    A list with one matrix per view.}
#'    \item{\code{expectations}:}{ expected values of the different variables of the model. }
#'    A list of matrices, one per variable. The most relevant are "W" for weights and "Z" for factors.
#'    \item{\code{training_stats}:}{ list with training statistics such as evidence lower bound (ELBO), number of active factors, etc.}
#'    \item{\code{data_options}:}{ list with the data processing options such as whether to center or scale the data.}
#'    \item{\code{training_options}:}{ list with the training options such as maximum number of iterations, tolerance for convergence, etc.}
#'    \item{\code{model_options}:}{ list with the model options such as likelihoods, whether an intercept factors was learnt, etc.}
#'    \item{\code{dimensions}:}{ list with the relevant dimensionalities of the model. N for the number of samples, 
#'    M for the number of views, D for the number of features of each view and K for the number of infered latent factors.}
#'    \item{\code{status}:}{Auxiliary variable indicating whether the model has been trained.}
#'}
#' @name BioFAModel
#' @rdname BioFAModel
#' @aliases BioFAModel-class
#' @exportClass BioFAModel
setClass("BioFAModel", 
         slots=c(input_data       = "MultiAssayExperiment",
                 training_data    = "list",
                 imputed_data     = "list",
                 parameters       = "list",
                 expectations     = "list", 
                 training_stats   = "list",
                 training_options = "list",
                 data_options     = "list",
                 model_options    = "list",
                 dimensions       = "list",
                 status           = "character")
)

# Printing method
setMethod("show", "BioFAModel", function(object) {
  
  if(!.hasSlot(object, "dimensions") | length(object@dimensions) == 0)
    stop("Error: dimensions not defined")
  if(!.hasSlot(object, "status") | length(object@status) == 0)
    stop("Error: status not defined")
  
  if (object@status == "trained") {
    nfactors <- object@dimensions[["K"]]
    if (object@model_options$learn_intercept) { nfactors <- nfactors - 1 }
    cat(sprintf("Trained BioFAModel with the following characteristics: \n Number of views: %d \n View names: %s \n Number of features per view: %s \n Number of groups: %d \n Group names: %s \n Number of samples per group: %s \n Number of factors: %d \n",
                object@dimensions[["M"]], paste(views_names(object),  collapse=" "), paste(as.character(object@dimensions[["D"]]), collapse=" "), 
                object@dimensions[["H"]], paste(groups_names(object), collapse=" "), paste(as.character(object@dimensions[["N"]]), collapse=" "), 
                nfactors))
  } else {
    cat(sprintf("Untrained BioFAModel model with the following characteristics: \n Number of views: %d \n View names: %s \n Number of features per view: %s \n Number of samples: %d \n",
                object@dimensions[["M"]], paste(views_names(object),  collapse=" "), paste(as.character(object@dimensions[["D"]]), collapse=" "),
                object@dimensions[["H"]], paste(groups_names(object), collapse=" "), paste(as.character(object@dimensions[["N"]]), collapse=" ")))
  }
  cat("\n")
})
