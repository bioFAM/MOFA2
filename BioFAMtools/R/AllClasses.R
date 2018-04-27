
##########################################################
## Define a general class to store a MOFA trained model ##
##########################################################

#' @title Class to store a Multi-Omics Factor Analysis (MOFA) model
#' @description
#' The \code{BioFAModel} is an S4 class used to store all relevant data to analyse a MOFA model.
#' @section Slots:
#'  \itemize{
#'    \item{\code{InputData}:}{ the input data before being parsed to Training Data. 
#'    Either a MultiAssayExperiment object or a list of matrices, one per view.}
#'    \item{\code{TrainData}:}{ the parsed data used to fit the MOFA model. 
#'    A list with one matrix per view.}
#'    \item{\code{ImputedData}:}{ the parsed data with the missing values imputed using the MOFA model. 
#'    A list with one matrix per view.}
#'    \item{\code{Expectations}:}{ expected values of the different variables of the model. }
#'    A list of matrices, one per variable. The most relevant are "W" for weights and "Z" for factors.
#'    \item{\code{TrainStats}:}{ list with training statistics such as evidence lower bound (ELBO), number of active factors, etc.}
#'    \item{\code{DataOpts}:}{ list with the data processing options such as whether to center or scale the data.}
#'    \item{\code{TrainOpts}:}{ list with the training options such as maximum number of iterations, tolerance for convergence, etc.}
#'    \item{\code{ModelOpts}:}{ list with the model options such as likelihoods, whether an intercept factors was learnt, etc.}
#'    \item{\code{Dimensions}:}{ list with the relevant dimensionalities of the model. N for the number of samples, 
#'    M for the number of views, D for the number of features of each view and K for the number of infered latent factors.}
#'    \item{\code{Status}:}{Auxiliary variable indicating whether the model has been trained.}
#'}
#' @name BioFAModel
#' @rdname BioFAModel
#' @aliases BioFAModel-class
#' @exportClass BioFAModel
setClass("BioFAModel", 
         slots=c(InputData = "MultiAssayExperiment",
                 TrainData = "list",
                 ImputedData = "list",
                 Parameters = "list",
                 Expectations = "list", 
                 TrainStats = "list",
                 TrainOpts = "list",
                 DataOpts = "list",
                 ModelOpts = "list",
                 Dimensions = "list",
                 Status = "character")
)

# Printing method
setMethod("show", "BioFAModel", function(object) {
  
  if(!.hasSlot(object,"Dimensions") | length(object@Dimensions) == 0)
    stop("Error: Dimensions not defined")
  if(!.hasSlot(object,"Status") | length(object@Status) == 0)
    stop("Error: Status not defined")
  
  if (object@Status == "trained") {
    nfactors <- object@Dimensions[["K"]]
    if (object@ModelOpts$learnIntercept==TRUE) { nfactors <- nfactors-1 }
    cat(sprintf("Trained BioFAModel with the following characteristics: \n Number of views: %d \n View names: %s \n Number of features per view: %s \n Number of samples: %d \n Number of factors: %d \n",
                object@Dimensions[["M"]], paste(viewNames(object),collapse=" "), paste(as.character(object@Dimensions[["D"]]),collapse=" "), object@Dimensions[["N"]], nfactors))
  } else {
    cat(sprintf("Untrained BioFAModel model with the following characteristics: \n Number of views: %d \n View names: %s \n Number of features per view: %s \n Number of samples: %d \n",
                object@Dimensions[["M"]], paste(viewNames(object),collapse=" "), paste(as.character(object@Dimensions[["D"]]),collapse=" "), object@Dimensions[["N"]]))
    cat("\n")
  }
})
