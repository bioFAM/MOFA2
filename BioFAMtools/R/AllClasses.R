
##########################################################
## Define a general class to store a MOFA trained model ##
##########################################################

#' @title Class to store a biofam model
#' @description
#' The \code{BioFAModel} is an S4 class used to store all relevant data to analyse a BioFAModel
#' @section Slots:
#'  \itemize{
#'    \item{\code{input_data}:}{ TO SPECIFY, the input data before being parsed to Training Data. 
#'    }
#'    \item{\code{training_data}:}{ TO SPECIFY the parsed data used to fit the BioFAModel 
#'    }
#'    \item{\code{imputed_data}:}{ TO SPECIFY the parsed data with the missing values imputed using the BioFAModel 
#'    }
#'    \item{\code{expectations}:}{ expected values of the different nodes of the model. }
#'    
#'    \item{\code{training_stats}:}{ list with training statistics such as evidence lower bound (ELBO) and number of factors per iteration}
#'    \item{\code{data_options}:}{ list with the data processing options such as whether to center or scale the data.}
#'    \item{\code{training_options}:}{ list with the training options such as maximum number of iterations, tolerance for convergence, etc.}
#'    \item{\code{model_options}:}{ list with the model options such as likelihoods, whether the intercept was learnt, etc.}
#'    \item{\code{dimensions}:}{ list with the relevant dimensionalities of the model. 
#'    N for the number of samples, 
#'    P for the number of sample groups, 
#'    D for the number of features of each view,
#'    M for the number of feature_grups (views), 
#'    K for the number of factors.}
#'    \item{\code{status}:}{Auxiliary variable indicating whether the model has been trained.}
#'}
#' @name BioFAModel
#' @rdname BioFAModel
#' @aliases BioFAModel-class
#' @exportClass BioFAModel
setClass("BioFAModel", 
         slots=c(input_data       = "list",
                 training_data    = "list",
                 imputed_data     = "list",
                 expectations     = "list", 
                 training_stats   = "list",
                 training_options = "list",
                 data_options     = "list",
                 model_options    = "list",
                 dimensions       = "list",
                 on_disk          = "logical",
                 cache            = "list",
                 status           = "character")
)

# Printing method
setMethod("show", "BioFAModel", function(object) {
  
  if (!.hasSlot(object, "dimensions") | length(object@dimensions) == 0)
    stop("Error: dimensions not defined")
  if (!.hasSlot(object, "status") | length(object@status) == 0)
    stop("Error: status not defined")
  
  if (object@status == "trained") {
    nfactors <- object@dimensions[["K"]]
    cat(sprintf("Trained BioFAModel with the following characteristics: \n Number of views: %d \n View names: %s \n Number of features per view: %s \n Number of sample groups: %d \n Sample groups names: %s \n Number of samples per group: %s \n Number of factors: %d \n",
                object@dimensions[["M"]], paste(views_names(object),  collapse=" "), paste(as.character(object@dimensions[["D"]]), collapse=" "), 
                object@dimensions[["P"]], paste(groups_names(object), collapse=" "), paste(as.character(object@dimensions[["N"]]), collapse=" "), 
                nfactors))
  } else {
    cat(sprintf("Untrained BioFAModel model with the following characteristics: \n Number of views: %d \n View names: %s \n Number of features per view: %s \n Number of sample groups: %d \n Sample groups names: %s \n Number of samples per group: %s \n",
                object@dimensions[["M"]], paste(views_names(object),  collapse=" "), paste(as.character(object@dimensions[["D"]]), collapse=" "),
                object@dimensions[["P"]], paste(groups_names(object), collapse=" "), paste(as.character(object@dimensions[["N"]]), collapse=" ")))
  }
  cat("\n")
})


