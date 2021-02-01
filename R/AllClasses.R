
##########################################################
## Define a general class to store a MOFA trained model ##
##########################################################

#' @title Class to store a mofa model
#' @description
#' The \code{MOFA} is an S4 class used to store all relevant data to analyse a MOFA model
#' @slot data The input data
#' @slot intercepts Feature intercepts
#' @slot samples_metadata Samples metadata
#' @slot features_metadata Features metadata.
#' @slot imputed_data The imputed data.
#' @slot expectations expected values of the factors and the loadings.
#' @slot dim_red non-linear dimensionality reduction manifolds.
#' @slot training_stats model training statistics.
#' @slot data_options Data processing options.
#' @slot training_options Model training options.
#' @slot stochastic_options Stochastic variational inference options.
#' @slot model_options Model options.
#' @slot mefisto_options  Options for the use of MEFISO
#' @slot dimensions Dimensionalities of the model: 
#'    M for the number of views, 
#'    G for the number of groups,
#'    N for the number of samples (per group),
#'    C for the number of covariates per sample,
#'    D for the number of features (per view),
#'    K for the number of factors.
#' @slot on_disk Logical indicating whether data is loaded from disk.
#' @slot cache Cache.
#' @slot status Auxiliary variable indicating whether the model has been trained.
#' @slot covariates optional slot to store sample covariate for training in MEFISTO
#' @slot covariates_warped optional slot to store warped sample covariate for training in MEFISTO
#' @slot interpolated_Z optional slot to store interpolated factor values (used only with MEFISTO)
#' @name MOFA
#' @rdname MOFA
#' @aliases MOFA-class
#' @exportClass MOFA

setClassUnion("listOrNULL",members = c("list","NULL"))
setClass("MOFA", 
        slots=c(
            data                = "list",
            covariates          = "listOrNULL",
            covariates_warped   = "listOrNULL",
            intercepts          = "list",
            imputed_data        = "list",
            interpolated_Z      = "list",
            samples_metadata    = "list",
            features_metadata   = "list",
            expectations        = "list", 
            training_stats      = "list",
            data_options        = "list",
            model_options       = "list",
            training_options    = "list",
            stochastic_options  = "list",
            mefisto_options      = "list",
            dimensions          = "list",
            on_disk             = "logical",
            dim_red             = "list",
            cache               = "list",
            status              = "character"
        )
)

# Printing method
setMethod("show", "MOFA", function(object) {
  
  if (!.hasSlot(object, "dimensions") || length(object@dimensions) == 0)
    stop("Error: dimensions not defined")
  if (!.hasSlot(object, "status") || length(object@status) == 0)
    stop("Error: status not defined")
  
  if (object@status == "trained") {
    nfactors <- object@dimensions[["K"]]
    if(!.hasSlot(object, "covariates") || is.null(object@covariates)) {
      cat(sprintf("Trained MOFA with the following characteristics: \n Number of views: %d \n Views names: %s \n Number of features (per view): %s \n Number of groups: %d \n Groups names: %s \n Number of samples (per group): %s \n Number of factors: %d \n",
                  object@dimensions[["M"]], paste(views_names(object),  collapse=" "), paste(as.character(object@dimensions[["D"]]), collapse=" "),
                  object@dimensions[["G"]], paste(groups_names(object), collapse=" "), paste(as.character(object@dimensions[["N"]]), collapse=" "),
                  nfactors))
    } else {
      cat(sprintf("Trained MEFISTO with the following characteristics: \n Number of views: %d \n Views names: %s \n Number of features (per view): %s \n Number of groups: %d \n Groups names: %s \n Number of samples (per group): %s \n Number of covariates per sample: %d \n Number of factors: %d \n",
                  object@dimensions[["M"]], paste(views_names(object),  collapse=" "), paste(as.character(object@dimensions[["D"]]), collapse=" "),
                  object@dimensions[["G"]], paste(groups_names(object), collapse=" "), paste(as.character(object@dimensions[["N"]]), collapse=" "),
                  object@dimensions[["C"]], nfactors))
    }
  } else {
    if(!.hasSlot(object, "covariates") || is.null(object@covariates)) {
      cat(sprintf("Untrained MOFA model with the following characteristics: \n Number of views: %d \n Views names: %s \n Number of features (per view): %s \n Number of groups: %d \n Groups names: %s \n Number of samples (per group): %s \n ",
                  object@dimensions[["M"]], paste(views_names(object),  collapse=" "), paste(as.character(object@dimensions[["D"]]), collapse=" "),
                  object@dimensions[["G"]], paste(groups_names(object), collapse=" "), paste(as.character(object@dimensions[["N"]]), collapse=" ")))
    } else {
      cat(sprintf("Untrained MEFISTO model with the following characteristics: \n Number of views: %d \n Views names: %s \n Number of features (per view): %s \n Number of groups: %d \n Groups names: %s \n Number of samples (per group): %s \n Number of covariates per sample: %d \n ",
                  object@dimensions[["M"]], paste(views_names(object),  collapse=" "), paste(as.character(object@dimensions[["D"]]), collapse=" "),
                  object@dimensions[["G"]], paste(groups_names(object), collapse=" "), paste(as.character(object@dimensions[["N"]]), collapse=" "),
                  object@dimensions[["C"]]))
    }
  }
  cat("\n")
})


