#' @title One-liner wrapper to create a MOFA2 object ready for training
#' @name mofa2_from_mae
#' @description \code{mofa2_from_mae} generates a MOFA2 model from a MultiAssayExperiment with default options that is ready for training.
#' This is a one-line wrapper that combines the use of \code{\link[create_mofa_from_MultiAssayExperiment]{create_mofa_from_MultiAssayExperiment}} followed by \code{\link[prepare_mofa]{prepare_mofa}}.
#' @param mae a \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}}
#' @param assay.names List of assays to include in MOFA model (required).
#' @param groups One of the variable names of coldata from which groups should be derived. Default is NULL, only relevant when using the multi-group framework.
#' @param extract_metadata Boolean specifying whether sample metadata should be propagated from the MultiAssayExperiment to the MOFA2 object  (default is FALSE).
#' @param ... additional parameters for MOFA, named as the options of \code{\link[prepare_mofa]{prepare_mofa}}.
#' 
#' @return 
#' Returns an untrained \code{\link[MOFA2:MOFA]{MOFA}} with specified options
#' filled in the corresponding slots
#'
#' 
#' @examples
#' # Load package and import dataset
#' library(mia)
#' data("HintikkaXOData", package = "mia")
#' mae <- HintikkaXOData
#' 
#' # Prepare basic model with selected assays
#' prep_model <- mofa2_from_mae(mae, assay.names = c("counts", "nmr", "signals"))
#' 
#' # Specify grouping variable and extract metadata
#' prep_model <- mofa2_from_mae(mae, assay.names = c("counts", "nmr", "signals"),
#'                     groups = "Diet", extract_metadata = TRUE)
#'
#' # Modify MOFA options with corresponding arguments
#' prep_model <- mofa2_from_mae(mae, assay.names = c("counts", "nmr", "signals"),
#'                     num_factors = 5)
mofa2_from_mae <- function(mae, assay.names, groups = NULL, extract_metadata = FALSE, ...) {
  
  # Select assays of each experiment for MOFA
  mae <- .select_assays(mae, assay.names)
  
  # Create MOFA from selected experiments
  mofa_obj <- create_mofa_from_MultiAssayExperiment(
    mae,
    groups = groups,
    extract_metadata = extract_metadata
  )
  
  # Make a list of arguments for MOFA
  mofa_args <- list(
    object = mofa_obj,
    data_options = .set_opts(get_default_data_options(mofa_obj), ...),
    model_options = .set_opts(get_default_model_options(mofa_obj), ...),
    training_options = .set_opts(get_default_training_options(mofa_obj), ...),
    mefisto_options = .set_opts(get_default_mefisto_options(mofa_obj), ...)
  )
  
  # Add stochastic options if stochastic is turned on
  if ( mofa_args[["training_options"]][["stochastic"]] ){
    mofa_args[["stochastic_options"]] <- .set_opts(get_default_stochastic_options(mofa_obj), ...)
  }
  
  # Prepare MOFA
  mofa_obj <- do.call("prepare_mofa", mofa_args)
  
  return(mofa_obj)
}
  
########################## HELP FUNCTIONS ##########################

# Select assays to be included in MOFA
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SummarizedExperiment assay assays assayNames
.select_assays <- function(mae, assay.names) {
  # Give corresponding experiment names to assay.names
  names(assay.names) <- names(experiments(mae))
  # For every experiment in MAE
  for ( exp in names(experiments(mae)) ){
    # Keep only selected assay.type from a given experiment
    assays(mae[[exp]]) <- list(assay(mae[[exp]], assay.names[[exp]]))
    # Update assay names
    assayNames(mae[[exp]]) <- assay.names[[exp]]
  }
  return(mae)
}

# Combine custom options found in ... with default options
.set_opts <- function(default, ...) {
  # For every option in a set (data, model, train, ...)
  for ( opt in names(default) ){
    # If that option is found among arguments
    if ( opt %in% names(list(...)) ){
      # Replace default with value specified in arguments
      default[[opt]] <- list(...)[[opt]]
    }
  }
  return(default)
}
