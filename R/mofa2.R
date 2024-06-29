#' @title One-liner wrapper to create a MOFA2 object from multiple input objects and ready for training
#' @name mofa2
#' @description
#' This is a one-line wrapper that combines the use of \code{\link[create_mofa]{create_mofa}} followed by \code{\link[prepare_mofa]{prepare_mofa}}.
#' Please read the documentation of the corresponding functions for more details.
#' @param data input data. See \code{\link[create_mofa]{create_mofa}} for details on input data formats.
#' @param assays Assays to include in MOFA model. Required for \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}}, \code{\link[SingleCellExperiment:SingleCellExperiment]{SingleCellExperiment}} and \code{\link[SeuratObject:CreateSeuratObject]{Seurat}} objects.
#' @param groups One of the variable names of sample metadata from which groups should be derived. Only relevant when using the multi-group framework with one of the three formats above. Default is NULL.
#' @param extract_metadata Boolean specifying whether sample metadata from the input object should be propagated to the MOFA2 object  (default is TRUE).
#' @param ... additional parameters for MOFA, named as the options of \code{\link[prepare_mofa]{prepare_mofa}}.
#' 
#' @return 
#' Returns an untrained \code{\link[MOFA2:MOFA]{MOFA}} with specified options
#' filled in the corresponding slots
#'
#' 
#' @examples

#' # Example with MultiAssayExperiment
#' library(mia)
#' data("HintikkaXOData", package = "mia")
#' mae <- HintikkaXOData
#' mofa_obj <- mofa2(mae, assays = c("counts", "nmr"))
#' mofa_obj <- mofa2(mae, assays = c("counts", "nmr"), num_factors = 5)

#' # Example with long data.frame format (two views and two groups)
#' file <- system.file("extdata", "test_data.RData", package = "MOFA2")
#' load(file) 
#' mofa_obj <- mofa2(dt)

#' # Example with list of matrices
#' mtx <- make_example_data()$data
#' mofa_obj <- mofa2(mtx)
mofa2 <- function(data, assays = NULL, groups = NULL, extract_metadata = TRUE, ...) {
  
  # Select assays of each experiment for MOFA
  # TO-DO: IF INPUT DATA IS MAE, SINGLECELLEXPERIMENT OR SEURAT, MAKE USE ASSAYS IS NOT NUL
  # mofa_obj <- create_mofa(data, assays, groups, extract_metadata, ...)
  mofa_obj <- create_mofa(data, assays, groups, extract_metadata)
  
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

###########################
###### HELP FUNCTIONS #####
###########################

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
