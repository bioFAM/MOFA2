#' Create MOFA from MAE
#'
#' \code{mofa2} produces a prepared MOFA2 model from a MAE.
#'
#' @param mae a \code{\link[MultiAssayExperiment:MultiAssayExperiment]{MultiAssayExperiment}}
#'
#' @param experiments Name or index of experiments selected from \code{mae}.
#' 
#' @param assay.types List of assays to include in MOFA model (default is "counts").
#' 
#' @param groups One of the variable names of coldata from which groups should be derived (default is NULL).
#' 
#' @param extract_metadata Boolean specifying whethet metadata should be extracted (default is FALSE).
#' 
#' @param ... additional parameters for MOFA, named as the options
#' of \code{\link[MOFA2:prepare_mofa]{prepare_mofa}}.
#' \itemize{
#' \item \code{altexps}: Specifies alternative experiments for each experiment.
#' Only applicable for experiments in \code{SingleCellExperiment} format. If
#' \code{NULL}, the option is disabled. For disabling option for only certain
#' experiments, use \code{NA} for that experiment.
#' }
#' 
#' @return 
#' Returns an untrained \code{\link[MOFA2:MOFA]{MOFA}} with specified options
#' filled in the corresponding slots
#'
#' @name create_mofa_from_mae
#' 
#' @examples
#' \dontrun{
#' # Load package and import dataset
#' library(mia)
#' data("HintikkaXOData", package = "mia")
#' mae <- HintikkaXOData
#' 
#' # Prepare basic model with selected assays
#' prep_model <- mofa2(
#'     mae,
#'     experiments = names(mae),
#'     assay.types = c("counts", "nmr", "signals")
#'     )
#' 
#' # Specify grouping variable and extract metadata
#' prep_model <- mofa2(
#'     mae,
#'     experiments = names(mae),
#'     assay.types = c("counts", "nmr", "signals"),
#'     groups = "Diet", extract_metadata = TRUE
#'     )
#'
#' # Modify MOFA options with corresponding arguments
#' prep_model <- mofa2(
#'     mae,
#'     experiments = names(mae),
#'     assay.types = c("counts", "nmr", "signals"),
#'     num_factors = 5, stochastic = TRUE
#'     )
#' 
#' # Agglomerate and transform microbiome counts data
#' mae[[1]] <- agglomerateByRanks(mae[[1]])
#' mae[[1]] <- transformAssay(
#'     mae[[1]], method = "clr", pseudocount = 1,
#'     altexp = altExpNames(mae[[1]])
#'     )
#' # Create MOFA model
#' prep_model <- mofa2(
#'     mae,
#'     experiments = c(1, 2),
#'     altexps = c("Genus", NA),
#'     assay.types = c("clr", "nmr")
#'     )
#' }
NULL

#' @rdname create_mofa_from_mae
#' @export
setGeneric("mofa2", signature = c("mae"),
           function(mae, ...) standardGeneric("mofa2")
)

#' @rdname create_mofa_from_mae
#' @export
setMethod("mofa2", signature = c(mae = "MultiAssayExperiment"),
          function(mae, experiments, assay.types,
                   groups = NULL, extract_metadata = FALSE, ...){
            # Select experiments from MAE
            mae <- .select_experiments(mae, experiments)
            # Select alternative experiments
            mae <- .select_altexps(mae, ...)
            # Select assays of each experiment for MOFA
            mae <- .select_assays(mae, assay.types)
            # Create MOFA from selected experiments
            mdl <- create_mofa_from_MultiAssayExperiment(
              mae,
              groups = groups,
              extract_metadata = extract_metadata
            )
            # Make a list of arguments for MOFA
            mofa_args <- list(
              object = mdl,
              data_options = .set_opts(get_default_data_options(mdl), ...),
              model_options = .set_opts(get_default_model_options(mdl), ...),
              training_options = .set_opts(get_default_training_options(mdl), ...),
              mefisto_options = .set_opts(get_default_mefisto_options(mdl), ...)
            )
            # Add stochastic options if stochastic is turned on
            if ( mofa_args[["training_options"]][["stochastic"]] ){
              mofa_args[["stochastic_options"]] <- .set_opts(get_default_stochastic_options(mdl), ...)
            }
            # Prepare MOFA
            prep_mdl <- do.call("prepare_mofa", mofa_args)
            return(prep_mdl)
          }
)

########################## HELP FUNCTIONS ##########################

# Select experiments from MAE
#' @importFrom MultiAssayExperiment experiments
.select_experiments <- function(mae, experiments) {
  # Check that the value is correct
  is_name <- is.character(experiments) && length(experiments) > 0 &&
    length(experiments) <= length(experiments(mae)) &&
    all(experiments %in% names(mae))
  is_index <- is.numeric(experiments) && all(experiments%%1==0) &&
    length(experiments) > 0 &&
    length(experiments) <= length(experiments(mae)) &&
    all(experiments>0 & experiments<=length(experiments(mae)))
  if( !(is_name || is_index) ){
    stop("'experiments' must specify names of index of experiments of 'mae'")
  }
  # Subset experiments
  mae <- mae[,, experiments]
  # Check that all objects are SE
  all_SE <- lapply(experiments(mae), function(x) is(x, "SummarizedExperiment")) |>
    unlist() |> all()
  if( !all_SE ){
    stop("All experiments must be SummarizedExperiment objects.")
  }
  return(mae)
}

# Select optionally alternative experiments
.select_altexps <- function(mae, altexps = NULL, ...) {
  # Check that value is correct
  is_name <- is.character(altexps) && length(altexps) > 0 &&
    length(altexps) <= length(experiments(mae))
  is_index <- is.numeric(altexps) && all(altexps%%1==0) &&
    length(altexps) > 0 &&
    length(altexps) <= length(experiments(mae))
  if( !( is.null(altexps) || is_name || is_index) ){
    stop("'altexps' must be NULL or specify aletnative experiments for each ",
         "experiment.")
  }
  # If specified, select altExps from experiments
  if( !is.null(altexps) ){
    if( !require("SingleCellExperiment") ){
      stop("To enable 'altexps' option, 'SingleCellExperiment' package must ",
           "be installed.")
    }
    names(altexps) <- names(mae)
    for ( exp in names(mae) ){
      # Get altExp if it is not NA, which disables altExp for single experiment
      if( !is.na(altexps[[exp]]) ){
        if( !is(mae[[exp]], "SingleCellExperiment") ){
          stop("Experiment '", exp, "' must be SingleCellExperiment object.")
        }
        # Check that alExp can be found
        is_name <- is.character(altexps[[exp]]) &&
          altexps[[exp]] %in% altExpNames(mae[[exp]])
        is_index <- is.numeric(altexps[[exp]]) && all(altexps[[exp]]%%1==0) &&
          altexps[[exp]]>0 && altexps[[exp]]<=length(alExps(mae[[exp]]) )
        if( !( is_name || is_index ) ){
          stop("'", altexps[[exp]], "' does not specify altExp from ",
               "experiment '", exp, "'.")
        }
        mae[[exp]] <- altExp(mae[[exp]], altexps[[exp]])
      }
    }
  }
  return(mae)
}

# Select assays to be included in MOFA
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SummarizedExperiment assay assays assayNames
.select_assays <- function(mae, assay.types) {
  # Check that value is correct
  is_name <- is.character(assay.types) && length(assay.types) > 0 &&
    length(assay.types) <= length(experiments(mae))
  if( !is_name ){
    stop("'assay.type' must specify name of assays. The lenght must equal to ",
         "'experiments'.")
  }
  # Give corresponding experiment names to assay.types
  names(assay.types) <- names(mae)
  # For every experiment in MAE
  for ( exp in names(mae) ){
    # Check that assay exists
    if( !assay.types[[exp]] %in% assayNames(mae[[exp]]) ){
      stop("Cannot find assay '", assay.types[[exp]], "' from experiment '",
           exp, "'.")
    }
    # Keep only selected assay.type from a given experiment
    assays(mae[[exp]]) <- assays(mae[[exp]])[ assay.types[[exp]] ]
  }
  return(mae)
}

# Combine custom options found in ... with default options
.set_opts <- function(default, ...) {
  user_options <- list(...)
  set_options <- intersect(names(default), names(user_options))
  default[set_options] <- user_options[set_options]
  return(default)
}
