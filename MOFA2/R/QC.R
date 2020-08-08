#' @title Quality control
#' @name quality_control
#' @description Function to do quality control on a \code{\link{MOFA}} object.
#' @param object a trained \code{\link{MOFA}} object.
#' @param verbose logical indicating whether to generate a verbose output.
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("exdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Do quality control
#' model <- quality_control(model, verbose = TRUE)
quality_control <- function(object, verbose = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Check views names
  if (verbose == TRUE) message("Checking views names...")
  stopifnot(!is.null(views_names(object)))
  stopifnot(!duplicated(views_names(object)))
  if (any(grepl("/", views_names(object)))) {
    stop("Some of the views names contain `/` symbol, which is not supported.
  This can be fixed e.g. with:
    views_names(object) <- gsub(\"/\", \"-\", views_names(object))")
  }
  
  # Check groups names
  if (verbose == TRUE) message("Checking groups names...")
  if (any(grepl("/", groups_names(object)))) {
    stop("Some of the groups names contain `/` symbol, which is not supported.
  This can be fixed e.g. with:
    groups_names(object) <- gsub(\"/\", \"-\", groups_names(object))")
  } 
  stopifnot(!is.null(groups_names(object)))
  stopifnot(!duplicated(groups_names(object)))
  
  # Check samples names
  if (verbose == TRUE) message("Checking samples names...")
  stopifnot(!is.null(samples_names(object)))
  stopifnot(!duplicated(unlist(samples_names(object))))
  
  # Check dimensionalities in the input data
  if (verbose == TRUE) message("Checking dimensions...")
  N <- object@dimensions$N
  D <- object@dimensions$D
  for (i in views_names(object)) {
    for (j in groups_names(object)) {
      stopifnot(ncol(object@data[[i]][[j]]) == N[[j]])
      stopifnot(nrow(object@data[[i]][[j]]) == D[[i]])
      stopifnot(length(colnames(object@data[[i]][[j]])) == N[[j]])
      stopifnot(length(rownames(object@data[[i]][[j]])) == D[[i]])
    }
  }
  
  # Check that there are no features with complete missing values (across all views)
  if (isTRUE(object@data_options[["loaded"]])) {
    if (verbose == TRUE) message("Checking there are no features with complete missing values...")
    for (i in views_names(object)) {
      if (!(is(object@data[[i]][[1]], "dgCMatrix") || is(object@data[[i]][[1]], "dgTMatrix"))) {
        tmp <- as.data.frame(sapply(object@data[[i]], function(x) rowMeans(is.na(x)), simplify = TRUE))
        if (any(unlist(apply(tmp, 1, function(x) mean(x==1)))==1))
          warning("You have features which do not contain a single observation in any group, consider removing them...")
      }
    }
  }
    
  # Sanity checks that are exclusive for an untrained model  
  if (object@status == "untrained") {
    
    # Check features names
    if (verbose == TRUE) message("Checking features names...")
    tmp <- lapply(object@data, function(x) unique(lapply(x,rownames)))
    for (x in tmp) stopifnot(length(x)==1)
    for (x in tmp) if (any(duplicated(x[[1]]))) stop("There are duplicated features names within the same view. Please rename")
    all_names <- unname(unlist(tmp))
    duplicated_names <- unique(all_names[duplicated(all_names)])
    if (length(duplicated_names)>0) 
      warning("There are duplicated features names across different views. We will add the suffix *_view* only for those features 
            Example: if you have both TP53 in mRNA and mutation data it will be renamed to TP53_mRNA, TP53_mutation")
    for (i in names(object@data)) {
      for (j in names(object@data[[i]])) {
        tmp <- which(rownames(object@data[[i]][[j]]) %in% duplicated_names)
        if (length(tmp)>0) {
          rownames(object@data[[i]][[j]])[tmp] <- paste(rownames(object@data[[i]][[j]])[tmp], i, sep="_")
        }
      }
    }
    
  # Sanity checks that are exclusive for a trained model  
  } else if (object@status == "trained") {
    
    # Check expectations
    stopifnot(all(c("W", "Z") %in% names(object@expectations)))
    if (verbose == TRUE) message("Checking expectations...")
    stopifnot(all(sapply(object@expectations$W, is.matrix)))
    stopifnot(all(sapply(object@expectations$Z, is.matrix)))
    
    # Check for intercept factors
    if (isTRUE(object@data_options[["loaded"]])) { 
      if (verbose == TRUE) message("Checking for intercept factors...")
      if (!is.null(object@data)) {
        factors <- do.call("rbind",get_factors(object))
        r <- suppressWarnings( t(do.call('rbind', lapply(object@data, function(x) 
          abs(cor(colMeans(do.call("cbind",x),na.rm=TRUE),factors, use="pairwise.complete.obs"))
        ))) )
        intercept_factors <- which(rowSums(r>0.75)>0)
        if (length(intercept_factors)) {
            warning(sprintf("Factor(s) %s are strongly correlated with the total number of expressed features for at least one of your omics. Such factors appear when there are differences in the total 'levels' between your samples, *sometimes* because of poor normalisation in the preprocessing steps.\n",paste(intercept_factors,collapse=", ")))
        }
      }
    }
  
    # Check for correlated factors
    if (verbose == TRUE) message("Checking for highly correlated factors...")
    Z <- do.call("rbind",get_factors(object))
    op <- options(warn=-1) # suppress warnings
    
    noise <- matrix(rnorm(n=length(Z), mean=0, sd=1e-10), nrow(Z), ncol(Z))
    tmp <- cor(Z+noise); diag(tmp) <- NA
    options(op) # activate warnings again
    if (max(tmp,na.rm=TRUE)>0.5) {
      warning("The model contains highly correlated factors (see `plot_factor_cor(MOFAobject)`). \nWe recommend that you train the model with less factors and that you let it train for a longer time.\n")
    }

  }
  
  return(object)  
}
