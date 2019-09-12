#' @title Quality control
#' @name quality_control
#' @description Function to do quality control on a \code{\link{MOFA}} object. \cr
#' @param object a trained \code{\link{MOFA}} object.
#' @param verbose logical indicating whether to generate a verbose output.
#' @export
#'
quality_control <- function(object, verbose = FALSE) {
  
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Check views names
  if (verbose == TRUE) message("Checking views names...")
  stopifnot(!is.null(views(object)))
  stopifnot(!duplicated(views(object)))
  if (any(grepl("/", views(object)))) {
    stop("Some of the views names contain `/` symbol, which is not supported.
  This can be fixed e.g. with:
    views(object) <- gsub(\"/\", \"-\", views(object))")
  }
  
  # Check groups names
  if (verbose == TRUE) message("Checking groups names...")
  if (any(grepl("/", groups(object)))) {
    stop("Some of the groups names contain `/` symbol, which is not supported.
  This can be fixed e.g. with:
    groups(object) <- gsub(\"/\", \"-\", groups(object))")
  } 
  stopifnot(!is.null(groups(object)))
  stopifnot(!duplicated(groups(object)))
  
  # Check samples names
  if (verbose == TRUE) message("Checking samples names...")
  stopifnot(!is.null(samples(object)))
  stopifnot(!duplicated(unlist(samples(object))))
  
  # Check dimensionalities in the input data
  N <- object@dimensions$N
  D <- object@dimensions$D
  for (i in views(object)) {
    for (j in groups(object)) {
      stopifnot(ncol(object@data[[i]][[j]]) == N[[j]])
      stopifnot(nrow(object@data[[i]][[j]]) == D[[i]])
      stopifnot(length(colnames(object@data[[i]][[j]])) == N[[j]])
      stopifnot(length(rownames(object@data[[i]][[j]])) == D[[i]])
    }
  }
  
  # Check that there are no features with complete missing values (across all views)
  if (verbose == TRUE) message("Checking there are no features with complete missing values...")
  for (i in views(object)) {
    tmp <- as.data.frame(sapply(object@data[[i]], function(x) rowMeans(is.na(x)), simplify = T))
    if (any(unlist(apply(tmp, 1, function(x) mean(x==1)))==1))
      warning("You have features which do not contain a single observation in any group, consider removing them...")
  }
    
  # Check that the likelihoods match the data distribution
  # if (verbose == TRUE) message("Checking likelihooods...")
  # predicted_lik <- .inferLikelihoods(object)
  # for (view in viewNames(object)) {
  #   lk <- object@ModelOptions$likelihood[view]
  #   if (lk != predicted_lik[view])
  #     message(sprintf("Warning, view %s should follow a %s distribution rather than %s ", view, predicted_lik[view], lk))
  # }
  
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
  }
  
  return(object)  
}
