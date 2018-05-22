
infer_likelihoods <- function(object) {
  likelihood <- rep(x="gaussian", times=object@Dimensions$M)
  names(likelihood) <- views_names(object)
  
  for (view in views_names(object)) {
    data <- get_training_data(object, view)[[1]]
    # if (all(data %in% c(0,1,NA))) {
    if (length(unique(data[!is.na(data)]))==2) {
      likelihood[view] <- "bernoulli"
    } else if (all(data[!is.na(data)]%%1==0)) {
      likelihood[view] <- "poisson"
    }
  }
  
  return(likelihood)
}

.update_old_model <- function(object) {
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")  
  
  # Update node names
  if ("SW" %in% names(object@expectations)) {
    # object@model_options$schedule[object@model_options$schedule == "SW"] <- "W" # schedule is depreciated from model_options
    names(object@expectations)[names(object@expectations) == "SW"] <- "W"
    colnames(object@training_stats$elbo_terms)[colnames(object@training_stats$elbo_terms)=="SW"] <- "W"
  }
  if ("SZ" %in% names(object@expectations)) {
    names(object@expectations)[names(object@expectations) == "SZ"] <- "Z"
    colnames(object@training_stats$elbo_terms)[colnames(object@training_stats$elbo_terms)=="SZ"] <- "Z"
  }

  # Update expectations
  if (is.list(object@expectations$Z[[1]]) & ("E" %in% names(object@expectations$Z[[1]]))) {
    # Multi-view nodes
    for (m in views_names(object)) {
      for (node in object@model_options$nodes$multiview_node) {
        if (node %in% names(object@expectations)){
          object@expectations[[node]][[m]] <- object@expectations[[node]][[m]]$E
        }
      }
    }
    # Multi-group nodes
    for (p in groups_names(object)) {
      for (node in object@model_options$nodes$multigroup_nodes) {
        if (node %in% names(object@expectations)){
          object@expectations[[node]][[p]] <- object@expectations[[node]][[p]]$E
        }
      }
    }
    # Multi-view & multi-group nodes
    for (m in views_names(object)) {
      for (p in groups_names(object)) {
        for (node in object@model_options$nodes$twodim_nodes) {
          object@expectations[[node]][[m]][[p]] <- object@expectations[[node]][[m]][[p]]$E
        }
      }
    }
  }
  
  
  # update learn_mean to learn_intercept
  if ("learn_mean" %in% names(object@model_options)) {
    tmp <- names(object@model_options)
    tmp[tmp=="learn_mean"] <- "learn_intercept"
    names(object@model_options) <- tmp
  }
  object@model_options$learn_intercept <- as.logical(object@model_options$learn_intercept)

  # Set the status as trained if it wasn't set before
  if((!.hasSlot(object, "status")) | (length(object@status) == 0))
    object@status <- "trained"

  
  return(object)
}

# Set view names and group names for nested list objects (e.g. Y)
.name_views_and_groups <- function(nested_list, view_names, group_names) {
  names(nested_list) <- view_names
  for (view in view_names) { names(nested_list[[view]]) <- group_names }
  nested_list
}

# Function to find factors that act like an intercept term for the sample, 
# which means that they capture global mean effects
find_intercept_factors <- function(object, cor_threshold = 0.8) {
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")  
  
  data <- get_training_data(object)
  factors <- get_factors(object, include_intercept = F)
  
  r <- lapply(data, function(x) abs(cor(apply(x,2,mean),factors, use="complete.obs")))
  for (i in names(r)) {
    if (any(r[[i]]>cor_threshold))
      cat(paste0("Warning: factor ",which(r[[i]]>cor_threshold)," is capturing a size factor effect in ", i, " view, which indicates that input data might not be properly normalised...\n"))
  }
}


subset_augment <- function(mat, pats) {
  pats <- unique(pats)
  mat <- t(mat)
  aug_mat <- matrix(NA, ncol=ncol(mat), nrow=length(pats))
  aug_mat <- mat[match(pats,rownames(mat)),,drop=FALSE]
  rownames(aug_mat) <- pats
  colnames(aug_mat) <- colnames(mat)
  return(t(aug_mat))
}


detect_passengers <- function(object, views = "all", groups = "all", factors = "all", r2_threshold = 0.03) {
  
  # Sanity checks
  if (class(object) != "BioFAModel") stop("'object' has to be an instance of BioFAModel")
  
  # Define views
  if (paste0(views, sep="", collapse="") == "all") { 
    views <- views_names(object) 
  } else {
    stopifnot(all(views %in% views_names(object)))  
  }
  M <- length(views)

  # Define groups
  if (paste0(groups, sep="", collapse="") == "all") { 
    groups <- groups_names(object) 
  } else {
    stopifnot(all(groups %in% groups_names(object)))  
  }
  H <- length(groups)
  
  # Define factors
  factors <- as.character(factors)
  if (paste0(factors, collapse="") == "all") { 
    factors <- factors_names(object)
  } else {
    stopifnot(all(factors %in% factors_names(object)))  
  }
  
  # Collect relevant data
  Z <- get_factors(object)
  
  # Identify factors unique to a single view by calculating relative R2 per factor
  r2 <- calculate_variance_explained(object, views = views, groups = groups, factors = factors)$r2_per_factor
  browser()
  unique_factors <- unique(unlist(lapply(groups, function(p) names(which(rowSums(r2[[p]]>=r2_threshold)==1)) )))
  
  # Mask samples that are unique in the unique factors
  missing <- lapply(get_training_data(object, views, groups), function(views) {
    lapply(views, function(group) {
      samples_names(object)[apply(group, 2, function(x) all(is.na(x)))]
    })
  })

  missing <- .name_views_and_groups(missing, views_names(object), groups_names(object))
  for (fctr in unique_factors) {
    # view <- names(which(r2[fctr,]>=r2_threshold))
    for (p in groups) {
      view <- colnames(r2[[p]][,which(r2[[p]][fctr,]>=r2_threshold), drop=F])
      if (!is.null(view)) {
        missing_samples <- missing[[view]][[p]]
        if (length(missing_samples) > 0) {
          Z[[p]][missing_samples, fctr] <- NA
        }
      }
    }
  }
  
  # Replace the latent matrix
  object@expectations$Z <- Z
  
  return(object)
  
}


flip_factor <- function(model, factor){
  for(groupnm in names(model@expectations$Z)) {
    model@expectations$Z[[groupnm]][,factor] <- - model@expectations$Z[[groupnm]][,factor]
  }
  for(viewnm in names(model@expectations$W)) {
    model@expectations$W[[viewnm]][,factor] <- -model@expectations$W[[viewnm]][,factor]
  }
return(model)
}




.check_and_get_views <- function(object, views) {
  if (is.numeric(views)) {
    stopifnot(all(views <= object@dimensions$M))
    views_names(object)[views] 
  } else {
    if (paste0(views, sep = "", collapse = "") == "all") { 
      views_names(object)
    } else {
      stopifnot(all(views %in% views_names(object)))
      views
    }
  }
}


.check_and_get_groups <- function(object, groups) {
  if (is.numeric(groups)) {
    stopifnot(all(groups <= object@dimensions$M))
    groups_names(object)[groups] 
  } else {
    if (paste0(groups, sep = "", collapse = "") == "all") { 
      groups_names(object)
    } else {
      stopifnot(all(groups %in% groups_names(object)))
      groups
    }
  }
}



.rep_string <- function(times, string, collapse = "") {
  paste(replicate(times, string), collapse = collapse)
}

.pad_left_with <- function(len, string, with = "") {
  wlen <- nchar(with)
  len  <- max(len - wlen, 0)
  paste0(with, paste(replicate(len, " "), collapse = ""), string)
}

.pad_left <- function(len, string) {
  .pad_left_with(len, string, with = "")
}

# Center and paste
.cpaste <- function(vals, cwidth, collapse = "") {
  vals <- sapply(vals, function(e) {
    e <- toString(e)
    lendiff <- cwidth - nchar(e)
    if (lendiff > 1) {
      paste0(.rep_string(ceiling(lendiff / 2), " "),
             e,
             .rep_string(floor(lendiff / 2), " "))
    } else {
      e
    }
  })
  paste(vals, collapse = collapse)
}

# Fancy printing method
vis <- function(object) {

  stopifnot(class(object) == "BioFAModel")
  
  if (!.hasSlot(object, "dimensions") | length(object@dimensions) == 0)
    stop("Error: dimensions not defined")
  if (!.hasSlot(object, "status") | length(object@status) == 0)
    stop("Error: status not defined")

  vis_lines <- ""

  lpad <- max(sapply(views_names(object), function(v) nchar(v)))
  wlim <- max(sapply(groups_names(object), function(v) nchar(v)))
  igr_sp <- .rep_string(5, " ")
  s <- 8             # extra lpadding shift
  w <- max(8, wlim)  # width of one block (minus 2 walls)
  hat    <- paste0(" ", .rep_string(w, "_"), " ")
  walls  <- paste0("|", .rep_string(w, " "), "|")
  ground <- paste0("|", .rep_string(w, "_"), "|")

  cat("
         \U2588︎\U2588︎\U2588︎\U2588︎\U2588︎     \U2588︎\U2588︎   \U2588\U2588︎\U2588︎\U2588︎\U2588︎
biofam   \U2588︎\U2588︎\U2588︎\U2588︎\U2588︎  =  \U2588︎\U2588︎ x \U2588︎\U2588︎\U2588︎\U2588︎\U2588︎
         \U2588︎\U2588︎\U2588︎\U2588︎\U2588︎     \U2588︎\U2588︎   
      ")

  groups_line    <- .pad_left(lpad + s, .cpaste(groups_names(object), w+2, collapse = igr_sp))
  nsamples_line  <- .pad_left(lpad + s, .cpaste(get_dimensions(object)$N, w+2, collapse = igr_sp))
  vis_lines      <- c(vis_lines, groups_line, nsamples_line)  

  for (m in 1:length(views_names(object))) {
    toprect_line   <- .pad_left(lpad + s, paste(.rep_string(get_dimensions(object)$P, hat, collapse = igr_sp)))
    midrect_line   <- .pad_left(lpad + s, paste(.rep_string(get_dimensions(object)$P, walls, collapse = igr_sp)))
    dfeatures_line <- .pad_left_with(lpad + s, 
                                     paste(.rep_string(get_dimensions(object)$P, walls, collapse = igr_sp)), 
                                     with = paste(c(views_names(object)[m], .cpaste(get_dimensions(object)$D[m], s)), collapse = ""))
    botrect_line   <- .pad_left(lpad + s, paste(.rep_string(get_dimensions(object)$P, ground, collapse = igr_sp)))

    vis_lines      <- c(vis_lines, toprect_line, midrect_line, dfeatures_line, botrect_line)  
  }

  cat(paste(vis_lines, collapse = "\n"))
  
  cat("\n\n")
}