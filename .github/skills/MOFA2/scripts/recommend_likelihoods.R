#!/usr/bin/env Rscript

# Propose a MOFA2 likelihood for each view by measuring its distribution.
#
# Why this is a script and not a prose rule: the choice is not "binary ->
# bernoulli, counts -> poisson". It depends on width and sparsity together. A
# narrow, dense binary view (somatic mutations) is modelled well by bernoulli,
# which is what the CLL tutorial does. A very wide binary view (scATAC peaks)
# makes bernoulli slow at that width, and the 10x Multiome tutorial instead
# TF-IDF-normalises peaks to a continuous scale and models them as gaussian.
# Those two both "look binary" to a naive check. Measuring the data and showing
# the numbers is more reliable than any sentence.
#
# Scope: the likelihood is a model option, so it is proposed here. How a view
# should be transformed is not — that depends on the assay and on processing
# upstream of this matrix. Where supplying a view on a different scale would
# change the answer, the output says so and leaves that choice to the analyst.
#
# Every proposal is printed with the measurements behind it. Treat the output as
# something to sanity-check, not a decision.
#
# Usage:
#   Rscript recommend_likelihoods.R --data <file>
#
#   --data   .rds holding a named list of matrices (features x samples),
#            a MOFA object, or an .hdf5 MOFA model

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
  self  <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  lines <- readLines(self, warn = FALSE)[-1]
  start <- which(grepl("^#", lines))[1]
  end   <- start + which(!grepl("^#", lines[start:length(lines)]))[1] - 2
  cat(sub("^# ?", "", lines[start:end]), sep = "\n")
  quit(status = 0)
}

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1]
}
data_path <- get_arg("--data")
if (is.null(data_path)) stop("--data is required. Use --help for usage.", call. = FALSE)

# Thresholds, kept explicit so they can be audited and adjusted. They are
# deliberately soft — the script reports which side of each line a view falls on
# rather than treating any one as decisive.
NARROW_MAX_WIDTH <- 1000   # binary views at or below this are bernoulli candidates
WIDE_WIDTH       <- 5000   # views at or above this are "wide" (feature-heavy)
SPARSE_FRAC      <- 0.90   # fraction of zeros above which a view is "very sparse"

# ---- load views (shared shape with profile_mofa_data.R) --------------------

load_views <- function(path) {
  if (grepl("\\.hdf5$|\\.h5$", path, ignore.case = TRUE)) {
    if (!requireNamespace("MOFA2", quietly = TRUE)) stop("MOFA2 needed to read .hdf5 models")
    obj <- MOFA2::load_model(path, load_data = TRUE)
    return(lapply(MOFA2::get_data(obj), function(v) do.call(cbind, v)))
  }
  x <- readRDS(path)
  if (inherits(x, "MOFA")) return(lapply(MOFA2::get_data(x), function(v) do.call(cbind, v)))
  if (is.list(x) && all(vapply(x, function(v) is.matrix(v) || inherits(v, "Matrix"), logical(1)))) {
    if (is.null(names(x))) names(x) <- paste0("view", seq_along(x))
    return(x)
  }
  stop("--data must be a list of matrices, a MOFA object, or an .hdf5 model", call. = FALSE)
}

views <- load_views(data_path)

# ---- measure one view ------------------------------------------------------

# Sample nonzero/all values cheaply from dense or sparse matrices.
view_values <- function(v) {
  if (inherits(v, "Matrix")) {
    nz <- v@x                                  # stored (nonzero) entries
    n_total <- prod(dim(v))
    n_zero  <- n_total - length(nz)
    list(nonzero = nz[is.finite(nz)], sparsity = n_zero / n_total, sparse = TRUE)
  } else {
    vals <- as.vector(v); vals <- vals[is.finite(vals)]
    list(nonzero = vals[vals != 0], sparsity = mean(vals == 0, na.rm = TRUE),
         all = vals, sparse = FALSE)
  }
}

profile_view <- function(v, name) {
  vv <- view_values(v)
  vals <- if (vv$sparse) vv$nonzero else vv$all      # for range/integer checks
  # unique-value count: capped sample for speed on huge views
  samp <- if (length(vals) > 1e5) vals[sample(length(vals), 1e5)] else vals
  uniq <- unique(samp)

  is_binary   <- all(uniq %in% c(0, 1))
  nonneg      <- all(vals >= 0)
  all_integer <- all(vals == round(vals))
  has_neg     <- any(vals < 0)
  # mean-variance trend on non-negative data => count-like, wants stabilisation
  mv_rho <- NA_real_
  if (nonneg && !is_binary && nrow(v) >= 20) {
    keep <- if (nrow(v) > 2000) sample(nrow(v), 2000) else seq_len(nrow(v))
    sub  <- as.matrix(v[keep, , drop = FALSE])
    fm <- rowMeans(sub, na.rm = TRUE); fv <- apply(sub, 1, var, na.rm = TRUE)
    ok <- is.finite(fm) & is.finite(fv) & fm > 0
    if (sum(ok) > 10) mv_rho <- suppressWarnings(cor(fm[ok], fv[ok], method = "spearman"))
  }

  list(name = name, width = nrow(v), n = ncol(v), sparsity = vv$sparsity,
       n_unique = length(uniq), is_binary = is_binary, nonneg = nonneg,
       all_integer = all_integer, has_neg = has_neg,
       range = range(vals), mv_rho = mv_rho)
}

# ---- turn measurements into a recommendation -------------------------------

recommend <- function(p) {
  wide   <- p$width >= WIDE_WIDTH
  sparse <- p$sparsity >= SPARSE_FRAC

  if (p$is_binary) {
    # Sparsity does NOT argue against bernoulli: somatic mutations are rare
    # (the CLL Mutations view is ~93% zeros) and bernoulli is exactly right for
    # them. What argues against it is WIDTH — tens of thousands of binary
    # features (scATAC peaks) make bernoulli slow, and the 10x Multiome tutorial
    # supplies such views on a continuous (TF-IDF) scale instead.
    if (p$width <= NARROW_MAX_WIDTH) {
      return(list(lik = "bernoulli",
        why = sprintf("binary, narrow (%d features%s) — bernoulli models it directly, as the CLL tutorial does for somatic mutations. Sparsity is fine here; rare 0/1 events are what bernoulli is for",
                      p$width, if (sparse) sprintf(", %.0f%% zeros", 100 * p$sparsity) else "")))
    }
    return(list(lik = "bernoulli / gaussian",
      why = sprintf("binary and wide (%d features%s) — bernoulli is valid but slow at this width. Gaussian applies to a continuous version of the view; the 10x Multiome tutorial TF-IDF-normalises scATAC peaks for that reason",
                    p$width, if (sparse) sprintf(", %.0f%% zeros", 100 * p$sparsity) else "")))
  }

  if (p$nonneg && p$all_integer && p$n_unique > 2) {          # count data
    trend <- if (!is.na(p$mv_rho) && p$mv_rho > 0.5)
      sprintf(" Mean-variance trend present (rho=%.2f), against gaussian's homoscedasticity assumption.", p$mv_rho) else ""
    # poisson is documented for small counts (getting_started_R.Rmd); flag the
    # range rather than ruling either likelihood in or out.
    small_counts <- p$range[2] <= 30
    pois <- if (small_counts) " Counts are small, the range poisson is documented for."
            else sprintf(" Poisson is documented for small counts; this view reaches %g.", p$range[2])
    scale_note <- " Tutorial practice for sequencing counts is library size correction plus variance stabilisation, then gaussian."
    if (wide && sparse) {
      return(list(lik = if (small_counts) "poisson / gaussian" else "gaussian (check input)",
        why = sprintf("wide (%d features), sparse (%.0f%% zeros) counts — single-cell-like.%s%s%s",
                      p$width, 100 * p$sparsity, pois, scale_note, trend)))
    }
    return(list(lik = if (small_counts) "poisson / gaussian" else "gaussian (check input)",
      why = sprintf("count data (integers, range %g-%g).%s%s%s",
                    p$range[1], p$range[2], pois, scale_note, trend)))
  }

  if (p$has_neg) {
    return(list(lik = "gaussian",
      why = sprintf("continuous with negative values (range %.1f to %.1f) — already transformed (log / z-score / M-value / CLR). Gaussian is the natural fit",
                    p$range[1], p$range[2])))
  }

  if (p$nonneg && p$range[2] <= 1 && p$n_unique > 20) {
    return(list(lik = "gaussian (check input)",
      why = "continuous but bounded in [0,1] — proportions, or methylation beta/rate values. Gaussian has unbounded support; the scMethylation cortex tutorial models M-values rather than rates. Confirm this is the intended input scale"))
  }

  list(lik = "gaussian",
    why = sprintf("continuous (range %.2f to %.2f), many distinct values — gaussian", p$range[1], p$range[2]))
}

# ---- report ----------------------------------------------------------------

cat("\n=== LIKELIHOOD PROPOSALS ===\n")
cat("Measurements first, proposal second. The likelihood is a model option; the\n")
cat("scale a view is supplied on is your decision, so some views list both.\n\n")

profiles <- Map(profile_view, views, names(views))
recs <- lapply(profiles, recommend)

tab <- data.frame(
  view        = vapply(profiles, `[[`, "", "name"),
  features    = vapply(profiles, `[[`, 0, "width"),
  samples     = vapply(profiles, `[[`, 0, "n"),
  distinct    = vapply(profiles, `[[`, 0, "n_unique"),
  zeros       = sprintf("%.0f%%", 100 * vapply(profiles, `[[`, 0, "sparsity")),
  kind        = vapply(profiles, function(p)
                  if (p$is_binary) "binary"
                  else if (p$nonneg && p$all_integer && p$n_unique > 2) "counts"
                  else if (p$has_neg) "continuous(+/-)"
                  else "continuous(+)", ""),
  recommend   = vapply(recs, `[[`, "", "lik"),
  row.names = NULL, stringsAsFactors = FALSE)
print(tab, right = FALSE)

cat("\nReasoning:\n")
for (i in seq_along(profiles)) {
  cat(sprintf("\n  [%s] -> %s\n", profiles[[i]]$name, recs[[i]]$lik))
  cat("     ", recs[[i]]$why, "\n", sep = "")
}

# Copy-pasteable code: gaussian baseline, then only the views that differ.
cat("\n=== SUGGESTED CODE ===\n")
cat("# Auto-inference of likelihoods is currently DISABLED in the package\n")
cat("# (the .infer_likelihoods() call is commented out in prepare_mofa.R), so the\n")
cat("# default is gaussian for EVERY view. Set non-gaussian views explicitly:\n\n")
cat("model_opts <- get_default_model_options(MOFAobject)\n")
for (i in seq_along(profiles)) {
  lik <- recs[[i]]$lik; nm <- profiles[[i]]$name
  if (identical(lik, "bernoulli")) {
    cat(sprintf('model_opts$likelihoods["%s"] <- "bernoulli"\n', nm))
  } else if (grepl("^bernoulli / ", lik)) {
    cat(sprintf('# model_opts$likelihoods["%s"] <- "bernoulli"   # or leave gaussian if the view is supplied continuous — see reasoning\n', nm))
  } else if (grepl("^poisson / ", lik)) {
    cat(sprintf('# model_opts$likelihoods["%s"] <- "poisson"     # or leave gaussian if the view is supplied continuous — see reasoning\n', nm))
  } else if (grepl("check input", lik)) {
    cat(sprintf('# "%s": gaussian by default — confirm the input scale is the intended one (see reasoning)\n', nm))
  }
}
cat("\n")