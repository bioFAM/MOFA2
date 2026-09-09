#!/usr/bin/env Rscript

# Check whether each view looks ready for MOFA, or still needs preprocessing.
#
# MOFA assumes each (gaussian) view is roughly continuous, library-size
# corrected, and free of a strong mean-variance trend. When those don't hold the
# model still trains — it just spends its first factor on library size or a
# technical axis, which is the single most common way a MOFA result goes quietly
# wrong. This script measures the tell-tale signs on the raw input so they can be
# caught before training rather than diagnosed after.
#
# It cannot know provenance — whether a batch effect was regressed out, what
# normalisation was already applied. It reports what the numbers show and what
# that implies for the model. Which step (if any) to take, and how, is the
# analyst's decision: that depends on the assay and on processing upstream of
# this matrix, neither of which is visible here.
#
# The load-bearing measurement is the library-size check: it correlates each
# sample's total signal with the leading axis of variation. A high correlation on
# raw data means "Factor 1 will be library size unless that effect is removed".
#
# Usage:
#   Rscript check_preprocessing.R --data <file>
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

LIBSIZE_WARN <- 0.7   # |rho| above this: total signal aligns with the leading axis

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

# Correlate the leading axis of variation with per-sample total signal.
# PC1 over samples is cheap: with N samples, the N x N cross-product has the
# sample scores as its leading eigenvector, so cost scales with samples, not
# features. Features are capped and variance-ranked to bound work on wide views.
libsize_alignment <- function(v) {
  n <- ncol(v)
  if (n < 5) return(NA_real_)
  totals <- if (inherits(v, "Matrix")) Matrix::colSums(v, na.rm = TRUE) else colSums(v, na.rm = TRUE)
  keep_rows <- if (nrow(v) > 3000) {
    rv <- if (inherits(v, "Matrix")) {
      m <- Matrix::rowMeans(v); Matrix::rowMeans(v^2) - m^2
    } else apply(v, 1, var, na.rm = TRUE)
    order(rv, decreasing = TRUE)[seq_len(3000)]
  } else seq_len(nrow(v))
  X <- as.matrix(v[keep_rows, , drop = FALSE])
  X[!is.finite(X)] <- NA
  X <- X - rowMeans(X, na.rm = TRUE)          # centre each feature
  X[is.na(X)] <- 0
  C <- crossprod(X)                            # n x n over samples
  pc1 <- eigen(C, symmetric = TRUE)$vectors[, 1]
  suppressWarnings(abs(cor(pc1, totals, method = "spearman")))
}

profile_view <- function(v, name) {
  if (inherits(v, "Matrix")) {
    nz <- v@x; nz <- nz[is.finite(nz)]
    n_total <- prod(dim(v)); sparsity <- (n_total - length(v@x)) / n_total
    vals <- nz                                 # nonzero entries suffice for the checks below
    rng  <- if (length(nz)) range(nz) else c(0, 0); rng[1] <- min(rng[1], 0)
    na_frac <- 0
  } else {
    vals <- as.vector(v); vals <- vals[is.finite(vals)]
    sparsity <- mean(vals == 0); rng <- range(vals)
    na_frac <- mean(is.na(v))
  }
  samp <- if (length(vals) > 1e5) vals[sample(length(vals), 1e5)] else vals
  uniq <- unique(samp)

  is_binary   <- all(uniq %in% c(0, 1))
  nonneg      <- all(vals >= 0)
  all_integer <- all(vals == round(vals))
  has_neg     <- any(vals < 0)

  mv_rho <- NA_real_
  if (nonneg && !is_binary && nrow(v) >= 20) {
    keep <- if (nrow(v) > 2000) sample(nrow(v), 2000) else seq_len(nrow(v))
    sub  <- as.matrix(v[keep, , drop = FALSE])
    fm <- rowMeans(sub, na.rm = TRUE); fv <- apply(sub, 1, var, na.rm = TRUE)
    ok <- is.finite(fm) & is.finite(fv) & fm > 0
    if (sum(ok) > 10) mv_rho <- suppressWarnings(cor(fm[ok], fv[ok], method = "spearman"))
  }

  list(name = name, width = nrow(v), n = ncol(v), na_frac = na_frac,
       sparsity = sparsity, is_binary = is_binary, nonneg = nonneg,
       all_integer = all_integer, has_neg = has_neg, range = rng,
       mv_rho = mv_rho, libsize = libsize_alignment(v))
}

# ---- verdict ---------------------------------------------------------------

verdict <- function(p) {
  notes <- character(0); flag <- "OK"

  raw_counts <- p$nonneg && p$all_integer && !p$is_binary && p$n <= p$width * 10

  if (p$is_binary) {
    notes <- c(notes, "Binary data — no normalisation needed; decide the likelihood with recommend_likelihoods.R.")
  } else if (raw_counts) {
    flag <- "NEEDS WORK"
    notes <- c(notes, sprintf("Looks like RAW COUNTS (non-negative integers, range %g-%g). The default likelihood is gaussian, and depth differences between samples are still in the data. Is this the intended model input?",
                              p$range[1], p$range[2]))
  } else if (p$has_neg) {
    notes <- c(notes, sprintf("Continuous with negative values (range %.1f to %.1f) — consistent with an already-applied transform (log / z-score / M-value / CLR).",
                              p$range[1], p$range[2]))
  } else if (p$range[2] <= 1) {
    flag <- "CHECK"
    notes <- c(notes, "Values bounded in [0,1] — rates or proportions. Gaussian has unbounded support. Is this the intended input scale?")
  } else {
    notes <- c(notes, sprintf("Continuous (range %.2f to %.2f) — looks ready.", p$range[1], p$range[2]))
  }

  if (!is.na(p$mv_rho) && p$mv_rho > 0.5 && !p$is_binary) {
    if (flag == "OK") flag <- "CHECK"
    notes <- c(notes, sprintf("Strong mean-variance trend (rho=%.2f): higher-mean features are more variable. The gaussian likelihood carries a homoscedasticity assumption that a variance-mean relationship conflicts with (MEFISTO preprocessing vignette).",
                              p$mv_rho))
  }

  if (!is.na(p$libsize) && p$libsize >= LIBSIZE_WARN && !p$is_binary) {
    if (flag == "OK") flag <- "CHECK"
    notes <- c(notes, sprintf("Leading axis of variation aligns with per-sample total signal (rho=%.2f). For sequencing data that is the signature of un-removed library size effects: Factor 1 would capture differences in total signal per sample, and more subtle sources of variation would be downweighted (MOFA2 FAQ). Confirm whether library size was already corrected for. If the per-sample total is itself meaningful here (e.g. global drug sensitivity, total abundance), this may instead be real biology — judge by the view.",
                              p$libsize))
  }

  list(flag = flag, notes = notes)
}

# ---- report ----------------------------------------------------------------

cat("\n=== PREPROCESSING CHECK ===\n")
cat("Measured signs of un-preprocessed data, and what they imply for the model.\n")
cat("Cannot see provenance (batch correction, prior normalisation) — confirm those\n")
cat("separately. Which step, if any, a flag calls for is your decision.\n")
cat("Flags: OK / CHECK / NEEDS WORK.\n\n")

for (nm in names(views)) {
  p <- profile_view(views[[nm]], nm)
  v <- verdict(p)
  cat(sprintf("[%s]  %d features x %d samples", nm, p$width, p$n))
  if (p$na_frac > 0) cat(sprintf(", %.0f%% NA", 100 * p$na_frac))
  if (p$sparsity > 0.1) cat(sprintf(", %.0f%% zeros", 100 * p$sparsity))
  cat(sprintf("\n   ==> %s\n", v$flag))
  for (nt in v$notes) cat("   - ", nt, "\n", sep = "")
  cat("\n")
}

cat("Reminder: the two things this cannot measure —\n")
cat("  1. whether batch/technical effects were regressed out (must precede feature selection)\n")
cat("  2. what each view's prior normalisation actually was\n")
cat("still need to be established with the analyst.\n\n")