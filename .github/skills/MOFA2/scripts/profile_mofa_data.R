#!/usr/bin/env Rscript

# Profile multi-omics data and propose a MOFA2 configuration.
#
# The point of this script is to derive the settings that data shape determines
# on its own, so the analyst is only asked about things that genuinely require
# judgement (normalisation state, batch structure, the biological question).
# Every proposal is printed with the observation that produced it, so nothing is
# taken on trust.
#
# Usage:
#   Rscript profile_mofa_data.R --data <file> [--metadata <file>] [--groups <col>]
#
#   --data      .rds holding a named list of matrices (features x samples),
#               an untrained or trained MOFA object, or a .hdf5 MOFA model
#   --metadata  optional .csv/.tsv sample metadata; needs a `sample` column
#   --groups    optional metadata column defining groups
#   --no-gpu    skip the cupy/GPU probe (it imports Python, which can be slow)

suppressWarnings(suppressMessages({
  args <- commandArgs(trailingOnly = TRUE)
}))

if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
  self  <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
  lines <- readLines(self, warn = FALSE)[-1]          # drop the shebang
  start <- which(grepl("^#", lines))[1]               # skip the blank line after it
  end   <- start + which(!grepl("^#", lines[start:length(lines)]))[1] - 2
  cat(sub("^# ?", "", lines[start:end]), sep = "\n")
  quit(status = 0)
}

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}

data_path  <- get_arg("--data")
meta_path  <- get_arg("--metadata")
group_col  <- get_arg("--groups")
probe_gpu  <- !("--no-gpu" %in% args)

if (is.null(data_path)) stop("--data is required. Use --help for usage.", call. = FALSE)

# ---------------------------------------------------------------- load data --

# Normalises whatever was supplied into a named list of feature x sample matrices.
load_views <- function(path) {
  if (grepl("\\.hdf5$|\\.h5$", path, ignore.case = TRUE)) {
    if (!requireNamespace("MOFA2", quietly = TRUE)) stop("MOFA2 needed to read .hdf5 models")
    obj <- MOFA2::load_model(path, load_data = TRUE)
    return(list(views = lapply(MOFA2::get_data(obj), function(v) do.call(cbind, v)), object = obj))
  }
  x <- readRDS(path)
  if (inherits(x, "MOFA")) {
    return(list(views = lapply(MOFA2::get_data(x), function(v) do.call(cbind, v)), object = x))
  }
  if (is.list(x) && all(vapply(x, function(v) is.matrix(v) || inherits(v, "Matrix"), logical(1)))) {
    if (is.null(names(x))) names(x) <- paste0("view", seq_along(x))
    return(list(views = x, object = NULL))
  }
  stop("--data must be a list of matrices, a MOFA object, or a .hdf5 model", call. = FALSE)
}

loaded <- load_views(data_path)
views  <- loaded$views

metadata <- NULL
if (!is.null(meta_path)) {
  sep <- if (grepl("\\.tsv$|\\.txt$", meta_path, ignore.case = TRUE)) "\t" else ","
  metadata <- utils::read.delim(meta_path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE)
}

# ------------------------------------------------------------- observations --

D <- vapply(views, nrow, numeric(1))                       # features per view
all_samples <- unique(unlist(lapply(views, colnames)))
N <- length(all_samples)
M <- length(views)

overlap <- vapply(views, function(v) sum(all_samples %in% colnames(v)) / N, numeric(1))
missing_frac <- vapply(views, function(v) {
  if (inherits(v, "Matrix")) return(NA_real_)
  mean(is.na(v))
}, numeric(1))

# What kind of numbers are in each view? Drives the likelihood conversation.
classify_view <- function(v) {
  vals <- if (inherits(v, "Matrix")) v@x else as.vector(v)
  vals <- vals[is.finite(vals)]
  if (!length(vals)) return("empty")
  u <- unique(vals[seq_len(min(length(vals), 1e5))])
  if (all(u %in% c(0, 1)))                      return("binary")
  if (all(u >= 0) && all(u == round(u)))        return("counts")
  if (all(u >= 0) && all(u <= 1))               return("proportions/rates")
  "continuous"
}
kinds <- vapply(views, classify_view, character(1))

# Features never observed anywhere contribute nothing and trigger warnings later.
dead_features <- vapply(views, function(v) {
  if (inherits(v, "Matrix")) return(NA_real_)
  sum(rowMeans(is.na(v)) == 1)
}, numeric(1))

cat("\n=== DATA PROFILE ===\n\n")
cat(sprintf("%d views, %d samples total\n\n", M, N))
prof <- data.frame(
  view          = names(views),
  features      = as.integer(D),
  samples       = vapply(views, ncol, integer(1)),
  sample_overlap = sprintf("%.0f%%", overlap * 100),
  missing       = ifelse(is.na(missing_frac), "-", sprintf("%.1f%%", missing_frac * 100)),
  data_kind     = kinds,
  row.names     = NULL
)
print(prof, right = FALSE)

# ------------------------------------------------- candidate MEFISTO covars --

covariate_candidates <- character(0)
if (!is.null(metadata)) {
  hint <- "time|day|week|month|age|stage|pseudotime|hour|depth|dist|coord|^x$|^y$|row|col"
  for (nm in setdiff(names(metadata), "sample")) {
    col <- metadata[[nm]]
    if (!is.numeric(col)) next
    n_distinct <- length(unique(col[!is.na(col)]))
    # A continuous covariate needs enough distinct values to define smoothness
    # over; a numeric column with 3 levels is a grouping variable, not an axis.
    if (n_distinct >= 5 && (n_distinct > 0.1 * N || grepl(hint, nm, ignore.case = TRUE))) {
      covariate_candidates <- c(covariate_candidates, nm)
    }
  }
  cat("\nMetadata: ", ncol(metadata), " columns, ", nrow(metadata), " rows\n", sep = "")
  if (length(covariate_candidates)) {
    cat("Candidate continuous covariates: ", paste(covariate_candidates, collapse = ", "), "\n", sep = "")
  }
  if (!is.null(group_col) && group_col %in% names(metadata)) {
    cat("Groups from `", group_col, "`: ", length(unique(metadata[[group_col]])), " levels\n", sep = "")
  }
}

# ------------------------------------------------------------- GPU / cupy ----

gpu_available <- FALSE
if (probe_gpu && requireNamespace("reticulate", quietly = TRUE)) {
  gpu_available <- tryCatch({
    cp <- reticulate::import("cupy", delay_load = FALSE)
    cp$cuda$runtime$getDeviceCount() > 0
  }, error = function(e) FALSE)
}

# ------------------------------------------------------------- proposals -----

cat("\n=== PROPOSED CONFIGURATION ===\n")
cat("Each line states the observation first, then what it implies.\n")
cat("These are proposals, not decisions — review before use.\n\n")

say <- function(setting, why) cat(sprintf("  %-42s %s\n", setting, why))

# Sample size gates everything downstream.
if (N < 15) {
  cat("!! N = ", N, " samples. MOFA needs ~15+ to decompose variance meaningfully.\n",
      "   Consider whether factor analysis is the right tool here at all.\n\n", sep = "")
}

num_factors <- if (N <= 25) 5 else if (N <= 1e3) 15 else if (N <= 1e4) 20 else 25
say(sprintf("model_opts$num_factors <- %d", num_factors),
    sprintf("# N = %d; package heuristic. Adjust for purpose:", N))
cat(sprintf("  %-42s %s\n", "", "#   exploration ~10 usable; imputation: more"))

if (N > 1e4) {
  say("train_opts$stochastic <- TRUE",
      sprintf("# N = %d; full-batch updates get expensive", N))
  say("stoch_opts <- get_default_stochastic_options(obj)",
      "# then tune batch_size / learning_rate")
}
if (N > 1e5) {
  say("data_opts$use_float32 <- TRUE", "# N > 1e5; MOFA2 sets this itself, stated for clarity")
}

# View width imbalance is the quiet killer: the widest view wins the factors.
if (M > 1 && max(D) / min(D) > 10) {
  say(sprintf("# view widths %s", paste(sprintf("%s=%d", names(D), D), collapse = ", ")),
      "")
  say("train_opts$weight_views <- TRUE",
      sprintf("# widest/narrowest = %.0fx; wide views dominate factors", max(D) / min(D)))
  cat("     Reducing the wide view is the alternative lever, but which features\n",
      "     to keep is a data-specific call.\n", sep = "")
}

# Likelihood: gaussian unless there is a real reason, because the alternatives
# are approximations rather than exact inference.
for (i in seq_len(M)) {
  if (kinds[i] == "counts") {
    say(sprintf("# view '%s' looks like raw counts", names(views)[i]),
        "")
    cat("     Confirm whether counts are the intended model input, and set the\n",
        "     likelihood with recommend_likelihoods.R.\n", sep = "")
  } else if (kinds[i] == "binary") {
    say(sprintf("model_opts$likelihoods['%s'] <- 'bernoulli'", names(views)[i]),
        "# view is 0/1")
  }
}

if (length(covariate_candidates)) {
  cat("\n  MEFISTO is worth considering:\n")
  cat("     Metadata column(s) ", paste(sprintf("`%s`", covariate_candidates), collapse = ", "),
      " look continuous.\n", sep = "")
  cat("     If samples are ordered along one of these (time, developmental stage,\n",
      "     spatial position), MEFISTO models smooth variation along it:\n", sep = "")
  say("obj <- set_covariates(obj, covariates = '<column>')", "")
  if (N > 1000) {
    say("mefisto_opts$sparseGP <- TRUE",
        sprintf("# N = %d; exact GP scales cubically in covariate points", N))
    say("mefisto_opts$frac_inducing <- 0.5", "# lower = faster, less exact")
  }
  cat("     This changes the whole analysis — confirm with the analyst before adopting.\n")
}

if (!is.null(group_col) && !is.null(metadata) && group_col %in% names(metadata)) {
  cat("\n  Multi-group requested. Confirm the question being asked:\n")
  cat("     Multi-group centers features WITHIN each group, so it cannot find a\n",
      "     factor separating the groups. Use it to ask which axes of variation are\n",
      "     shared vs. group-specific. For group separation, use plain MOFA and test\n",
      "     factors against the label afterwards.\n", sep = "")
}

if (gpu_available) {
  say("train_opts$gpu_mode <- TRUE", "# cupy import succeeded and a CUDA device is visible")
} else if (probe_gpu && (N > 1e4)) {
  cat("\n  No usable GPU detected (cupy import failed or no CUDA device).\n")
  cat("     Large model — a GPU would help. See references/troubleshooting.md.\n")
}

if (any(!is.na(dead_features) & dead_features > 0)) {
  cat("\n  Features with no observation in any sample:\n")
  for (i in which(!is.na(dead_features) & dead_features > 0)) {
    cat(sprintf("     %s: %d features — drop them before create_mofa()\n",
                names(views)[i], dead_features[i]))
  }
}

cat("\n=== STILL NEEDS A HUMAN ===\n")
cat("  Data shape cannot answer these; ask the analyst:\n")
cat("   1. What normalisation has each view already had?\n")
cat("   2. Are there batch/technical covariates, and were they regressed out\n")
cat("      BEFORE feature selection?\n")
cat("   3. What is the model for — biological axes, imputation, or features for\n")
cat("      a downstream classifier?\n\n")