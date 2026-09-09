# MEFISTO

MEFISTO extends MOFA with a Gaussian process prior over a **continuous sample
covariate** — time, developmental stage, or spatial coordinates. Instead of
treating samples as exchangeable, it models factors as smooth functions of that
covariate.

What that buys you: factors that vary smoothly along the covariate are favoured
over noisy ones; the model reports *how* smooth each factor is; factor values can
be interpolated to unobserved covariate values; and with several groups, the
covariate axis can be warped to align them.

## When to use it

Use MEFISTO when samples are genuinely *ordered* along a continuous axis and
smoothness along that axis is a meaningful assumption — time courses,
developmental series, spatial transcriptomics, longitudinal microbiome sampling.

Do not use it just because a numeric metadata column exists. A numeric column
with a handful of levels is a grouping variable.

## Setup

```r
obj <- create_mofa(data)
obj <- set_covariates(obj, covariates = "time")        # metadata column
obj <- set_covariates(obj, covariates = c("x", "y"))   # spatial coordinates

mef_opts <- get_default_mefisto_options(obj)
obj <- prepare_mofa(obj, data_options = data_opts, model_options = model_opts,
                    training_options = train_opts, mefisto_options = mef_opts)
obj <- run_mofa(obj, outfile = "mefisto.hdf5", use_basilisk = TRUE)
```

Option table in `options.md`. The two that matter most in practice:

- **`sparseGP = TRUE`** plus a lower `frac_inducing` when there are many distinct
  covariate values. Exact GP inference scales cubically in covariate points, so
  this is the difference between hours and minutes above ~1000 samples.
- **`model_groups`** controls whether covariance structure across groups is
  modelled. Automatically `FALSE` with a single group. Turning it off with
  multiple groups is a useful speed-up when group relationships aren't of
  interest.

## Downstream

### Factors against the covariate

```r
plot_factors_vs_cov(model, factors = "all", color_by = "group")
plot_data_vs_cov(model, factor = 1, view = 1, features = 10, sign = "positive")
```

The first is the main MEFISTO plot: it shows each factor as a function of the
covariate rather than as a distribution over samples. The second drops to the
feature level, showing how the top-weighted features themselves move along the
covariate.

### Variance explained by the covariate

```r
plot_variance_explained_by_covariates(model)
```

Separates variance a factor explains that is attributable to smooth variation
along the covariate from variance that is not. A factor with high total variance
but little covariate-attributable variance is real structure that simply isn't
temporal or spatial.

### Smoothness and sharedness

```r
plot_smoothness(model)        # per factor: how smooth along the covariate
plot_sharedness(model)        # per factor: how shared across groups
get_lengthscales(model)       # underlying GP lengthscales
get_scales(model)
```

**Smoothness** runs 0–1. Near 1 means the factor varies smoothly along the
covariate — genuinely temporal or spatial. Near 0 means it captures variation
unrelated to the covariate, which is not a failure: it is the model telling you
that axis is something else.

**Sharedness** describes whether a factor behaves the same way across groups.
Low sharedness with high smoothness means a process that occurs in each group but
on a different schedule — often the interesting case, and the motivation for
warping.

```r
plot_group_kernel(model)      # covariance structure between groups
```

### Interpolation

```r
model <- interpolate_factors(model, new_values = seq(0, 20, 0.5))
plot_interpolation_vs_covariate(model, covariate = 1, factors = 1,
                                only_mean = FALSE, show_observed = TRUE)
get_interpolated_factors(model)
```

Predicts factor values at unobserved covariate points. Setting `new_values` in
`mefisto_options` *before* training gives uncertainty estimates too;
`interpolate_factors()` after the fact gives means only. If uncertainty matters,
plan for it before training.

Extrapolation beyond the observed range is possible and increasingly unreliable
the further out you go — GP uncertainty grows, which is exactly why
`only_mean = FALSE` is worth plotting.

### Alignment / warping

For multi-group designs where the same process runs on different schedules —
different species' developmental timelines, patients progressing at different
rates:

```r
mef_opts$warping     <- TRUE
mef_opts$warping_ref <- "reference_group"
# then train, and afterwards:
plot_alignment(model)
get_covariates(model, warped = TRUE)
```

Warping learns a monotonic transformation of the covariate per group, aligning
them to the reference. `plot_alignment()` shows the learned mapping and is worth
inspecting directly — an implausible warp (extreme compression, non-monotonic
appearance) means the alignment has latched onto noise.

`warping_open_begin` / `warping_open_end` (both `TRUE` by default) allow groups
to be matched over only part of the range, which is right when groups were
sampled over different windows.

## Available Tutorials that apply MEFISTO

Source at <https://github.com/bioFAM/MEFISTO_tutorials>, rendered from
<https://biofam.github.io/MOFA2/MEFISTO.html>. `MEFISTO_temporal.Rmd` also ships
as a package vignette (`vignette("MEFISTO_temporal", "MOFA2")`).

| Tutorial | Shows |
|---|---|
| `MEFISTO_temporal.Rmd`, `MEFISTO_spatial.Rmd` | Minimal simulated templates — the clean starting point |
| `MEFISTO_temporal_Poisson.Rmd` | Non-gaussian likelihood with covariates |
| `MEFISTO_ST.Rmd` | Real spatial transcriptomics; model trained in Python, loaded in R |
| `evodevo_tutorial.Rmd` | Warping across species, smoothness and sharedness |
| `microbiome_tutorial.Rmd` | Interpolation, weight aggregation to genus level |
| `scnmt_mefisto_vignette.Rmd` | MEFISTO combined with imputation and GSEA |