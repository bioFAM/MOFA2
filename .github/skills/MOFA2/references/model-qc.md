# Model QC

Training always succeeds. Convergence tells you the optimiser stopped, not that
the model is meaningful. Run these checks before interpreting anything, and read
the variance plot with the user rather than for them — most of the diagnosis is
in what they know about their data.

## The standard pass

```r
model <- load_model("model.hdf5")

get_dimensions(model)                                    # N, D, M, G, K
plot_variance_explained(model, x = "view", y = "factor") # the central diagnostic
plot_variance_explained(model, plot_total = TRUE)        # total per view
plot_factor_cor(model)                                   # factor redundancy
```

`plot_variance_explained()` takes `x`, `y` and `split_by` from
`{"view", "factor", "group"}`, so multi-group models can be laid out as
`x = "group", y = "factor", split_by = "view"`.

## Reading the variance plot

A healthy model shows *structure*: some factors explaining variance across
several views (shared axes), others concentrated in one (view-specific biology),
and a tail of near-zero factors that can be dropped.

| What you see | Usual cause | What to do |
|---|---|---|
| Factor 1 dominates *one* view, and that view's leading axis tracks per-sample totals | Depth / library-size effect in that view | Report it and check that view's input state — `data-prep.md`, `check_preprocessing.R` |
| Factor 1 dominates across *all* views | A shared sample-level axis — often real biology, sometimes something that affected every assay (batch, sample quality) | Not diagnostic by itself: correlate it with metadata before concluding anything |
| A dominant factor tracks a known technical variable | That structure is still in the data — MOFA has no batch term, so it is free to become a factor | Report it; if the user removes it, that has to precede feature selection, and the model needs retraining |
| One view holds nearly all variance | View width imbalance | `weight_views = TRUE` compensates in-model; reducing the wide view is the analyst's call |
| Flat profile, nothing above a few percent | Over-filtered data, too few samples, or no shared structure | Reconsider filtering; possibly MOFA is the wrong tool |
| Many factors, all tiny | `num_factors` too high | Harmless — drop them with `subset_factors()` |
| Two factors highly correlated in `plot_factor_cor` | Redundancy, often from a strong non-linear axis split across factors | Usually fine to interpret jointly; consider fewer factors |

The last row deserves care. Because MOFA is linear, a strong non-linear
trajectory often appears as two or three correlated factors describing one
underlying process. That is a property of the data, not a bug — but it should be
described as one process, not three.

## Sanity checks against the input

The most valuable check is whether a factor is visible in the raw data:

```r
plot_data_heatmap(model, factor = 1, view = "mRNA", features = 25,
                  denoise = TRUE, cluster_rows = TRUE)
plot_data_scatter(model, factor = 1, view = "mRNA", features = 6,
                  add_lm = TRUE, color_by = "condition")
```

If the top-weighted features show no visible gradient along the factor in the
observed data, the factor is not describing something real about the input,
whatever the variance numbers say.

## Robustness across seeds

Standard variational inference in MOFA2 is deterministic — factors are
initialised by PCA on the concatenated data — so repeated runs with the same
options give the same answer. Robustness testing therefore means varying
something: the seed under *stochastic* inference, the feature set, or a
subsample of the data.

```r
models <- lapply(c(1, 2, 3), function(s) {
  train_opts$seed <- s
  run_mofa(prepare_mofa(obj, training_options = train_opts,
                        data_options = data_opts, model_options = model_opts),
           outfile = sprintf("model_seed%d.hdf5", s), use_basilisk = TRUE)
})

compare_factors(models)     # correlation of factors across models
compare_elbo(models)        # final ELBO per model
best <- select_model(models)
```

`compare_factors()` is the useful one: factors recovered consistently across
perturbations are the ones worth building a story on. Note that a factor which
disappears under downsampling is not necessarily an artefact — small but real
sources of variation need the full data to be detectable.

## Sign and scale, when reporting

Two things users reliably over-read, worth stating whenever showing results:

- **Factor values** are meaningful only relatively. The ordering of samples and
  the sign relative to other samples carry information; the absolute number does
  not. Interpretation is the same as for a principal component.
- **Signs are arbitrary per run.** A factor and all its weights can flip
  wholesale without changing the solution, because the model is rotation- and
  reflection-invariant. Compare signs within one model, never across models.

## Dropping and renaming

```r
r2 <- get_variance_explained(model)$r2_per_factor[[1]]
keep <- which(rowMeans(r2) > 1)                  # >1% variance on average
model <- subset_factors(model, keep)

factors_names(model) <- paste0("Factor", seq_along(keep))
views_names(model)   <- c("RNA", "Methylation")
```