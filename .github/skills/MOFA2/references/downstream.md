# Downstream analysis

Once the model passes QC, characterisation is cheap and reversible — this is
where to be proactive and suggest the next step unprompted.

## Attach metadata first

```r
samples_metadata(model) <- metadata_df    # needs a `sample` column
```

Everything below that colours or correlates by a variable depends on this.

## The per-factor loop

Work one factor at a time. The order matters: establish *where* a factor lives
before asking what it means, and check the raw data before believing either.

### 1. Where does it live?

```r
plot_variance_explained(model, x = "view", y = "factor")
calculate_variance_explained(model)                        # numeric
calculate_variance_explained_per_sample(model)             # per-sample
```

A factor explaining variance in several views is a shared axis; one confined to a
single view is modality-specific biology (or a modality-specific artefact).

### 2. What does it track?

```r
correlate_factors_with_covariates(model,
                                  covariates = c("condition", "age", "batch"),
                                  plot = "log_pval")

plot_factor(model, factors = 1, color_by = "condition", add_violin = TRUE,
            dodge = TRUE, add_boxplot = TRUE)

summarise_factors(model, df = metadata_df)    # mean factor value per level
```

Check technical variables here too, not just biological ones. A factor
correlating with `batch` at this stage is the cheapest bad news you will get all
analysis.

### 3. What drives it?

```r
plot_top_weights(model, view = "mRNA", factors = 1, nfeatures = 20)
plot_weights(model, view = "mRNA", factors = 1, nfeatures = 10,
             text_size = 4, scale = TRUE)

plot_weights_heatmap(model, view = "mRNA", factors = 1:5)
plot_weights_scatter(model, factors = c(1, 2), view = "mRNA")
```

Weight interpretation: magnitude is the strength of association with the factor;
sign is the direction. A positive weight means the feature is higher in samples
with positive factor values. Features near zero are unrelated to that factor.

### 4. Does the raw data agree?

```r
plot_data_heatmap(model, factor = 1, view = "mRNA", features = 25,
                  denoise = TRUE, cluster_rows = TRUE, cluster_cols = FALSE)

plot_data_scatter(model, factor = 1, view = "mRNA", features = 6,
                  add_lm = TRUE, color_by = "condition")
```

Do not let this step be skipped. It is what separates a factor describing real
covariation from one that merely has plausible-looking weights. `denoise = TRUE`
shows the model's reconstruction rather than raw values — useful for seeing the
pattern, but always look at the un-denoised version too, or you are checking the
model against itself.

### 5. What does it mean biologically?

`run_enrichment()` — see `enrichment.md`.

## Factor space as an embedding

```r
plot_factors(model, factors = c(1, 2), color_by = "condition", dot_size = 2)

model <- run_umap(model, n_neighbors = 30, min_dist = 0.3)
model <- run_tsne(model)
plot_dimred(model, method = "UMAP", color_by = "cell_type")

clusters <- cluster_samples(model, k = 5, factors = "all")
```

Running UMAP on MOFA factors rather than raw data is one of the main practical
uses of the model for single-cell work — the factors have already integrated the
modalities and removed most of the noise.

## Extracting for downstream use

```r
Z <- get_factors(model, factors = "all", as.data.frame = FALSE)   # list per group
W <- get_weights(model, views = "all", as.data.frame = FALSE)     # list per view
D <- get_data(model, as.data.frame = TRUE)                        # long format
```

`as.data.frame = TRUE` returns long format throughout, which is usually easier
for custom ggplot work.

Factors make good features for supervised models — low-dimensional, decorrelated
and interpretable. The CLL tutorial predicts clinical subgroups this way:

```r
Z <- get_factors(model)[[1]]
df <- data.frame(Z, outcome = metadata$outcome)
fit <- glm(outcome ~ ., data = df, family = "binomial")
```

Guard against circularity: if the factor was selected because it correlated with
the outcome, a model built on it is not independent evidence.

## Imputation and prediction

```r
model <- impute(model, views = "all", factors = "all")
imputed <- get_imputed_data(model, views = "mRNA")

pred <- predict(model, views = "all", factors = "all", add_intercept = TRUE)
```

`impute()` fills missing values from the factor decomposition; `predict()`
reconstructs the full data matrix from factors and weights. For imputation
specifically, train with more factors than you would for interpretation — small
sources of variation matter there.

## Per-sample view contributions

```r
model <- calculate_contribution_scores(model, views = "all", scale = TRUE)
```

Scores how much each view contributes to each sample's position in factor space —
useful for spotting samples driven by a single modality.

Note: `get_contribution_scores()` and `plot_contribution_scores()` exist in the
source but are **not exported**, so the scores must be read out of the object's
metadata directly after calling `calculate_contribution_scores()`.

## Subsetting

```r
model <- subset_factors(model, factors = 1:5)
model <- subset_views(model, views = c("mRNA", "meth"))
model <- subset_groups(model, groups = "group1")
model <- subset_samples(model, samples = keep_samples)
model <- subset_features(model, view = "mRNA", features = keep_features)
```

## Plot customisation

Every `plot_*` function returns a ggplot2 object, so styling is ordinary ggplot:

```r
plot_factor(model, factors = 1, color_by = "condition") +
  scale_fill_manual(values = my_palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Factor 1 by condition")
```

Several functions also take `return_data = TRUE`, which gives the underlying
data frame for fully custom figures.