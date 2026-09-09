# Building the MOFA object

`create_mofa(data, assays = NULL, groups = NULL, extract_metadata = TRUE, ...)`
dispatches on the class of `data`. The format-specific constructors can also be
called directly.

Whatever the route, follow it with `plot_data_overview(obj)` and read the result
with the user — sample overlap and view widths become visible there while they
are still cheap to fix.

## Naming rules, enforced at construction

The errors these raise are cryptic, so check them up front:

- View and group names must not contain `/`.
- Sample and feature names must be **ASCII** (accented characters in sample IDs
  are a recurring cause of failure) and **unique within a view**.
- Features appearing in several views get a view suffix appended automatically.
  Strip it when merging weights back against external annotation.

## List of matrices

The simplest route. **Features in rows, samples in columns** — the transpose of
what most single-cell tooling uses.

```r
data <- list(mRNA = mrna_mat, meth = meth_mat, drugs = drug_mat)
obj  <- create_mofa(data)

# multi-group: one group label per sample, in column order
obj  <- create_mofa(data, groups = metadata$condition)
```

Column names must be consistent across views; that is how samples are matched.
Samples absent from a view can simply be absent from its columns.

## Long data.frame

Best for complex designs, and for sparse data where a dense matrix would be
wasteful. Rows for missing observations can simply be omitted.

```r
# columns: sample, feature, view, value  (+ optional group)
obj <- create_mofa_from_df(df, extract_metadata = TRUE)
```

`extract_metadata = TRUE` pulls any additional columns into sample metadata,
which saves attaching them later.

## MultiAssayExperiment

```r
create_mofa_from_MultiAssayExperiment(mae,
                                      experiments      = NULL,
                                      alt_experiments  = NULL,
                                      assays           = NULL,
                                      groups           = NULL,
                                      extract_metadata = FALSE)
```

`experiments` selects which experiments become views (by name or index);
`assays` picks which assay within each. Note `extract_metadata` defaults to
`FALSE` here, unlike `create_mofa_from_df` — pass `TRUE` to carry `colData`
across, otherwise metadata must be attached afterwards.

## SingleCellExperiment

```r
create_mofa_from_SingleCellExperiment(sce,
                                      alt_experiments = c("main"),
                                      assays          = c("logcounts"),
                                      groups          = NULL,
                                      extract_metadata = FALSE)
```

`alt_experiments` is how a multi-modal SCE becomes multi-view: `"main"` refers to
the primary experiment, other entries to `altExp()` slots (ADT, ATAC…). Defaults
to `logcounts`, so normalise first.

## Seurat

```r
create_mofa_from_Seurat(seurat,
                        groups   = NULL,
                        assays   = NULL,
                        layer    = "data",
                        features = NULL,
                        extract_metadata = FALSE)
```

`assays` names the Seurat assays to use as views. `layer = "data"` takes the
normalised layer — this is what you usually want; `"counts"` would feed raw
counts to a gaussian likelihood. `features` restricts to a feature set, which is
the natural place to apply per-view HVG selection.

`groups` may name a column in the Seurat object's metadata.

## Attaching metadata

If not extracted at construction:

```r
samples_metadata(obj) <- metadata_df    # must contain a `sample` column
```

The `sample` column must match the sample names in the object exactly. Attach
metadata before downstream analysis so `plot_factor(color_by = ...)` and
`correlate_factors_with_covariates()` can see it.

## MEFISTO covariates

```r
obj <- set_covariates(obj, covariates = "time")            # metadata column
obj <- set_covariates(obj, covariates = c("x", "y"))       # spatial coordinates
obj <- set_covariates(obj, covariates = cov_matrix)        # covariates x samples
```

A supplied matrix has covariates in rows and samples in columns, with column
names matching the object's samples.

## Going the other way

```r
seurat <- add_mofa_factors_to_seurat(mofa_obj, seurat_obj)   # factors as a reduction
```

Puts the factors back into a Seurat object as a dimensional reduction, so
downstream Seurat workflows (clustering, UMAP) can run on the MOFA latent space.