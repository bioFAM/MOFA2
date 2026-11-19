## Changes in version 1.24.0

### New features

* `create_mofa_from_SingleCellExperiment()` can now load `altExps` as additional views
  (`alt_experiments`).
* `create_mofa_from_MultiAssayExperiment()` can now select experiments, assays and altExps
  (`experiments`, `assays`, `alt_experiments`).
* New `mofa2()` wrapper runs `create_mofa()` and `prepare_mofa()` in one call, taking only the
  options that differ from the defaults.
* New `plot_factor_mean_cor()` heatmap to identify factors related to intercept or library-size
  effects.

### Significant user-visible changes

* `create_mofa()` and the SingleCellExperiment and MultiAssayExperiment constructors have new
  arguments before `groups`; pass `groups` by name.
* `create_mofa_from_SingleCellExperiment()` argument `assay` is renamed to `assays`.
* Documentation and vignettes now recommend `reticulate::py_require()` to set up mofapy2.

### Bug fixes and minor improvements

* Updated to mofapy2 0.7.5, which fixes MEFISTO for some covariate grids and adds pandas 3
  compatibility.
* The warning about intercept factors now names the affected views and points to
  `plot_factor_mean_cor()`.
* `run_mofa()` no longer warns when mofapy2 is a newer patch release than the pinned version.
