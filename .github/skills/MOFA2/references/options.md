# MOFA2 options reference

Five option groups, all obtained with a `get_default_*_options(object)` call and
passed to `prepare_mofa()`. Defaults below are read from the package source, not
from the vignettes — see the discrepancy note at the end.

```r
data_opts   <- get_default_data_options(obj)
model_opts  <- get_default_model_options(obj)
train_opts  <- get_default_training_options(obj)
stoch_opts  <- get_default_stochastic_options(obj)   # only if stochastic = TRUE
mef_opts    <- get_default_mefisto_options(obj)      # only with covariates

obj <- prepare_mofa(obj,
                    data_options       = data_opts,
                    model_options      = model_opts,
                    training_options   = train_opts,
                    stochastic_options = stoch_opts,
                    mefisto_options    = mef_opts)
```

## Data options

| Option | Default | When to change |
|---|---|---|
| `scale_views` | `FALSE` | `TRUE` when views have very different variances and you don't want the high-variance one dominating |
| `scale_groups` | `FALSE` | `TRUE` when groups have different ranges/variances |
| `center_groups` | `TRUE` | Leave alone. This is what makes multi-group inference work |
| `use_float32` | `TRUE` | Halves memory at negligible accuracy cost; MOFA2 messages about it above 1e5 samples |

## Model options

| Option | Default | Notes |
|---|---|---|
| `likelihoods` | `"gaussian"` per view | Must be set explicitly for non-gaussian views — auto-inference is disabled. `"poisson"` / `"bernoulli"` rely on Gaussian approximations → `data-prep.md` |
| `num_factors` | sample-size heuristic (below) | The setting users ask about most |
| `spikeslab_factors` | `FALSE` | Sparsity prior on factors |
| `spikeslab_weights` | `FALSE` in code | **Docs say `TRUE` — see note below** |
| `ard_factors` | `FALSE`, auto-`TRUE` if >1 group | Group-wise ARD sparsity |
| `ard_weights` | `TRUE` | View-wise ARD sparsity; how factors become view-specific |

**`num_factors` heuristic**, applied at `get_default_model_options()`:

| N (total samples) | Default |
|---|---|
| ≤ 25 | 5 |
| 26 – 1,000 | 15 |
| 1,001 – 10,000 | 20 |
| > 10,000 | 25 |

There is no correct value; it depends on purpose. For major axes of biological
variation, ~10 usable factors is typical and over-specifying is safe — factors
with near-zero variance explained are visible in `plot_variance_explained()` and
can be dropped with `subset_factors()`. For imputation, ask for more, since small
sources of variation carry real signal there. Only change the sparsity and ARD
priors if you know the underlying model; the defaults are well chosen.

### The package's own guardrails

`prepare_mofa()` prints warnings the user will see and ask about. Know what
triggers each, so you can explain it rather than guess:

| Warning | Condition |
|---|---|
| "very large" factor count, training will be slow | `num_factors > 50` |
| "very small" sample count for this many factors | total N < 4 × `num_factors` |
| little power to learn factors for these group(s) | any group with < 10 samples |
| little power to learn factors for these view(s) | any view with < 15 features |
| recommends more stringent feature selection | any view with > 1e4 features |

Two details in the sample-count warning are worth knowing. Its suggested ceiling
("should not exceed ~N") is computed from the **smallest group** (`min(N)/4`)
while the trigger uses the **total** across groups — so with unbalanced groups
the two numbers refer to different things. And the default heuristic trips its
own check: for a total N between 26 and 59 the default is 15 factors, but the
warning fires below 60 and suggests N/4. Neither is a bug to fix; both are worth
explaining when a user asks why the defaults warned at them.

## Construction idiom: four steps vs `mofa2()`

Default to the explicit four-step form shown at the top of this file. It is what
every published tutorial uses, and it puts the option lists in variables the user
can print, inspect and modify — which is what makes a configuration reviewable.

The `mofa2()` one-liner collapses `create_mofa()` and `prepare_mofa()` and takes
options through `...`:

```r
MOFAobject <- mofa2(data, num_factors = 10)
```

Use it only where the analysis is deliberately compact and runs at or very near
defaults — a quick exploratory pass, a reproducible example, a teaching snippet.
The trade-off is visibility: options set through `...` never become inspectable
objects, so a reader cannot see how the model was configured. Fine for three
defaults, bad for a real analysis with eight tuned settings.

`mofa2()` is a **recent addition**, so it may be absent from an installed
version. Confirm it is available before generating it (`exists("mofa2")` after
loading MOFA2); the four-step form works everywhere.

## Training options

| Option | Default | Notes |
|---|---|---|
| `maxiter` | `1000` | Rarely the binding constraint; convergence usually triggers first |
| `convergence_mode` | `"fast"` | `"fast"` to explore, `"medium"`/`"slow"` for a final model |
| `drop_factor_threshold` | `-1` (off) | Set e.g. `0.01` to drop factors below 1% variance during training |
| `verbose` | `FALSE` | `TRUE` when diagnosing convergence |
| `startELBO` | `1` | First iteration at which the ELBO is computed |
| `freqELBO` | `5` | ELBO computation frequency |
| `stochastic` | `FALSE` | See below |
| `gpu_mode` | `FALSE` | Requires cupy + CUDA |
| `gpu_device` | `NULL` | Which device, when several |
| `seed` | `42` | Change it to test robustness across runs |
| `outfile` | `NULL` | Always set this — training is the expensive step |
| `weight_views` | `FALSE` | `TRUE` weights the ELBO by view width; a corrective for imbalance |
| `save_interrupted` | `FALSE` | `TRUE` keeps a partial model if training is interrupted |

## Stochastic options

Only read when `train_opts$stochastic <- TRUE`. Worth it for very large sample
counts (roughly >10k), and most useful alongside a GPU.

| Option | Default | Notes |
|---|---|---|
| `batch_size` | `0.5` | Fraction of samples per update; lower = faster, noisier |
| `learning_rate` | `1.0` | Starting rate |
| `forgetting_rate` | `0.5` | Decay of the learning rate |
| `start_stochastic` | `1` | First stochastic iteration |

Stochastic inference reintroduces run-to-run randomness that standard variational
inference does not have (MOFA2 initialises factors by PCA on the concatenated
data, which makes standard VI deterministic). Compare against a standard run
before trusting SVI results on a new dataset.

## MEFISTO options

Only relevant once `set_covariates()` has been called. See `mefisto.md` for the
analysis workflow.

| Option | Default | Notes |
|---|---|---|
| `scale_cov` | `FALSE` | Scale covariates |
| `start_opt` | `20` | First iteration to optimise GP hyperparameters |
| `n_grid` | `20` | Grid points for the hyperparameter search |
| `opt_freq` | `10` | Optimisation frequency |
| `model_groups` | `TRUE`, auto-`FALSE` if 1 group | Model covariance structure across groups |
| `sparseGP` | `FALSE` | `TRUE` for many covariate points — exact GP inference scales cubically |
| `frac_inducing` | `0.75` | Fraction of samples used as inducing points; lower = faster, less exact |
| `warping` | `FALSE` | Align covariates across groups; needs a multi-group design |
| `warping_freq` | `20` | Optimisation frequency for warping |
| `warping_ref` | first group | Reference group to align to |
| `warping_open_begin` / `_open_end` | `TRUE` | Allow unmatched start/end of the covariate range |
| `warping_groups` | `NULL` | Custom grouping for warping |
| `new_values` | `NULL` | Covariate values to interpolate/extrapolate to during training |

## Known documentation discrepancy

`spikeslab_weights` is set to `FALSE` in
`R/prepare_mofa.R` (`get_default_model_options`), but the roxygen block above it —
and therefore `man/get_default_model_options.Rd` and the
`getting_started_R.Rmd` vignette — states the default is `TRUE`.

Trust the code. If a user's model behaves unlike the documentation suggests, this
is a candidate explanation. Worth flagging to them so they set the option
explicitly rather than relying on either source.