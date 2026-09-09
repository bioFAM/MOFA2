# Troubleshooting and scale-up

MOFA2's R package is a front end; training happens in the Python package
`mofapy2`, reached through reticulate. Nearly all installation problems are about
that bridge, not about MOFA itself.

There is no OS-specific configuration in MOFA2 — the package contains no
platform branching, and reticulate or basilisk absorbs the differences. So
diagnose by *environment* (RStudio vs. terminal vs. cluster vs. container), not
by operating system, and don't invent macOS- or Windows-specific advice.

## Connecting R to Python

Three routes make `mofapy2` available, in the order the MOFA2 documentation
recommends them:

**1. Let reticulate provision it** (simplest; needs reticulate ≥ 1.41). Declare
the dependency before training and reticulate builds an isolated environment,
downloading an interpreter if none is suitable:

```r
reticulate::py_require("mofapy2", python = "3.12")
model <- run_mofa(obj, outfile = "model.hdf5")
```

This is also how to pin an exact version, which the version check below cares
about: `py_require("mofapy2==0.7.5")`.

**2. Let basilisk handle it** — a dedicated, version-pinned environment built on
first use:

```r
model <- run_mofa(obj, outfile = "model.hdf5", use_basilisk = TRUE)
```

Note it defaults to `FALSE`, so it must be passed explicitly. Because basilisk
runs Python in a *separate process*, this is the route to choose when MOFA2 is
used alongside other R packages with conflicting Python dependencies. Its
environment is pinned to `mofapy2 0.7.5` with python 3.12.10, numpy 1.26.4,
scipy 1.12.0, pandas 2.2.1, h5py 3.10.0, scikit-learn 1.4.0, dtw-python 1.3.1 —
worth knowing when results or speed differ from a hand-built environment. It
does **not** contain cupy, so GPU training cannot use it (see below).

**3. Use an existing Python installation** — the most configuration. Point
reticulate at the interpreter *before* loading MOFA2, since reticulate binds on
first use and will not switch afterwards:

```r
library(reticulate)
use_python("/path/to/python", required = TRUE)   # or use_condaenv("myenv")
library(MOFA2)
py_config()                                      # confirm the binding
py_module_available("mofapy2")                   # confirm mofapy2 is visible
```

If an interpreter is already bound, no call will move it — **restart R**. Both
the website and the package's own error messages say "in a new R session" for
this reason.

## `run_mofa()` fails

| Symptom | Cause | Fix |
|---|---|---|
| `ModuleNotFoundError: mofapy2` | Wrong interpreter, or not installed | `py_require("mofapy2")` in a fresh session, or install into the interpreter reticulate actually uses |
| "mofapy2 X is a major/minor version behind/ahead of 0.7.5" — an **error** | Version check against the pinned target | `py_require("mofapy2==0.7.5")` in a new session, or `use_basilisk = TRUE` |
| "Could not parse the mofapy2 version … Only numeric versions are supported" | A PEP 440 dev version such as `0.7.4.dev0` | Expected when installing mofapy2 from GitHub; install a released version or use basilisk |
| Works in terminal, fails in R/RStudio | RStudio picked a different Python | `py_config()`, then bind explicitly and restart |
| Session crashes with no error message | Can be basilisk itself, a numpy/BLAS conflict, or memory exhaustion | Try *both* directions: with and without basilisk; single-thread (below); or reduce data size |
| Fails only on a cluster | Interactive vs. batch environment differ | Set `RETICULATE_PYTHON` in the job script |

The version check is strict by design: `run_mofa()` **errors** on any major or
minor mismatch against the pinned `mofapy2` version, warns when a patch behind,
and merely notes a patch ahead.

For a crash with no message, the official troubleshooting page suggests
disabling multi-threading — note this is the opposite of the speed advice below,
and applies to crashes rather than slowness:

```r
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
```

**Non-ASCII names**: `create_mofa()` rejects non-ASCII sample or feature names,
and the message is not obvious if you don't know to look for it. Accented
characters in patient IDs are the usual source.

**`/` in view or group names**: also rejected, because names map onto HDF5 paths.

```r
views_names(obj)  <- gsub("/", "-", views_names(obj))
groups_names(obj) <- gsub("/", "-", groups_names(obj))
```

## Training runs but is very slow

In rough order of impact:

**1. Feature count.** It drives runtime more than anything else, and reducing a
wide view also rebalances it against the others — but which features to keep is
the analyst's call. Report the widths and the runtime consequence.

**2. Check BLAS threading.** MOFA's numerical work is numpy; if numpy isn't
linked against a threaded BLAS, it runs on one core.

```bash
python -c "import numpy; numpy.show_config()"
export MKL_NUM_THREADS=8        # Intel MKL
export OPENBLAS_NUM_THREADS=8   # OpenBLAS
```

Anaconda's numpy is linked to MKL by default.

**3. Loosen convergence while exploring.** `convergence_mode = "fast"` is the
default for good reason; save `"medium"`/`"slow"` for the final model.

**4. Stochastic inference**, for very large sample counts (roughly >10k):

```r
train_opts$stochastic <- TRUE
stoch_opts <- get_default_stochastic_options(obj)
stoch_opts$batch_size <- 0.25        # smaller = faster, noisier
```

SVI reintroduces run-to-run randomness that standard VI does not have, so compare
against a standard run on a subsample before trusting it on new data.

**5. `use_float32 = TRUE`** halves memory at negligible accuracy cost. Already
the default, and MOFA2 messages about it above 1e5 samples.

## GPU training

Needs an NVIDIA GPU and `cupy` matching the local CUDA version. The basilisk
environment has no cupy, so this route requires `use_basilisk = FALSE` —
otherwise `gpu_mode = TRUE` is silently ineffective.

```bash
pip install cupy-cuda12x     # match your CUDA
python -c "import cupy; print(cupy.cuda.runtime.getDeviceCount())"
```

or when using reticulate:

```r
reticulate::py_require(c("mofapy2==0.7.5", "cupy-cuda12x"), python_version = "3.12")
  # or cupy-cuda13x
```

```r
train_opts$gpu_mode   <- TRUE
train_opts$gpu_device <- 0       # if several
```

Worth it for large models, especially with stochastic inference. Not worth
setting up for a few hundred samples — CPU training there is minutes. The bundled
`scripts/profile_mofa_data.R` probes for cupy and only raises GPU as an option
when it is actually available.

## Interrupted training

```r
train_opts$save_interrupted <- TRUE
```

Keeps a partial model if training is stopped. Worth setting for long runs on
shared or time-limited compute.

## Reloading a model

```r
model <- load_model("model.hdf5", load_data = TRUE, sort_factors = TRUE)
model <- load_model("model.hdf5", on_disk = TRUE)     # large data, keep on disk
```

`on_disk = TRUE` uses HDF5Array so the data is not held in memory — necessary for
very large models, at the cost of slower access.

If a model was trained in Python without feature names, add them back after
loading (the MEFISTO spatial-transcriptomics tutorial does this):

```r
features_names(model) <- feature_name_list
samples_names(model)  <- sample_name_list
```

## Escalation

When a problem isn't covered here: the MOFA2 FAQ
(<https://biofam.github.io/MOFA2/faq.html>), the troubleshooting page
(<https://biofam.github.io/MOFA2/troubleshooting.html>), and the issue tracker at
<https://github.com/bioFAM/MOFA2>. Include `sessionInfo()` and `py_config()`
output in any report — installation issues are almost never reproducible without
both.