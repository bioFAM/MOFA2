---
layout: default
title: Troubleshooting
---

### I do not have R 4.0

The R package `MOFA2` works with R>=3. You need to clone the repository:
```
git clone https://github.com/bioFAM/MOFA2
```
and edit the `Depends` option in the `DESCRIPTION` file to your R version. Then install the R package using 
```
R CMD INSTALL MOFA2
```

### I get the following error when running `run_mofa`

```
AttributeError: 'module' object has no attribute 'core.entry_point

Error in py_module_import(module, convert = convert) :
 ModuleNotFoundError: No module named 'mofapy2'
```

This only happens on the default `reticulate` path: `mofapy2` is missing from the Python installation reticulate picked. Restart R (an R session binds to one interpreter, on first use), then either let reticulate provide it,

```r
reticulate::py_require("mofapy2")
```

point it at your own installation with `use_python("YOUR_PYTHON_PATH", required = TRUE)` (check with `py_config()`), or skip the setup with `run_mofa(MOFAobject, use_basilisk = TRUE)`. See the [installation notes](installation.html#notes-on-the-connection-of-r-to-python).



### I get the following error when installing the R package

```
ERROR: dependencies 'XXX', 'YYY' are not available for package 'MOFA2'
```
You probably tried to install them using `install.packages()`. These packages should be installed from Bioconductor.
