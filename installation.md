---
layout: default
title: Installation
---

The core of MOFA is implemented in the Python package `mofapy2`, but we recommend to use the R package `MOFA2` which provides an interface to train a MOFA model with R and run the downstream analysis and takes care of setting up all python dependencies. Alternatively, if you prefer to use Python the package `mofax` can be used for downstream analysis in Python, see also our FAQ section.

## Stable release (easiest)

You can install the stable release from Bioconductor ([link](http://www.bioconductor.org/packages/release/bioc/html/MOFA2.html)): 
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MOFA2")
```


## Developmental version
To use the latest features of MOFA you can install the software from GitHub:
<!--
BiocManager::install(version='devel')
BiocManager::install("MOFA2")
-->
```r
devtools::install_github("bioFAM/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))
```

You will also need to make `mofapy2` available to R — see [below](#notes-on-the-connection-of-r-to-python).
If you'd like the development version of `mofapy2` to match, install it directly from GitHub:
```r
reticulate::py_require("git+https://github.com/bioFAM/mofapy2")
```

## Notes on the connection of R to Python

The connection between R and Python is done via [reticulate](https://rstudio.github.io/reticulate). There are three ways to make `mofapy2` available to the R package:

- **Let reticulate provision it** (simplest, requires reticulate >= 1.41). Declare the dependency before training, and reticulate sets up an isolated Python environment — downloading a suitable Python interpreter if none is available:
  ```r
  reticulate::py_require("mofapy2")
  MOFAobject <- run_mofa(MOFAobject)
  ```
- **Let [basilisk](https://bioconductor.org/packages/release/bioc/html/basilisk.html) handle it.** `run_mofa(MOFAobject, use_basilisk = TRUE)` uses a dedicated, version-pinned Python environment, created the first time you call it this way. It is also the option to choose if you use `MOFA2` alongside other R packages with conflicting Python dependencies, since basilisk runs Python in a separate process.
- **Use an existing Python installation.** Install the `mofapy2` package and its dependencies manually with `pip install mofapy2` (from the Unix terminal), then select that installation with [`reticulate::use_python()`](https://rstudio.github.io/reticulate/reference/use_python.html) (or `reticulate::use_condaenv()`) - this option needs the most configuration.


Note that the connection of R and python is the source of most problems when running MOFA, see our [troubleshooting](troubleshooting.html) page or reach us if you have issues.


## Using MOFA2 with older R versions

We recommend using R (>= 4.0) with `MOFA2`. If you want to use it with older R versions, you can install `MOFA2` as

```r
remotes::install_github("bioFAM/MOFA2", ref = "R36", build_opts = c("--no-resave-data --no-build-vignettes"))
```
Note, that this is only maintained intermittently and you will need to manually install the python package as described above and possibly configure the `reticulate` interface.


## Installation using Docker image

If you use Docker, you can build an image using the provided Dockerfile:

```
docker build -t mofa2 .
```

You will then be able to use R or Python from the container. 

```
docker run -ti --rm -v $DATA_DIRECTORY:/data mofa2 R
#                   ^
#                   |
#                    use `-v` to map a folder on your machine to a container directory
```

The command above will launch R with `MOFA2` and its dependencies installed while mounting `$DATA_DIRECTORY` to the container.

You can also pull [the pre-build image from dockerhub](https://hub.docker.com/r/gtca/mofa2).
