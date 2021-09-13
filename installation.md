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

This uses [basilisk](https://bioconductor.org/packages/release/bioc/html/basilisk.html) to automatically set up a Python environment and install all required dependencies. 


## Developmental version
To use the latest features of MOFA (see [NEWS](https://biofam.github.io/MOFA2/news.html)) you can install the latest version from GitHub:
<!--
BiocManager::install(version='devel')
BiocManager::install("MOFA2")
-->
```r
devtools::install_github("bioFAM/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))
```

If you do so, you have to manually install the Python dependencies using pip (from the Unix terminal). Importantly, this has to be done before the R installation. 
```r
pip install mofapy2
```
In addition, it is very likely that you will have to connect R to Python manually using the `reticulate` interface (see paragraph below).

## Notes on the connection of R to Python

The connection between R and Python is dona via [reticulate](ttps://rstudio.github.io/reticulate). Latest version of `MOFA2` use [basilisk](https://bioconductor.org/packages/release/bioc/html/basilisk.html) to automatically set up a Python environment and install all required dependencies. Alternatively, you can install the python pacakge `mofapy2` manually as described above and specify to use this installation when running MOFA. Note that this sometimes this needs [configuration](https://rstudio.github.io/reticulate/reference/use_python.html) and it is the source of most problems in the `MOFA2` R package, specially when you have multiple versions of Python installed. See our FAQ section or reach us if you have issues.

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
