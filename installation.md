---
layout: default
title: Installation
---

The core of MOFA is implemented in Python. However, the training and downstream analysis can also be run with R

### Python dependencies (for both Python and R users)

Python dependencies can be installed using pip (from the Unix terminal). This has to be done before opening the R terminal

```r
pip install mofapy2
```

### R package (only for R users)

MOFA2 R package can be installed in R:

```r
devtools::install_github("bioFAM/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))
```

The connection between R and Python is dona via [reticulate](ttps://rstudio.github.io/reticulate). Sometimes this needs [configuration]((https://rstudio.github.io/reticulate/reference/use_python.html)) and it is the source of most problems in the MOFA R package, specially when you have multiple versions of Python installed. See our FAQ section or reach us if you have issues.

### Installation using Docker image

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

The command above will launch R with MOFA2 and its dependencies installed while mounting `$DATA_DIRECTORY` to the container.

You can also pull [the pre-build image from dockerhub](https://hub.docker.com/r/gtca/mofa2).
