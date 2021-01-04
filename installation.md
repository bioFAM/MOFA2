---
layout: default
title: Installation
---

The core of MOFA is implemented in Python. However, the training and downstream analysis can also be run with R.

### Python dependencies (for Python users)

Python dependencies can be installed using pip (from the Unix terminal). This has to be done before opening the R terminal

```r
pip install mofapy2
```

### R package (for R users)

MOFA2 R package is available from Bioconductor ([link](http://www.bioconductor.org/packages/release/bioc/html/MOFA2.html))
and can be installed using:
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MOFA2")
```

#### Development version
To use the latest features of MOFA (see [NEWS](https://biofam.github.io/MOFA2/news.html), e.g.[MEFISTO](https://biofam.github.io/MOFA2/MEFISTO.html) for temporal and spatial data)) you can install the developmental version from Bioconductor or install the R package from GitHub:

```r
BiocManager::install(version='devel')
BiocManager::install("MOFA2")
```
or
```r
devtools::install_github("bioFAM/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))
```

#### Notes on the connection of R to Python
The connection between R and Python is dona via [reticulate](ttps://rstudio.github.io/reticulate). Latest version of MOFA2 use [basilisk](https://bioconductor.org/packages/release/bioc/html/basilisk.html) to automatically set up a Python environment and install all required dependcies. Alternatively, you can install the python pacakge `mofapy2` manually as described above and specify to use this installation when running mofa. Note that this sometimes this needs [configuration]((https://rstudio.github.io/reticulate/reference/use_python.html)) and it is the source of most problems in the MOFA R package, specially when you have multiple versions of Python installed. See our FAQ section or reach us if you have issues.

#### Using MOFA2 with older R versions

We recommend using R (>= 4.0) with MOFA2. If you want to use it with older R versions, you can install MOFA2 as

```r
remotes::install_github("bioFAM/MOFA2", ref = "R36", build_opts = c("--no-resave-data --no-build-vignettes"))
```
Note, that this is only maintained intermittently and you will need to manually install the python package as described above and possibly configure the `reticulate` interface.


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
