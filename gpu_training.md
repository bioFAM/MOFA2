---
layout: default
title: MOFA
---

## Overview

MOFA training can be massively speed up by using GPUs. We have implemented all computations using [CuPy](https://cupy.dev/), an open-source array library for GPU-accelerated computing with Python.

GPU training can be accessed from R via the training options `gpu_mode` and `gpu_device`. This needs a Python environment with CuPy, which the basilisk environment does not provide — so run with `use_basilisk = FALSE`. The simplest setup is to let reticulate provision it, choosing the CuPy build matching your CUDA version:

```r
reticulate::py_require(c("mofapy2", "cupy-cuda12x"))  # or cupy-cuda13x
```

## Installation

To use GPU need to make sure that you have a functional NVIDIA GPU with the right version of CUDA and CuPy installed. Please check the [CuPy installation instructions](https://docs.cupy.dev/en/stable/install.html)

### Docker file

[Frederik Ziebell](https://github.com/frederikziebell) kindly prepared a Dockerfile, available [here](https://github.com/bioFAM/MOFA2/issues/95). If you try it please leave your feedback in the github issue.

### Conda environment

First, create a file called `mofa_conda_env_gpu.yml` with the following content 
```
name: mofa_env_gpu
channels:
  - conda-forge
dependencies:
  - cupy
  - cudatoolkit=11.0 # note that you might have to change this depending on your desired cuda version
  - dtw-python
  - pip
  - pip:
    - mofapy2==0.6.7
```
and then run:
```
conda env create --file mofa_conda_env_gpu.yml
```

## Basic example

We provide some simple scripts to test whether your GPU set up is working:
- [Script to test basic linear alebra operations using CuPy](https://github.com/bioFAM/mofapy2/blob/master/mofapy2/notebooks/test_cupy.py)
- [Script to train a MOFA model with a toy data set](https://github.com/bioFAM/mofapy2/blob/master/mofapy2/notebooks/run_mofa_cpu_vs_gpu.py)
