---
layout: default
title: Troubleshooting
---


### I get the following error when running `run_mofa` 

```
AttributeError: 'module' object has no attribute 'core.entry_point

Error in py_module_import(module, convert = convert) :
 ModuleNotFoundError: No module named 'mofapy2'
```
First thing: restart R and try again. If the error still holds, this means that either you did not install the mofapy2 Python package (see instructions above), or you have multiple Python installations and R is not detecting the correct one where mofapy2 is installed. You need to find out the right Python interpreter (which usually will be the one you get when running `which python` in the terminal) and specify the following at the beginning of your R script:
```
library(reticulate)
use_python("YOUR_PYTHON_PATH", required=TRUE)
```
You can also use `use_conda` instead of `use_python` if you work with conda environments. For details read more about the [reticulate](https://rstudio.github.io/reticulate/) package.


### I get the following error when installing the R package

```
ERROR: dependencies 'XXX', 'YYY' are not available for package 'MOFA2'
```
You probably tried to install them using `install.packages()`. These packages should be installed from Bioconductor.
