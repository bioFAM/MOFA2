---
layout: default
title: NEWS
---

**`MOFA2 1.3.4` (latest)**:  
- Added a more flexible alignment option in MEFISTO to align distinct sets of groups instead of individual groups
- `mofapy2` has been updated to version 0.6.4


**`MOFA2 1.2.0`**:  
<!-- - Added contribution scores -->
- Improve interoperability with `Seurat` and `SingleCellExperiment`
- MOFA factors can be saved to a `Seurat` object using `add_mofa_factors_to_seurat`
- Automatically extract metadata from Seurat and SingleCellExperiment objects
- Improve memory usage and training time by optionally replacing float64 arrays with `float32` arrays (specified as `data_options` in the `prepare_mofa` step)
- `mofapy2` has been updated to version 0.6.0 and now it has its [own repository](https://github.com/bioFAM/mofapy2)


**`MOFA2 1.1.7`**:  
- Added [MEFISTO](https://www.biorxiv.org/content/10.1101/2020.11.03.366674v1) into `MOFA2`
- MOFA2 Package available via [Bioconductor](http://bioconductor.org/packages/release/bioc/html/MOFA2.html)
- Improving Python interface with [basilisk](http://www.bioconductor.org/packages/release/bioc/html/basilisk.html)
- Sample metadata can be incorporated to the `MOFAobject` before and after training using the `samples_metadata` function



**The `MOFA` package was deprecated and replaced by `MOFA2 1.0.0`**:  
The following new features are now available:

* **Multi-group functionality**: intuitively, this functionality breaks the assumption of independent samples and allows for are predefined groups of samples (i.e. different conditions, batches, cohorts, etc.). Importantly, the model is not focused on capturing the differential changes between the groups (you will find no factors that *separate* the groups). The aim of the multi-group framework is to discover which sources of variability are shared between the different groups and which ones are exclusive to a single group. Note that this functionality is optional.

* **Improved downstream visualisations**: see the documentation in the Tutorials section for a comprehensive list of available functions.

* **No need for model selection**: The old package used random parameter initialisation, which led to (slightly) different solutions depending on the initial conditions. In `MOFA2` we initialise the factors using Principal Component Analysis on the concatenated data set, and the weights are initialised to zero. If using standard variational inference (not stochastic) this removes the randomness in the training algorithm.

* **Speed**: the training procedure is now 2-3x faster in standard CPUs.

* **GPU support**: the training procedure can be massively accelerated using GPUs. For this you have to install and configure the [CuPy package](https://cupy.chainer.org).