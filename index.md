---
layout: default
title: MOFA
---

MOFA is a factor analysis model that provides a **general framework for the integration of multi-omic data sets** in an unsupervised fashion.  
Intuitively, MOFA can be viewed as a versatile and statistically rigorous generalization of principal component analysis to multi-omics data. Given several data matrices with measurements of multiple -omics data types on the same or on overlapping sets of samples, MOFA infers an **interpretable low-dimensional representation in terms of a few latent factors**. These learnt factors represent the driving sources of variation across data modalities, thus facilitating the identification of cellular states or disease subgroups. 

<p align="center"> 
<img src="images/mofa_overview.png" style="width: 100%; height: 100%"/>
</p>

For more details you can read our two papers: 
- [MOFA v1, published in in Molecular Systems Biology](http://msb.embopress.org/cgi/doi/10.15252/msb.20178124)  
- [MOFA v2, published in Genome Biology](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1)

## What changed from MOFA to MOFA2?

In MOFA2 we added the following improvements:  

* **Multi-group functionality**: intuitively, this functionality breaks the assumption of independent samples and allows for are predefined groups of samples (i.e. different conditions, batches, cohorts, etc.). Importantly, the model is not focused on capturing the differential changes between the groups (you will find no factors that *separate* the groups). The aim of the multi-group framework is to discover which sources of variability are shared between the different groups and which ones are exclusive to a single group.

* **Improved downstream visualisations**: see the documentation in the Tutorials section for a comprehensive list of available functions.

* **No need for model selection**: In MOFA v1 we used random parameter initialisation, which led to (slightly) different solutions depending on the initial conditions. In MOFA v2 we initialise the factors using Principal Component Analysis on the concatenated data set, and the weights are initialised to zero. If using standard variational inference (not stochastic) this removes the randomness in the training algorithm.

* **Speed**: the training procedure is now 2-3x faster in standard CPUs.

* **GPU support**: the training procedure can be massively accelerated using GPUs. For this you have to install and configure the [CuPy package](https://cupy.chainer.org).
