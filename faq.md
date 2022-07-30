---
layout: default
title: FAQ
---

## FAQ on the data processing

**(Q) How do I normalise the data?**  
Proper normalisation is critical for the model to work. First, one needs to remove library size effects. For count-based data such as RNA-seq or ATAC-seq we recommend size factor normalisation + variance stabilisation. If this is not done correctly, the model will learn a very strong Factor 1 that will capture differences in the _total_ expression per sample, and more subtle sources of variation will be downweighted.  

**(Q) Should I do any filtering to the input data?**  
It is strongly recommended that you filter highly variable features (HVGs) per assay. When doing multi-group inference, you have to regress out the group effect before selecting HVGs.

**(Q) How many samples do I need?**  
Factor Analysis models are only useful with large sample sizes, at least more than 15.

**(Q) Should I remove undesired sources of variability (i.e. batch effects) before fitting the model?**  
Yes. If you have clear technical factors, we strongly encourage to regress it out a priori using for example a linear model like [limma](https://bioconductor.org/packages/release/bioc/html/limma.html). This is very important for the feature selection step. If you do not remove technical variability then MOFA will "focus" on capturing this huge variability driven by the technical factors, and smaller sources of (biological) variability could be missed.

**(Q) My data sets have different dimensionalities, does this matter?**  
Yes. Bigger data modalities will tend to be overrepresented in the factors. It is good practice to filter uninformative features (based for example on a minimum variance threshold) in order to have the different views within the same order of magnitudes. If this is unavoidable, take into account that the model has the risk of missing (small) sources of variation in the small data set.

**(Q) Does MOFA handle missing values?**  
Yes. It simply ignores them from the likelihood, there is no hidden imputation step. Matrix factorisation models are known to be very robust to the presence of missing values!


## FAQ on the transition of the implementation from `MOFA` to `MOFA2`

**(Q) What happened to `MOFA` and what has changed in `MOFA2`?**  
The R package `MOFA` and python package `mofapy` have been superseded by their new versions `MOFA2` and python package `mofapy2`. The new implementation inherits all the features of the old version and remains applicable for both bulk and single cell multi-omic data but provides additional functionalities (multi-group framework, MEFISTO) and is much faster. Have a look at our [News site](https://biofam.github.io/MOFA2/NEWS.html) for most recent developments and a list of changes in `MOFA2` compared to `MOFA`.


## FAQ on the downstream analysis

**(Q) How do I interpret the factors?**  
The MOFA factors capture the global sources of variability in the data. Mathematically, each factor ordinates cells along a one-dimensional axis centered at zero. The value per se is not interpretable, only the relative positioning of samples is important. Samples with different signs manifest opposite "effects" along the inferred axis of variation, with higher absolute value indicating a stronger effect. Note that the interpretation of factors is analogous to the interpretation of the principal components in PCA.

**(Q) How do I interpret the weights?**  
The weights provide a score for how strong each feature relates to each factor, hence allowing a biological interpretation of the latent factors. Features with no as- sociation with the factor have values close to zero, while genes with strong association with the factor have large absolute values. The sign of the weight indicates the direction of the effect: a positive weight indicates that the feature has higher levels in the cells with positive factor values, and vice versa.

<!-- **How can I do Gene Set Enrichment Analysis?**  
This is explained in the [GSEA vignette](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/GSEA.html) -->

<!-- **How can I assess the robustness of factors?** 
A procedure that can be applied to evaluate the robustness of factors is to downsample the number of samples and/or the number of features and inspect if the factors are consistently found. However, keep in mind that there could be cases where the full data set is required to detect small yet important sources of variation. Hence, lack of robustness under downsampling does not necessarily imply that a factor is not biologically meaningful. -->

## FAQ on the software

**(Q) Can I do MOFA only with Python?**  
You can use Python to train the model, see [this notebook](https://github.com/bioFAM/MOFA2/blob/master/mofapy2/notebooks/getting_started_python.ipynb) and [this template script](https://github.com/bioFAM/MOFA2/blob/master/template_script.py). However, we currently do not provide downstream analysis functions in Python (it is in our to-do list). For now we strongly recommend that you use the `MOFA2` R package for the analysis.

**(Q) Can I speed up the training procedure using CPU parallel processing?**  
MOFA uses [numpy](https://numpy.org/) for the mathematical operations. This library can be massively optimised by linking it to OpenBLAS or the Intel MKL libraries, which take advantage of multiple cores and multithreading. 

You can check which libraries you have linked to numpy using 
```
python -c "import numpy; print(numpy.show_config())"
```
Note that if you are using anaconda, numpy is automatically linked to MKL (see https://docs.anaconda.com/mkl-optimizations/).  
The next step is to define the environmental variables. For MKL you need to set `MKL_NUM_THREADS=N` where N is the number of cores.

**How can I use GPUs to speed up training?**  
The Python core of MOFA can take advantage of NVIDIA GPUs to massively speed up training. For this you have to install and configure the [CuPy package](https://cupy.chainer.org), which is an open-source matrix library accelerated with NVIDIA CUDA. 


## FAQ on the multi-group functionality

**(Q) How does the multi-group inference work?**  
The aim of the multi-group framework is not to capture differential changes in *mean* levels between the groups (as for example when doing differential RNA expression). The goal is to compare the sources of variability that drive each group. If your aim is to find a factor that "separates" the groups, you _DO NOT_ want to use the multi-group framework. In this setting, the features are centered per group before fitting the model.

**(Q) How do I define groups?**  
The selection of groups is hypothesis-driven, and typically motivated by the experimental design. There is no "right" or "wrong" definition of groups, but some definitions will be more useful than others. However, the user always needs to keep in mind that the aim of the multi-group framework is not to capture differential changes between the groups. The aim is to find out which sources of variability are shared between the different groups and which ones are exclusive to a single group. To achieve this, the group effect is regressed out from the data before fitting the model.

<!-- **(Q) How do I assess the quality/robustness of groups?**  
A quick approach to assess the validity of groups is to inspect the resulting variance explained plot. If the groups are too granular, the model will not recover significant amounts of variation. If the groups are not "interesting", this can result in a lack of "structure" in the variance explained plot (i.e. all factors being shared across all groups). 
More computationally intensive approaches can be used to assess the robustness of groups, including cross-validation and downsampling or bootstrapping samples within groups. -->


## FAQ on the model options

**(Q) How many factors should I learn?**  
The optimal number of factors depends on the aim of the analysis, the dimensions of the assays, the complexity of the data, there is no simple answer. In general, if the aim is to identify the major sources of biological variation one would typically consider the top 10 factors or so. In other tasks, such as imputation of missing values, even small sources of variation can be important and hence models should be trained with a large number of factors. 

**(Q) Can I include known covariates in the model?**  
We extensively tested this functionality and it was not yielding good results. The reason is that covariates are usually discrete labels that do not reflect the underlying molecular biology. For example, if you introduce age as a covariate, but the actual age is different from the molecular age, the model will simply learn a new factor that corresponds to this _latent_ molecular age, and it will drop the covariate from the model.  
We recommend that you learn the factors in a completely unsupervised manner and then relate them to the biological covariates a posteriori (see vignettes). If your covariate of interest is an important driver of variability, do not worry, MOFA will find it! 

<!-- **(5.4) The factors and weights have different values between runs. Is this expected?**  
This is normal and it happens because factor analysis models are rotation invariant. This means that you can rotate your factors and your weights and still find the same solution. This implies that the signs of the weight or the factors can NOT be compared across trials, only within a trial. -->

**(Q) What data modalities can MOFA cope with?**  
* Continuous data: modelled using a gaussian likelihood
* Binary data: modelled using a bernoulli likelihood
* Count data: using a poisson likelihood  

Importantly, the use of non-gaussian likelihoods require statistical approximations and are not as accurate as the gaussian likelihood. If your data can be safely transformed to match the gaussian likelihood assumptions, this is ALWAYS recommended. For example RNA-seq data is expected to be normalised and modelled with a gaussian distribution, do not input the counts directly.

<!-- **(Q) Do I need to do model selection?**  
Not anymore. In `MOFA` we did random parameter initialisation, which led to (slightly) different solutions depending on the initial conditions. In `MOFA2` we initialise the factors using Principal Component Analysis on the concatenated data set, and the weights are initialised to zero. If using standard variational inference (not stochastic) this removes the randomness in the training algorithms. -->
