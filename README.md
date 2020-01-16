# Multi-Omics Factor Analysis v2 (MOFA+)

MOFA is a factor analysis model that provides a **general framework for the integration of multi-omic data sets** in an unsupervised fashion.  
Intuitively, MOFA can be viewed as a versatile and statistically rigorous generalization of principal component analysis (PCA) to multi-omics data. Given several data matrices with measurements of multiple -omics data types on the same or on overlapping sets of samples, MOFA infers an **interpretable low-dimensional data representation in terms of (hidden) factors**. These learnt factors represent the driving sources of variation across data modalities, thus facilitating the identification of cellular states or disease subgroups.  

In MOFA v2 (MOFA+) we added the following improvements:
* **Multi-group functionality**: intuitively, this breaks the assumption of independent samples and allows inference across multiple groups, where groups are predefined sets of samples (i.e. different conditions, batches, cohorts, etc.).
* **Fast inference** using a stochastic variational framework: this can be powered by GPUs: enabling inference with very large data sets.

For more details you can read our papers: 
- MOFA: http://msb.embopress.org/cgi/doi/10.15252/msb.20178124
- MOFA+: https://www.biorxiv.org/content/10.1101/837104v1

<p align="center"> 
<img src="images/figure1a_mofa2.png" style="width: 50%; height: 50%"/>​
</p>


## Manual installation

The core of MOFA is implemented in Python. However, the whole procedure can be run with R and we provide the downstream analysis functions only in R.

#### Python dependencies 

Python dependencies can be installed using pip (from the Unix terminal)

```r
pip install mofapy2
```

#### R package

MOFA2 R package can be installed using R:

```r
devtools::install_github("bioFAM/MOFA2/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))
```

--------------

### Installation using Docker image

You can build an image using the provided Dockerfile:

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

## Tutorials/Vignettes

* [**Getting started**](https://github.com/bioFAM/MOFA2/blob/master/MOFA2/vignettes/getting_started.md): general overview and description of the method.
* [**Training a model in R**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/getting_started_R.html)
* [**Training a model in Python (jupyter notebook)**](https://github.com/bioFAM/MOFA2/blob/master/mofapy2/notebooks/getting_started_python.ipynb)
* [**Downstream analysis (in R)**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/downstream_analysis.html)
* [**Analysis of a multi-group scRNA-seq data set (in R)**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/scRNA_gastrulation.html): Figure 2 of the paper.
* [**Analysis of single-cell DNA methylation data (in R)**](https://github.com/bioFAM/MOFA2/blob/master/MOFA2/vignettes/scMethylation_cortex.html): Figure 3 of the paper, in preparation...
* [**Integration of single-cell multi-modal data (in R)**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/scNMT_gastrulation.html): Figure 4 of the paper.
* **Using fast stochastic variational inference**: in preparation...
* **Analysis of CITE-seq data**: in preparation...
* **Analysis of chronic lymphocytic leukaemia cohort for personalised medicine**: in preparation...
* [**Robustness analysis and model selection**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/model_selection.html)

The data and the pre-trained models can be downloaded [here](ftp://ftp.ebi.ac.uk/pub/databases/mofa)


## Contact
E-mail: Ricard Argelaguet (ricard@ebi.ac.uk)
Slack: we have a Slack group where we provide quick and personalised help, [this is the link](https://join.slack.com/t/mofahelp/shared_invite/enQtMjcxNzM3OTE3NjcxLWNhZmM1MDRlMTZjZWRmYWJjMGFmMDkzNDBmMDhjYmJmMzdlYzU4Y2EzYTI1OGExNzM2MmUwMzJkZmVjNDkxNGI).


## Frequently asked questions (FAQ)

### (1) FAQ on the transition from MOFA to MOFA+

**(1.1) Can MOFA+ be applied to bulk data?**  
MOFA+ remains 100% applicable to bulk data. 

**(1.2) How does the multi-group inference work in MOFA+?**  
Groups must be defined a priori, they are not learnt by the model. There is  flexibility on how to define them, they usually correspond to different conditions, cohorts, time points, etc.  
Importantly, the model is not focused on capturing differential changes *between* groups. The aim is to find out which sources of variability are shared between the different groups and which ones are unique to a single group. To achieve this, the group effect is regressed out from the data before fitting the model.

**(1.3) Do I need to have multiple groups to use MOFA+?**  
No. Unless provided, MOFA+ assumes that you do not have multi-group structure in your data. In this case the model simplifies to MOFA v1 (but significantly faster).

**(1.4) Does MOFA+ inherit previous features from MOFA v1?**  
Yes, pretty much everything: handling of missing data, non-gaussian likelihoods and sparsity in the weights.


### (2) FAQ on the data processing

**(2.1) How do I normalise the data?**  
Proper normalisation of the data is critical for the model to work. First, one needs to remove library size effects. For count-based data such as RNA-seq or ATAC-seq we recommend size factor normalisation + variance stabilisation. For microarray DNA methylation data, make sure that samples have no differences in the average intensity. If this is not done correctly, the model will learn a very strong Factor 1 that will capture this variability, and more subtle sources of variation will be harder to identify.  

**(2.2) Should I do any filtering to the input data?**  
It is strongly recommended that you filter highly variable features (HVGs) per assay. When doing multi-group inference, you have to regress out the group effect before selecting HVGs.

**(2.3) How many samples do I need?**  
Factor Analysis models are only useful with large sample sizes (let's say more than ~15).

**(2.4) Should I remove undesired sources of variability (i.e. batch effects) before fitting the model?**  
Yes. If you have clear technical factors, we strongly encourage to regress it out a priori using a simple linear model. The reason for this is that the model will "focus" on the huge variability driven by the technical factors, and smaller sources of variability could be missed.
You can regress out known covaraites using the function `regressCovariates`. See the corresponding documentation and the CLL vignette for details.

**(2.5) My data sets have different dimensionalities, does this matter?**  
Yes. Bigger data modalities will tend to be overrepresented in the factors. It is good practice to filter features (based for example on variance) in order to have the different dimensionalities within the same order of magnitudes. If this is unavoidable, take into account that the model has the risk of missing (small) sources of variation unique to the small data set.

**(2.6) Does MOFA handle missing values?**  
Yes. It simply ignores them from the likelihood, there is no hidden imputation step. Matrix factorisation models are known to be very robust to the presence of missing values!


### (3) FAQ on the software

**(3.1) I get one of the following errors when running MOFA:**  
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

**(3.2) I get the following error when installing the R package:**  
```
ERROR: dependencies 'XXX', 'YYY' are not available for package 'MOFA2'
```
You probably tried to install them using `install.packages()`. These packages should be installed from Bioconductor.

**(3.3) I hate R, can I do MOFA only with Python?**  
XXX
<!-- You can use Python to train the model, see [this template script](https://github.com/bioFAM/MOFA2/blob/master/template_run.py). However, we currently do not provide downstream analysis functions in Python. We strongly recommend that you use our MOFA2 R package for this. -->


### (4) FAQ on the model options

**(4.1) How many factors should I learn?**  
Similar to other latent variable models, this is a hard question to answer. It depends on the data set and the aim of the analysis. If you want to get an overview on the major sources of variability then use a small number of factors (K<=10). If you want to capture small sources of variability, for example to do imputation or eQTL mapping, then go for a large number of factors (K>25).

**(4.2) Can MOFA automatically learn the number of factors?**  
Yes, but the user needs to specify a minimum value of % variance explained. Then, MOFA will actively remove factors (during training) that explain less than the specified amount of variance.
If you have no idea on what to expect, it is better to start with a fixed number of factors and set the % variance threshold to 0.

**(4.3) Can I include known covariates in the model?**  
We extensively tested this functionality and it was not yielding good results. The reason is that covariates are usually discrete labels that do not reflect the underlying molecular biology. For example, if you introduce age as a covariate, but the actual age is different from the “molecular age”, the model will simply learn a new factor that corresponds to this “latent” molecular age, and it will drop the covariate from the model.  
We recommend that you learn the factors in a completely unsupervised manner and then relate them to the biological covariates (see vignettes). If your covariate of interest is an important driver of variability, do not worry, MOFA will find it! 

**(4.4) The weights have different values between runs. Is this expected?**  
This is normal and it happens because of two reasons. The first one is that the model does not always converge to the same exact solution (see below in the FAQ), although different model instances should be pretty similar. The second reason is that factor analysis models are rotation invariant. This means that you can rotate your factors and your weights and still find the same solution. This implies that the signs of the weight or the factors can NOT be compared across trials, only within a trial.

**(4.5) What data modalities can MOFA cope with?**  
* Continuous data: modelled using a gaussian likelihood
* Binary data: modelled using a bernoulli likelihood
* Count data: using a poisson likelihood.
Importantly, the use of non-gaussian likelihoods require further approximations and are not as accurate as the gaussian likelihood. Hence, if your data can be safely transformed to match the gaussian likelihood assumptions, this is ALWAYS recommended. For example RNA-seq data is expected to be normalised and modelled with a gaussian distribution, do not input the counts directly.

**(4.6) The model does not converge smoothly, and it oscillates between positive and negative deltaELBO values**  
First, check that you are using the right likelihood model (see above). Second, make sure that you have no features or samples that are full of missing values. Third, check that you have no features with zero (or very little) variance. If the problem does not disappear, please contact us via mail or the Slack group.

**(4.7) Does MOFA always converge to the same solutions?**  
No, as occurs in most complex Bayesian models, they are not guaranteed to always converge to the same solution.
In practice, however, we observed that the solutions are highly consistent, particularly for the top factors. However, we recommend doing a robustness analysis. This is done by training multiple model instances and check the correlation of the factors across the different solutions See the function `compare_models()`.


### (5) FAQ on the downstream analysis

**(5.1) How can I do Gene Set Enrichment Analysis?**  
First, you need to create your binary gene set matrix where rows are feature sets and columns are features (genes). We have manually processed some of Reactome and MSigDB gene sets for mouse and human. Contact us if you would like to use the data.  
Then, you will have to choose a local statistic per feature (the loading, by default), a global statistic per pathway (average loading, by default), and a statistical test. The most trustworthy one is a permutation test with a long number of iterations, but this is slow and a fast parametric tests is also available. However, note that it tends to inflate the p-values due to the correlation structure between related genes (see for example [Gatti2010](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-11-574)).
