# Multi-Omics Factor Analysis v2 (MOFA+)

MOFA is a factor analysis model that provides a **general framework for the integration of multi-omic data sets** in an unsupervised fashion.  
Intuitively, MOFA can be viewed as a versatile and statistically rigorous generalization of principal component analysis (PCA) to multi-omics data. Given several data matrices with measurements of multiple ‘omics data types on the same or on overlapping sets of samples, MOFA infers an **interpretable low-dimensional data representation in terms of (hidden) factors**. These learnt factors represent the driving sources of variation across data modalities, thus facilitating the identification of cellular states or disease subgroups.  

In MOFA v2 (MOFA+) we added the following improvements:
* Fast Stochastic variational inference framework amenable to GPU computations: enables inference with very large data sets
* Multi-group functionality: intuitively, this breaks the assumption of independent samples and allows inference across multiple groups, where groups are predefined sets of samples (i.e. different conditions, batches, cohorts, etc.).


For more details you can read our papers: 
- MOFA v1: http://msb.embopress.org/cgi/doi/10.15252/msb.20178124
- MOFA v2: XXX

<p align="center"> 
<img src="images/logo.png" style="width: 50%; height: 50%"/>​
</p>


## Installation

The core of MOFA is implemented in Python. However, the whole procedure can be run with R and we provide the downstream analysis functions only in R.

### Python dependencies 

Python dependencies can be installed using pip (from the Unix terminal)

```r
pip install mofapy2
```

Alternatively, they can be installed from R itself using the reticulate package:

```r
library(reticulate)
py_install("mofapy2", envname = "r-reticulate", method="auto")
```

### MOFA2 R package

Can be installed using R:

```r
devtools::install_github("bioFAM/MOFA2", build_opts = c("--no-resave-data"))
```

--------------

### Using Docker image

You can build an image with `mofa2py` python library and `MOFA2` R package using the provided Dockerfile:

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

## Usage

### Step 1: Prepare the data

- Data processing
- Filtering
- Regressing out technical variation


If you work with single-cell data, MOFA+ comes with interfaces to build and train a model directly from objects commonly used for scRNA-seq data analysis, namely [AnnData](https://github.com/theislab/anndata) ([scanpy](https://github.com/theislab/scanpy)) in Python and [Seurat](https://github.com/satijalab/seurat) in R.

See the vignette XXXX and the documentation for details


### Step 2: Fitting the model

### Step 3: Downstream analysis
- Disentangling variance explained across views and groups
- Visualisation of factors
- Visualisation of loadings
- Transfer learning and imputation
	

Downstream analysis: disentangle the variability between omics

## Tutorials/Vignettes
We currently provide the following vignettes:

* **Data processing and creation of MOFA object**: bazzz
* **Integration of heterogeneous scRNA-seq data**: foo.
* **Integration of heterogeneous DNA methylation data**: bar.
* **Integration of single-cell multi-modal data:**: baz.
* **Transfer learning and imputation:**: baz.
* **Model selection and robustness with simulated data**: bazz



## Frequently asked questions

## Contact
The package is maintained by Ricard Argelaguet (ricard@ebi.ac.uk) and Danila Bredikhin (danila.bredikhin@embl.de ). Please, reach us for problems, comments or suggestions. You can also contact us via a Slack group where we provide quick and personalised help, [this is the link](https://join.slack.com/t/mofahelp/shared_invite/enQtMjcxNzM3OTE3NjcxLTkyZmE5YzNiMDc4OTkxYWExYWNlZTRhMWI2OWNkNzhmYmNlZjJiMjA4MjNiYjI2YTc4NjExNzU2ZTZiYzQyNjY).  


