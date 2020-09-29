## What is MOFA?

MOFA is a factor analysis model that provides a **general framework for the integration of multi-omic data sets** in an unsupervised fashion.  
Intuitively, MOFA can be viewed as a versatile and statistically rigorous generalization of principal component analysis to multi-omics data. Given several data matrices with measurements of multiple -omics data types on the same or on overlapping sets of samples, MOFA infers an **interpretable low-dimensional representation in terms of a few latent factors**. These learnt factors represent the driving sources of variation across data modalities, thus facilitating the identification of cellular states or disease subgroups. 

For more details you can read our two papers: 
- MOFA v1: http://msb.embopress.org/cgi/doi/10.15252/msb.20178124  
- MOFA+: http://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1  

<p align="center"> 
<img src="images/figure1a_mofa2.png" style="width: 100%; height: 100%"/>
</p>


## Installation

The core of MOFA is implemented in Python. However, the whole procedure can be run with R and we provide the downstream analysis functions only in R.

### (1) Python dependencies (for training the model)

Python dependencies can be installed using pip (from the Unix terminal). This has to be done before opening the R terminal

```r
pip install mofapy2
```

### (2) R package (for downstream analysis)

MOFA2 R package can be installed using R:

```r
devtools::install_github("bioFAM/MOFA2/MOFA2", build_opts = c("--no-resave-data --no-build-vignettes"))
```

After this, if you have multiple versions of Python installed you may have to configure [reticulate](https://rstudio.github.io/reticulate/reference/use_python.html) to [connect R to the right Python binary](https://github.com/bioFAM/MOFA2#4-faq-on-the-software)
--------------

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

## Tutorials/Vignettes

### Learning the basics

* [**Getting started**](https://github.com/bioFAM/MOFA2/blob/master/MOFA2/vignettes/getting_started.md): general overview and description of the method.
* [**Training a model in R**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/getting_started_R.html)
* [**Training a model in Python (jupyter notebook)**](https://github.com/bioFAM/MOFA2/blob/master/mofapy2/notebooks/getting_started_python.ipynb)
* [**Downstream analysis in R**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/downstream_analysis.html)
* [**Downstream analysis in Python**:](https://github.com/gtca/mofax)
* [**Gene set enrichment analysis**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/GSEA.html): demonstrates how to do gene set enrichment analysis.

### Case examples

* [**(authors' favourite) Analysis of chronic lymphocytic leukaemia cohort for personalised medicine**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/CLL.html): a bulk multi-omics data set. Figure 2 and 3 of the MOFA v1 paper.
* [**(authors' favourite) Integrative analysis of the Chromium Single Cell Multiome ATAC + Gene Expression assay**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/10x_scRNA_scATAC.html): the new multi-modal protocol released by 10x Genomics.
* [**Analysis of a time course scRNA-seq data set using the multi-group framework**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/scRNA_gastrulation.html): Figure 2 of the MOFA+ paper.
* [**Integration of scNMT-seq data**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/scNMT_gastrulation.html): Figure 4 of the MOFA+ paper.
* [**Integration of SNARE-seq data: scRNA-seq and scATAC-seq from the same cell)**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/SNARE_seq.html)
* **Analysis of CITE-seq data**: still in preparation, reach us if you have questions...
* [**Analysis of multi-modal microbiome data**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/microbiome_vignette.html)
<!-- * [**Robustness analysis and model selection**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/model_selection.html) -->
* [**Demonstration of the stochastic inference algorithm (for very large data sets)**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/stochastic_inference.html)
<!-- * [**Analysis of single-cell DNA methylation data (in R)**](https://github.com/bioFAM/MOFA2/blob/master/MOFA2/vignettes/scMethylation_cortex.html): Figure 3 of the paper, in preparation... -->

## Web server

We provide a [Shiny-based web server](http://www.ebi.ac.uk/shiny/mofa/) to interactively explore MOFA models. Note that the web server only provides basic functionalities. For a comprehensive analysis please use the MOFA2 R package.  
You can also download the latest version from the corresponding [github repository](https://github.com/gtca/mofaplus-shiny/)

## Support or Contact

- **Slack (recommended)**: we have a Slack group where we provide quick and personalised help, [this is the link](https://join.slack.com/t/mofahelp/shared_invite/enQtMjcxNzM3OTE3NjcxLWNhZmM1MDRlMTZjZWRmYWJjMGFmMDkzNDBmMDhjYmJmMzdlYzU4Y2EzYTI1OGExNzM2MmUwMzJkZmVjNDkxNGI).
- e-mail: Ricard Argelaguet (ricard@ebi.ac.uk)

## Citation

    @article{Argelaguet2018,
        author = {Argelaguet, R. and Velten, B. and Arnol, D. and Dietrich, S. and Zenz, T. and Marioni, J. C. and Buettner, F. and Huber, W. and Stegle, O.},
        title = {Multi-Omics Factor Analysis-a framework for unsupervised integration of multi-omics data sets},
        journal = {Mol Syst Biol},
        volume = {14},
        number = {6},
        pages = {e8124},
        ISSN = {1744-4292 (Electronic) 1744-4292 (Linking)},
        DOI = {10.15252/msb.20178124},
        url = {https://www.ncbi.nlm.nih.gov/pubmed/29925568},
        year = {2018},
        type = {Journal Article}
    }

    @article{Argelaguet2020,
        Author = {Argelaguet, Ricard and Arnol, Damien and Bredikhin, Danila and Deloro, Yonatan and Velten, Britta and Marioni, John C. and Stegle, Oliver},
        Title = {MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data},
        Journal = {Genome Biology},
        Number = {1},
        Pages = {111},
        Volume = {21},
        Year = {2020}
    }


