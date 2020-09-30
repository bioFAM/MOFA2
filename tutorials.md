---
layout: default
title: Tutorials
---

<!-- ## Learning the basics -->

<!-- * [**Getting started**](https://github.com/bioFAM/MOFA2/blob/master/MOFA2/vignettes/getting_started.md): general overview and description of the method, for both R and Python users. -->

All .Rmd files can be found [here](https://github.com/bioFAM/MOFA2/tree/master/tutorials_extended)

All jupyter notebook files can be found ...

## Getting started using R

* [**Training a model in R**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/getting_started_R.html): using simple simulated data  

* [**Downstream analysis in R**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/downstream_analysis.html): using simple simulated data  

* [**Gene set enrichment analysis**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/GSEA.html): demonstrates how to do gene set enrichment analysis in R.  

* [**Demonstration of the stochastic inference algorithm**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/stochastic_inference.html): this is only useful for very large data sets and when having access to GPUs.


## Getting started using Python

* [**Training a model in Python**](https://github.com/bioFAM/MOFA2/blob/master/mofapy2/notebooks/getting_started_python.ipynb): a jupyter notebook demonstrating how to train a MOFA model using simple simulated data  

* **Downstream analysis in Python**: in preparation...


## Case examples using real data (in R)

* [**(authors' favourite) Analysis of chronic lymphocytic leukaemia cohort for personalised medicine**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/CLL.html): a bulk multi-omics data set. [Figure 2 and 3 of the MOFA v1 paper](https://www.embopress.org/doi/full/10.15252/msb.20178124#msb178124-fig-0002).  

* [**(authors' favourite) Integrative analysis of the Chromium Single Cell Multiome ATAC + Gene Expression assay**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/10x_scRNA_scATAC.html): this is the result of a collaboration between the MOFA team and the 10x Genomics R&D team to provide a downstream analysis pipeline for the new RNA+ATAC multi-modal technology.  

* [**Analysis of a time course scRNA-seq data set using the multi-group framework**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/scRNA_gastrulation.html): [Figure 2 of the MOFA+ paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1#Fig2). Demonstrates the multi-group functionality and how to train a MOFA model from a Seurat object.  

* [**Integration of scNMT-seq data  (single-cell multi-omics)**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/scNMT_gastrulation.html): [Figure 4 of the MOFA+ paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1#Fig4). Demonstrates the simultaneous multi-view and multi-group functionality using the [scNMT-seq mouse gastrulation atlas](https://www.nature.com/articles/s41586-019-1825-8).  

* [**Integration of SNARE-seq data (single-cell multi-omics)**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/SNARE_seq.html). Demonstrates how MOFA can be used for the analysis of paired scRNA+scATAC data (from the same cell) using a multi-modal Seurat object. This data set is very noisy and the results are not fantastic, we suggest you have a look at the Chromium Single Cell Multiome ATAC + Gene Expression vignette instead.  

* [**Analysis of multi-modal microbiome data**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/microbiome_vignette.html): we demonstrate how to systematically integrate viral, fungal and bacterial sequence data.
*
<!-- * [**Robustness analysis and model selection**](https://raw.githack.com/bioFAM/MOFA2/master/MOFA2/vignettes/model_selection.html) -->

<!-- * [**Analysis of single-cell DNA methylation data (in R)**](https://github.com/bioFAM/MOFA2/blob/master/MOFA2/vignettes/scMethylation_cortex.html): Figure 3 of the paper, in preparation... -->

* **Analysis of CITE-seq data**: still in preparation, reach us if you have questions... 
 

## Case examples (in Python)

In preparation...


## Template scripts

We provide two template scripts to train your model in [R](https://github.com/bioFAM/MOFA2/blob/master/template_script.R) and [Python](https://github.com/bioFAM/MOFA2/blob/master/template_script.py). Please see the tutorials for more details regarding the input data format, etc.

## Documentation

* [R package](https://github.com/bioFAM/MOFA2/blob/master/documentation.pdf)  
* [Python package]()  
