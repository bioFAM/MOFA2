---
layout: default
title: Tutorials
---

## Precorded talks

* [**MOFA overview**](https://www.youtube.com/watch?v=_BfHeZ0s2i0): precorded talk for the VIB workshop (Belgium, 2021), includes the model overview, intuition and a brief discussion of the [CLL application](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/CLL.html).

* [**Overview of single-cell multi-omics data integration**](https://www.youtube.com/watch?v=4Nt4oz0cfIk): precorded talk for a webinar, includes brief discussion on CITE-seq and 10x Multiome applications.

## Getting started using R

* [**Training a MOFA model in R**](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/getting_started_R.html): using simple simulated data  

* [**Downstream analysis in R**](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/downstream_analysis.html): using simple simulated data  

* [**Gene set enrichment analysis**](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/GSEA.html): demonstrates how to do gene set enrichment analysis in R.  

* [**Demonstration of the stochastic inference algorithm**](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/stochastic_inference.html): this is only useful for very large data sets and when having access to GPUs.


## Case examples using real data (in R)

* [**(authors' favourite) Analysis of chronic lymphocytic leukaemia cohort for personalised medicine**](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/CLL.html): a bulk multi-omics data set. Figure 2 and 3 of the [MOFA paper (https://www.embopress.org/doi/full/10.15252/msb.20178124#msb178124-fig-0002).  

* [**Integrative analysis of the Chromium Single Cell Multiome ATAC + Gene Expression assay**](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/10x_scRNA_scATAC.html): this is the result of a collaboration between the MOFA team and the 10x Genomics R&D team to provide a downstream analysis pipeline for the 10x Multiome kit.  

* [**Analysis of multi-modal microbiome data**](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/microbiome_vignette.html): we demonstrate how to systematically integrate viral, fungal and bacterial sequence data. Manuscript published in [mSystems](https://msystems.asm.org/content/6/2/e01148-20)

* [**Analysis of a time course scRNA-seq data set using the multi-group framework**](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/scRNA_gastrulation.html): Figure 2 of the [MOFA+ paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1#Fig2). Demonstrates the multi-group functionality and how to train a MOFA model from a Seurat object.  

* [**Integration of scNMT-seq data  (single-cell multi-omics)**](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/scNMT_gastrulation.html): Figure 4 of the [MOFA+ paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1#Fig4). Demonstrates the simultaneous multi-view and multi-group functionality using the [multi-modal mouse gastrulation atlas](https://www.nature.com/articles/s41586-019-1825-8).  

* [**Integration of SNARE-seq data (single-cell multi-omics)**](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/SNARE_seq.html). Demonstrates how MOFA can be used for the analysis of paired scRNA+scATAC data (from the same cell) from a Seurat object. This data set is very noisy and the results are not fantastic, we suggest you have a look at the [10x Multiome vignette](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/10x_scRNA_scATAC.html) instead.  

<!-- * [**Robustness analysis and model selection**](https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/model_selection.html) -->

<!-- * [**Analysis of single-cell DNA methylation data (in R)**](https://github.com/bioFAM/MOFA2/blob/master/MOFA2/vignettes/scMethylation_cortex.html): Figure 3 of the paper, in preparation... -->

<!-- * **Analysis of CITE-seq data**: still in preparation, reach us if you have questions...  -->
 
All .Rmd files for the tutorials can be found [here](https://github.com/bioFAM/MOFA2_tutorials/tree/master/R_tutorials)


## Getting started using Python

MOFA is interfaced through [muon](https://github.com/scverse/muon):

* [**Training MOFA with muon, CLL**](https://muon-tutorials.readthedocs.io/en/latest/CLL.html): the multimodal CLL dataset in the [MuData](https://github.com/scverse/mudata) format.

* [**Integration of multimodal single-cell data**](https://muon-tutorials.readthedocs.io/en/latest/single-cell-rna-atac/pbmc10k/3-Multimodal-Omics-Data-Integration.html): shows how [muon](https://github.com/scverse/muon) can be used to train and explore the model together with the [mofax](https://github.com/bioFAM/mofax) package.

* [**Training a model with mofapy2**](https://github.com/bioFAM/mofapy2/blob/master/mofapy2/notebooks/getting_started_python.ipynb): a Jupyter notebook demonstrating how to train a MOFA model using simple simulated data. It uses stand-alone [mofapy2](https://github.com/bioFAM/mofapy2) to train the model and [mofax](https://github.com/bioFAM/mofax) for downstream analysis.

* [**Training a model on AnnData**](https://github.com/bioFAM/mofax/blob/master/notebooks/training_pbmc10k.ipynb): demonstrates how to train a MOFA model on scRNA-seq data stored in [AnnData format](https://github.com/scverse/anndata).

* [**Downstream analysis in Python**](https://github.com/bioFAM/mofax/blob/master/notebooks/getting_started_pbmc10k.ipynb): demonstrates how to perform the downstream analysis of a MOFA model trained on scRNA-seq data, using [mofax](https://github.com/bioFAM/mofax).


## Template scripts

We provide template scripts to train your model:

* [R](https://github.com/bioFAM/MOFA2/blob/master/inst/scripts/template_script.R)
* [Python, using a list of data matrices](https://github.com/bioFAM/MOFA2/blob/master/inst/scripts/template_script_matrix.py). 
* [Python, using the long data frame format](https://github.com/bioFAM/MOFA2/blob/master/inst/scripts/template_script_dataframe.py)
