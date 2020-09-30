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
- [MOFA v2, published in Genome Biology](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1])

<!-- ## Web server
We provide a [Shiny-based web server](http://www.ebi.ac.uk/shiny/mofa/) to interactively explore MOFA models. Note that the web server only provides basic functionalities. For a comprehensive analysis please use the MOFA2 R package.  
You can also download the latest version from the corresponding [github repository](https://github.com/gtca/mofaplus-shiny/) -->


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


