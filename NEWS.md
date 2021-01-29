---
layout: default
title: NEWS
---

**MOFA2 1.2.0 (latest)**:  
<!-- - Added contribution scores -->
- Improve interoperability with `Seurat` and `SingleCellExperiment`
- MOFA factors can be saved to a `Seurat` object using `add_mofa_factors_to_seurat`
- Automatically extract metadata from Seurat and SingleCellExperiment objects
- Improve memory usage and training time by optionally replacing float64 arrays with `float32` arrays (specified as `data_options` in the `prepare_mofa` step)
- `mofapy2` has been updated to version 0.6.0 and now it has its [own repository](https://github.com/bioFAM/mofapy2)


**MOFA2 1.1.7**:  
- Added [MEFISTO](https://www.biorxiv.org/content/10.1101/2020.11.03.366674v1) into MOFA2
- MOFA2 Package available via [Bioconductor](http://bioconductor.org/packages/release/bioc/html/MOFA2.html)
- Improving Python interface with [basilisk](http://www.bioconductor.org/packages/release/bioc/html/basilisk.html)
- Sample metadata can be incorporated to the `MOFA` object before and after training using the `samples_metadata` function