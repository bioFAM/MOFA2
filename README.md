# bioFAM

Biology Factor Analysis Models (biofam) is a flexible framework for biological data analysis. Relying on a flexible graphical model, biofam enables users to model different assumptions for their data. Example usage is applying biofam to the datasets with multiple data types for the same samples or to multiple groups of samples or cells with same observed features (e.g. gene expression values).

## Getting Started

### Installing

To install biofam with `pip` run:

```
pip install git+https://github.com/bioFAM/biofam.git --user
```

The BioFAMtools R package can be installed with `devtools`:

```
devtools::install_github("bioFAM/biofam", subdir="BioFAMtools")
```

An alternative way is to clone the repository and install biofam library and BioFAMtools package from source:

```
git clone git@github.com:bioFAM/biofam.git
cd biofam

python setup.py install

R CMD INSTALL --build BioFAMtools
```

### Using Docker image

You can build an image with biofam python library and R package using the provided [Dockerfile](./Dockerfile):

```
docker build -t biofam .
```

## Usage

TODO: basic usage, studying factors, multi-omics tutorials.

### scRNA-seq data

bioFAM comes with interfaces to build and train a model directly from objects commonly used for scRNA-seq data analysis, namely [AnnData](https://github.com/theislab/anndata)([scanpy](https://github.com/theislab/scanpy)) in Python and [Seurat](https://github.com/satijalab/seurat) in R.

#### With scanpy

Use an AnnData object to build and train a MOFA model:

```{python}
mf.set_data_from_anndata(adata, "louvain")
```

For more information see this tutorial (TODO: tutorial on PBMC with scanpy).

#### With Seurat

```{r}
mf <- create_biofam(seurat_object, "louvain")
```

For more information see this tutorial (TODO: tutorial on PBMC with Seurat).

## Authors

In alphabetical order:

* Ricard Argelaguet
* Damien Arnol
* Danila Bredikhin
* Yonatan Deloro
* Britta Velten


[StatGenomics group on Twitter](https://twitter.com/statgenomics).


## License

TBA

