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

