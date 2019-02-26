context("Creating the model from different objects")
library(BioFAMtools)


test_that("a model can be created from a list of matrices", {
	m <- as.matrix(read.csv('matrix.csv'))
	expect_is(create_biofam(list("view1" = m)), "BioFAModel")
	expect_warning(create_biofam(list("view1" = m)))
	expect_error(create_biofam(m))
})

test_that("a model can be created from a Seurat object", {
	skip_if_not_installed("Seurat")
	library(Seurat)
	m <- readMM(url('https://github.com/satijalab/seurat/blob/master/tests/testdata/matrix.mtx?raw=true'))
	genes <- read.delim(url('https://github.com/satijalab/seurat/blob/master/tests/testdata/genes.tsv?raw=true'), sep='\t', header=FALSE)[,1]
	cells <- read.delim(url('https://github.com/satijalab/seurat/blob/master/tests/testdata/barcodes.tsv?raw=true'), sep='\t', header=FALSE)[,1]
	colnames(m) <- cells
	rownames(m) <- genes
	srt <- Seurat::CreateSeuratObject(m)
	expect_is(create_biofam(srt), "BioFAModel")
})

