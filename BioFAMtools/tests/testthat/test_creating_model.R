context("Creating the model from different objects")
library(BioFAMtools)


test_that("a model can be created from a list of matrices", {
	m <- as.matrix(read.csv('matrix.csv'))
	expect_is(create_biofam(list("view1" = m)), "BioFAModel")
	expect_warning(create_biofam(list("view1" = m)))
	expect_error(create_biofam(m))
})

test_that("a model can be created from a list of sparse matrices", {
	skip_if_not_installed("Matrix")
	# Generate a sparse matrix
	m <- matrix(rnorm(100 * 5), ncol = 5) %*% t(matrix(rnorm(5 * 50), ncol = 5))
	m[sample(1:nrow(m), 100, replace = TRUE), sample(1:ncol(m), 100, replace = TRUE)] <- 0
	library(Matrix)
	m <- Matrix(m, sparse = TRUE)
	# Test if a sparse matrix can be imported to the BioFAModel
	expect_is(create_biofam(list("view1" = m)), "BioFAModel")
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

