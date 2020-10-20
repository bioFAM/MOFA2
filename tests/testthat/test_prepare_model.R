context("Prepare the model from different objects")
library(MOFA2)


test_that("a MOFA model can be prepared from a list of matrices", {
	m <- as.matrix(read.csv('matrix.csv'))
	# Set feature names
	rownames(m) <- paste("feature_", seq_len(nrow(m)), paste = "", sep = "")
	# Set sample names
	colnames(m) <- paste("sample_", seq_len(ncol(m)), paste = "", sep = "")
	factor_model <- create_mofa(list("view1" = m))
	expect_is(prepare_mofa(factor_model), "MOFA")
})

test_that("a model can be created from a list of sparse matrices", {
	skip_if_not_installed("Matrix")

	# Generate a sparse matrix
	m <- matrix(rnorm(100 * 5), ncol = 5) %*% t(matrix(rnorm(5 * 50), ncol = 5))
	m[sample(1:nrow(m), 100, replace = TRUE), sample(1:ncol(m), 100, replace = TRUE)] <- 0
	library(Matrix)
	m <- Matrix(m, sparse = TRUE)

	# Set feature names
	rownames(m) <- paste("feature_", seq_len(nrow(m)), paste = "", sep = "")
	# Set sample names
	colnames(m) <- paste("sample_", seq_len(ncol(m)), paste = "", sep = "")
	# Initialise a model
	factor_model <- create_mofa(list("view1" = m))
	
	# Test if a sparse matrix can be used to prepare the MOFA model for training
	expect_is(prepare_mofa(factor_model), "MOFA")
})

test_that("a model can be created from a Seurat object", {
	skip_if_not_installed("Seurat")

	# Create a Seurat object from demo data
	library(Seurat)
	library(Matrix)
	m <- readMM('matrix.mtx')
	genes <- read.delim('genes.tsv', sep='\t', header=FALSE, stringsAsFactors=FALSE)[,2]
	cells <- read.delim('barcodes.tsv', sep='\t', header=FALSE, stringsAsFactors=FALSE)[,1]
	colnames(m) <- cells
	rownames(m) <- genes
	srt <- Seurat::CreateSeuratObject(m)

	# Prepare a model before training
	factor_model <- create_mofa(srt, features = genes)

	# Test if a Seurat object can be used to prepare the MOFA model for training
	expect_is(prepare_mofa(factor_model), "MOFA")
})

