context("Creating the model from different objects")
library(MOFA2)

test_that("a model can be created from a list of matrices", {
	m <- as.matrix(read.csv('matrix.csv'))
	expect_warning(create_mofa(list("view1" = m)))  # no feature names provided
	rownames(m) <- paste("feature", seq_len(nrow(m)), paste = "", sep = "")
	expect_is(create_mofa(list("view1" = m)), "MOFA")
	expect_error(create_mofa(m))
})

test_that("a model can be created from a list of sparse matrices", {
	skip_if_not_installed("Matrix")
	library(Matrix)
	# Generate a sparse matrix
	m <- matrix(rnorm(100 * 5), ncol = 5) %*% t(matrix(rnorm(5 * 50), ncol = 5))
	m[sample(1:nrow(m), 100, replace = TRUE), sample(1:ncol(m), 100, replace = TRUE)] <- 0
	m <- Matrix(m, sparse = TRUE)
	# Set feature names
	rownames(m) <- paste("feature_", seq_len(nrow(m)), paste = "", sep = "")
	# Set sample names
	colnames(m) <- paste("sample_", seq_len(ncol(m)), paste = "", sep = "")
	# Test if a sparse matrix can be imported to the MOFA
	expect_is(create_mofa(list("view1" = m)), "MOFA")
})

test_that("a model can be created from a Seurat object", {
	skip_if_not_installed("Seurat")
	library(Seurat)
	library(Matrix)
	m <- readMM('matrix.mtx')
	genes <- read.delim('genes.tsv', sep='\t', header=FALSE, stringsAsFactors=FALSE)[,2]
	cells <- read.delim('barcodes.tsv', sep='\t', header=FALSE, stringsAsFactors=FALSE)[,1]
	colnames(m) <- cells
	rownames(m) <- genes
	srt <- Seurat::CreateSeuratObject(m)
	expect_is(create_mofa(srt, features = genes), "MOFA")
})

test_that("a list of matrices per view is split correctly into a nested list of matrices according to samples groups", {
	n_groups <- 3
	# Create view 1
	m <- as.matrix(read.csv('matrix.csv'))
	rownames(m) <- paste("feature", seq_len(nrow(m)), paste = "", sep = "")
	colnames(m) <- paste("sample", seq_len(ncol(m)), paste = "", sep = "")
	# Add second view
	m2 <- m[1:(nrow(m)/3),]
	rownames(m2) <- paste("view2", rownames(m2), sep = "_")
	# Define multiple groups
	samples_groups <- sample(x = paste0("group", 1:n_groups), replace = TRUE, size = ncol(m))
	# Split the data
	data_split <- .split_data_into_groups(list("view1" = m, "view2" = m2), samples_groups)
	# Check group assignments
	for (g in 1:n_groups) {
		g_name <- paste0("group", g)
		expect_equal(colnames(data_split[[1]][[g_name]]), colnames(m)[which(samples_groups == g_name)])
		expect_equal(colnames(data_split[[2]][[g_name]]), colnames(m)[which(samples_groups == g_name)])
	}
})


