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
	skip_if_not_installed("SeuratObject")
	library(Seurat)
	library(Matrix)
	m <- readMM('matrix.mtx')
	genes <- read.delim('genes.tsv', sep='\t', header=FALSE, stringsAsFactors=FALSE)[,2]
	cells <- read.delim('barcodes.tsv', sep='\t', header=FALSE, stringsAsFactors=FALSE)[,1]
	colnames(m) <- cells
	rownames(m) <- genes
	srt <- SeuratObject::CreateSeuratObject(m)
	# only for testing purpose, should use scale.data
	expect_is(create_mofa(srt, features = genes, layer = "counts"), "MOFA")
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


test_that("a model can be created from a MultiAssayExperiment Object", {
	skip_if_not_installed("MultiAssayExperiment")
	library(MultiAssayExperiment)
	library(SummarizedExperiment)

	# Import and preprocessing of miniACC Data
	data(miniACC)
	miniACC <- intersectColumns(miniACC)
	mae_sub <- miniACC
	experiments(mae_sub) <- experiments(miniACC)[
		c("RNASeq2GeneNorm","RPPAArray","Mutations","miRNASeqGene")
	]
	#mae_sub <- miniACC[,,c("RNASeq2GeneNorm","RPPAArray","Mutations","miRNASeqGene")]
	
	
	# apply log1p to the RNASeq
	se <- mae_sub[["RNASeq2GeneNorm"]]
	assay(se, "log1p") <- log1p(assay(se, "exprs"))
	mae_sub[["RNASeq2GeneNorm"]] <- se

	# apply log1p to the miRNASeq
	se <- mae_sub[["miRNASeqGene"]]
	assay(se, "log1p") <- log1p(assay(se, "exprs"))
	mae_sub[["miRNASeqGene"]] <- se

	# create MOFA model
	model <- create_mofa(mae_sub,
                              assays = c("log1p", "exprs", "",'log1p'),
                              extract_metadata = TRUE)
	

	# do checks (??)

	expect_is(model, "MOFA")
})

test_that("a model can be created from a SingleCellExperiment Object", {
	skip_if_not_installed("SingleCellExperiment")
	skip_if_not_installed("MOFAdata")
	library(SingleCellExperiment)
	
	# Import CLL data
	utils::data("CLL_data", package = "MOFAdata")
	CLL_metadata <- data.table::fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt")
	IGHV <- CLL_metadata$IGHV 
	IGHV_filled <- ifelse(is.na(IGHV), 'NA', IGHV)
	CLL_metadata$IGHV_filled <- IGHV_filled

	# Create SCE Object
	sce <- SingleCellExperiment(assays = list(mRNA = CLL_data$mRNA), colData = CLL_metadata)
	altExps(sce) <- list(
	Drugs = SummarizedExperiment(list(expr = CLL_data$Drugs)),
	Methylation = SummarizedExperiment(list(expr = CLL_data$Methylation)),
	Mutations = SummarizedExperiment(list(expr = CLL_data$Mutations))
	)

	# create MOFA model
	MOFAobject <- create_mofa_from_SingleCellExperiment(
	sce, assay = "mRNA",
	alt_experiments = c("Drugs", "Methylation", "Mutations"), 
	alt_assays = c("expr","expr","expr"),
	groups = "IGHV_filled"
	)

	# do checks (??)
	expect_is(MOFAobject, "MOFA")
})