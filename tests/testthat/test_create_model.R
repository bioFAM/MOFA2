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
	m <- as(readMM('matrix.mtx'),'dgCMatrix')
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
	rownames(mae_sub[["RNASeq2GeneNorm"]]) <- paste0(rownames(mae_sub[["RNASeq2GeneNorm"]]), "_1")
	rownames(mae_sub[["RPPAArray"]]) <- paste0(rownames(mae_sub[["RPPAArray"]]), "_2")
	rownames(mae_sub[["Mutations"]]) <- paste0(rownames(mae_sub[["Mutations"]]), "_3")
	rownames(mae_sub[["miRNASeqGene"]]) <- paste0(rownames(mae_sub[["miRNASeqGene"]]), "_4")

	
	
	# apply log1p to the RNASeq
	se <- mae_sub[["RNASeq2GeneNorm"]]
	assay(se, "log1p") <- log1p(assay(se, "exprs"))
	mae_sub[["RNASeq2GeneNorm"]] <- se

	# apply log1p to the miRNASeq
	se <- mae_sub[["miRNASeqGene"]]
	assay(se, "log1p") <- log1p(assay(se, "exprs"))
	mae_sub[["miRNASeqGene"]] <- se

	# create MOFA model
	model <- create_mofa(
		mae_sub,
		assays = c("log1p", "exprs", "",'log1p'),
		extract_metadata = TRUE
	)
	# Warning: duplicated feature names
	

	# do checks
	# class check
	expect_is(model, "MOFA")

	#check dimensions
	expect_equal(get_dimensions(model)$M,3) # right number of views

	# right data matrix
	expect_equivalent(
		get_data(model, views = c("RPPAArray"))$RPPAArray$group1,
		assay(mae_sub[["RPPAArray"]])
	) 

	#check also with HintikkaXOData

})

test_that("a model can be created from a MultiAssayExperiment Object - HintikkaXOData data", {
	skip_if_not_installed("MultiAssayExperiment")
	library(MultiAssayExperiment)
	library(SummarizedExperiment)
	library(mia)

	#check also with HintikkaXOData
	data("HintikkaXOData")

	# Prepare data 
	mae <- HintikkaXOData
	altExp(mae[[1]], "asd") <- mae[[1]]

	MOFAobject <- create_mofa(
		mae,
		alt_experiments = list("asd", "main", "main"),
		assays = list("counts",NULL,"signals"),
	)

	# do checks
	# class check
	expect_is(MOFAobject, "MOFA")

	#check dimensions
	expect_equal(get_dimensions(MOFAobject)$M,2) # right number of views

	MOFAobject <- create_mofa(
		mae,
		experiments = c("microbiota","biomarkers"),
		alt_experiments = list("asd", "main"),
		assays = list("counts","signals"),
	)

	# do checks
	# class check
	expect_is(MOFAobject, "MOFA")

	#check dimensions
	expect_equal(get_dimensions(MOFAobject)$M,2) # right number of views
})

test_that("a model can be created from a SingleCellExperiment Object", {
	skip_if_not_installed("SingleCellExperiment")
	skip_if_not_installed("MOFAdata")
	skip_if_not_installed("data.table")
	library(SingleCellExperiment)
	library(data.table)
	
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
		sce,
		alt_experiments = c("Main","Drugs", "Methylation", "Mutations"), 
		assays = c("mRNA","expr","expr","expr"),
		groups = "IGHV_filled"
	)

	# do checks
	# class check
	expect_is(MOFAobject, "MOFA")
	# right feature metadata (sample metadata not implemented)
	expect_identical(colData(sce)$sample, CLL_metadata$sample)
	expect_identical(colData(sce)$treatedAfter, CLL_metadata$treatedAfter)
	expect_identical(colData(sce)$age, CLL_metadata$age)
	expect_identical(colData(sce)$IGHV_filled, CLL_metadata$IGHV_filled)

	# right data matrix - with groups
	expect_equivalent(
		get_data(MOFAobject, views = c("Main"),groups = c("0"))$Main$`0`,
		CLL_data$mRNA[,CLL_metadata$IGHV==0 & !is.na(CLL_metadata$IGHV)]
	) 
	expect_equivalent(
		get_data(MOFAobject, views = c("Main"),groups = c("1"))$Main$`1`,
		CLL_data$mRNA[,CLL_metadata$IGHV==1 & !is.na(CLL_metadata$IGHV)]
	)

	#check dimensions
	expect_equal(get_dimensions(MOFAobject)$G,3) # right number of groups

	# right data matrix - without groups
	# create MOFA model
	MOFAobject <- create_mofa_from_SingleCellExperiment(
		sce,
		alt_experiments = c("Main","Drugs", "Methylation", "Mutations"),
		assays = c("mRNA","expr","expr","expr")
	)
    expect_equivalent(
		get_data(MOFAobject, views = c("Main"))$Main$group1,
		CLL_data$mRNA
	)
})

test_that("mofa2 wrapper correctly initializes MOFAobject", {
	#check that manually setting the parameters creates the same model as the mofa2() wrapper

	# 1. get sample data
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
	rownames(mae_sub[["RNASeq2GeneNorm"]]) <- paste0(rownames(mae_sub[["RNASeq2GeneNorm"]]), "_1")
	rownames(mae_sub[["RPPAArray"]]) <- paste0(rownames(mae_sub[["RPPAArray"]]), "_2")
	rownames(mae_sub[["Mutations"]]) <- paste0(rownames(mae_sub[["Mutations"]]), "_3")
	rownames(mae_sub[["miRNASeqGene"]]) <- paste0(rownames(mae_sub[["miRNASeqGene"]]), "_4")

	# apply log1p to the RNASeq
	se <- mae_sub[["RNASeq2GeneNorm"]]
	assay(se, "log1p") <- log1p(assay(se, "exprs"))
	mae_sub[["RNASeq2GeneNorm"]] <- se

	# apply log1p to the miRNASeq
	se <- mae_sub[["miRNASeqGene"]]
	assay(se, "log1p") <- log1p(assay(se, "exprs"))
	mae_sub[["miRNASeqGene"]] <- se

	# 2. create MOFA instances
	MOFA_init <- create_mofa(
		mae_sub,
		assays = c("log1p", "exprs", "",'log1p'),
		extract_metadata = TRUE
	)

	data_opts <- get_default_data_options(MOFA_init)
	model_opts <- get_default_model_options(MOFA_init)
	model_opts$num_factors <- 4
	train_opts <- get_default_training_options(MOFA_init)
	train_opts$seed <- 42
	train_opts$convergence_mode <- "fast"

	MOFA_prep1 <- prepare_mofa(MOFA_init,
		data_options = data_opts,
		model_options = model_opts,
		training_options = train_opts
	)

	MOFA_prep2 <- mofa2(
		mae_sub,
		assays = c("log1p", "exprs", "",'log1p'),
		num_factors = 4,
		seed = 42,
		convergence_mode = "fast"
	)
	
	expect_equivalent(MOFA_prep1,MOFA_prep2)
	expect_identical(MOFA_prep1,MOFA_prep2)

	# check manual parameter values
	MOFA_prep2 <- mofa2(
		mae_sub,
		assays = c("log1p", "exprs", "",'log1p'),
		num_factors = 4,
		seed = 1337,
		convergence_mode = "fast"
	)
	expect_equal(MOFA_prep2@training_options$seed,1337)
	expect_equal(MOFA_prep2@training_options$convergence_mode,"fast")
	expect_equal(MOFA_prep2@training_options$maxiter,1000)
	expect_equal(MOFA_prep2@model_options$spikeslab_factors,FALSE)
	
	# add: stochastic options 
	# add: mefisto options
	# add: SCE with altexperiments?
})