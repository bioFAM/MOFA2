context("Creating the model from different object")
library(BioFAMtools)

test_that("a model can be created from a list of matrices", {
	m <- as.matrix(read.csv('matrix.csv'))
	expect_is(create_biofam(list("view1" = m)), "BioFAModel")
	expect_warning(create_biofam(list("view1" = m)))
	expect_error(create_biofam(m))
})

