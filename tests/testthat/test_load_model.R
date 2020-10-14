context("Loading the model")
library(MOFA2)

test_that("a pre-trained model can be loaded from disk", {
  filepath <- system.file("extdata", "model.hdf5", package = "MOFA2")
  expect_is(load_model(filepath), "MOFA")
})

