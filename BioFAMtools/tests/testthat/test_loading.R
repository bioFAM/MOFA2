context("Loading the model")
library(BioFAMtools)

test_that("a pre-trained model can be loaded from disk", {
  expect_is(load_model("test_biofam.hdf5"), "BioFAModel")
})

