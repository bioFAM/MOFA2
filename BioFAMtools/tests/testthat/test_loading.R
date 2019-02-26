context("Loading the model")
library(BioFAMtools)

test_that("a pre-trained model can be loaded from disk", {
  if (file.exists("test_biofam.hdf5"))
    expect_is(load_model("test_biofam.hdf5"), "BioFAModel")
})

