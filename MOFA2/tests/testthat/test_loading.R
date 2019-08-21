context("Loading the model")
library(MOFA2)

test_that("a pre-trained model can be loaded from disk", {
  expect_is(load_model("test_mofa2.hdf5"), "MOFA")
})

