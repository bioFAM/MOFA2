context("Making plots")
library(BioFAMtools)

test_that("plot data overview works", {
	expect_error(plot_data_overview(load_model("test_biofam.hdf5")), NA)
})

test_that("plot data heatmap", {
	expect_error(plot_data_heatmap(load_model("test_biofam.hdf5"), view = 1, factor = 1), NA)
})

test_that("plot data scatter", {
	expect_error(plot_data_scatter(load_model("test_biofam.hdf5"), view = 1, factor = 1), NA)
})

test_that("plot data ASCII in terminal", {
	expect_error(plot_ascii_data(load_model("test_biofam.hdf5")), NA)
})
