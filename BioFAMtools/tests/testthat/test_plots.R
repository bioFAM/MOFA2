context("Making plots")
library(BioFAMtools)

# Data plots

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


# Plotting weights

test_that("plot weights heatmap", {
	expect_error(plot_weights_heatmap(load_model("test_biofam.hdf5"), view = 1), NA)
})

test_that("plot weights scatter", {
	expect_error(plot_weight_scatter(load_model("test_biofam.hdf5"), view = 1, factors = 1:2), NA)
})

test_that("plot weights", {
	expect_error(plot_weights(load_model("test_biofam.hdf5"), factor = 1), NA)
})

test_that("plot top weights", {
	expect_error(plot_top_weights(load_model("test_biofam.hdf5"), view = 1, factor = 1), NA)
})

test_that("plot weights correlation", {
	expect_error(plot_weight_cor(load_model("test_biofam.hdf5"), view = 1), NA)
})

# Plotting factor values

test_that("plot factor values", {
	expect_error(plot_factor(load_model("test_biofam.hdf5")), NA)
})

test_that("plot factor values", {
	expect_error(plot_factors(load_model("test_biofam.hdf5"), factors = 1:2), NA)
})

test_that("plot factors correlation", {
	expect_error(plot_factor_cor(load_model("test_biofam.hdf5")), NA)
})