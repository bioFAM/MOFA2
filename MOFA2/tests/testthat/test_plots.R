context("Making plots")
library(MOFA2)

filepath <- system.file("extdata", "model.hdf5", package = "MOFA2")
test_mofa2 <- load_model(filepath)

# Data plots

test_that("plot data overview works", {
	expect_silent(p <- plot_data_overview(test_mofa2))
})

test_that("plot data heatmap", {
	expect_silent(p <- plot_data_heatmap(test_mofa2, view = 1, factor = 1, silent = TRUE))
})

test_that("plot data scatter", {
	expect_silent(p <- plot_data_scatter(test_mofa2, view = 1, factor = 1))
})

test_that("plot data ASCII in terminal", {
	expect_error(plot_ascii_data(test_mofa2), NA)
})


# Plotting weights

test_that("plot weights heatmap", {
	expect_silent(p <- plot_weights_heatmap(test_mofa2, view = 1, silent = TRUE))
})

test_that("plot weights", {
	# For multiple factors
	expect_silent(p <- plot_weights(test_mofa2, view = 1, factors = 1:2))
	# For one factor
	expect_silent(p <- plot_weights(test_mofa2, factors = 1))
})

test_that("plot top weights", {
	expect_silent(p <- plot_top_weights(test_mofa2, view = 1, factors = 1))
})



# Plotting factor values

test_that("plot factor values", {
	expect_silent(p <- plot_factor(test_mofa2))
})

test_that("plot factor values", {
	expect_silent(p <- plot_factors(test_mofa2, factors = 1:2))
})

test_that("plot factors correlation", {
	expect_error({plot_factor_cor(test_mofa2); dev.off()}, NA)
})
