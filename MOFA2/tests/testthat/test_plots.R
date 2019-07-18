# context("Making plots")
# library(MOFA2)
# 
# # Data plots
# 
# test_that("plot data overview", {
# 	expect_error(plot_data_overview(load_model("test_mofa.hdf5")), NA)
# })
# 
# test_that("plot data heatmap", {
# 	expect_error(plot_data_heatmap(load_model("test_mofa.hdf5"), view = 1, factor = 1), NA)
# })
# 
# test_that("plot data scatter", {
# 	expect_error(plot_embeddings(load_model("test_mofa.hdf5"), view = 1, factor = 1), NA)
# })
# 
# # test_that("plot data ASCII in terminal", {
# # 	expect_error(plot_ascii_data(load_model("test_mofa.hdf5")), NA)
# # })
# 
# 
# # Plotting weights
# 
# test_that("plot weights heatmap", {
# 	expect_error(plot_weights_heatmap(load_model("test_mofa.hdf5"), view = 1), NA)
# })
# 
# test_that("plot weights scatter", {
# 	expect_error(plot_weights_scatter(load_model("test_mofa.hdf5"), view = 1, factors = 1:2), NA)
# })
# 
# test_that("plot weights", {
# 	expect_error(plot_weights(load_model("test_mofa.hdf5"), factor = 1), NA)
# })
# 
# test_that("plot top weights", {
# 	expect_error(plot_top_weights(load_model("test_mofa.hdf5"), view = 1, factor = 1), NA)
# })
# 
# # Plotting factor values
# 
# test_that("plot factor values", {
# 	expect_error(plot_factors(load_model("test_mofa.hdf5")), NA)
# 	expect_error(plot_factors(load_model("test_mofa.hdf5"), factors = 1:2), NA)
# 	expect_error(plot_factors(load_model("test_mofa.hdf5"), factors = 1:2, color_by = 'group_name'), NA)
# 	expect_error(plot_factors(load_model("test_mofa.hdf5"), factors = 1:2, shape_by = 'group_name'), NA)
# 	expect_error(plot_factors(load_model("test_mofa.hdf5"), factors = 1:2, group_by = 'group_name'), NA)
# })
# 
# test_that("plot factor values", {
# 	expect_error(plot_factors(load_model("test_mofa.hdf5"), factors = 1:2), NA)
# 	expect_error(plot_factors(load_model("test_mofa.hdf5"), factors = 1:2, color_by = 'group_name'), NA)
# })
# 
# test_that("plot factors correlation", {
# 	expect_error(plot_factor_cor(load_model("test_mofa.hdf5")), NA)
# })
# 
