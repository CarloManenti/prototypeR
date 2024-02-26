load('data2test/sce.rda')

test_that("decomp_kaa works", {
  sce <- decomp_kaa(sce, n_components = 5, n_highly_variable_genes = 250, n_pcs = 10, n_waypoint_eigs = 3, verbose = FALSE, plots = FALSE)
  expect_equal(is.null(SingleCellExperiment::reducedDim(sce, 'kaa')), FALSE)
})

test_that('decomp_kaa throws an error', {
  expect_error(decomp_kaa(sce, n_components = 10, n_highly_variable_genes = 3))
})

