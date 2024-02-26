load('data2test/sce.rda')

test_that("decomp_cnmf renaming works", {
  sce <- decomp_cnmf(sce, n_components = 10, levels = 10, num_iterations = 50)
  expect_equal(SingleCellExperiment::reducedDimNames(sce)[[2]], 'cnmf')
})

test_that('decomp_cnmf throws an error', {
  expect_error(decomp_cnmf(sce, n_components = 0, levels = 1, num_iterations = 10))
})
