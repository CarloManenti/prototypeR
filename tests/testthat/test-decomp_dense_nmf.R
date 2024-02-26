load('data2test/sce.rda')

test_that("decomp_dense_nmf works", {
  sce <- decomp_dense_nmf(sce, n_components = 10)
  expect_equal(dim(SingleCellExperiment::reducedDim(sce, 'dense-nmf')), c(1000, 10))
})

test_that('decomp_dense_nmf throws an error', {
  expect_error(decomp_dense_nmf(sce, n_components = -1))
})
