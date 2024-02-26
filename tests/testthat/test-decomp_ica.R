load('data2test/sce.rda')

test_that("decomp_ica works", {
  sce <- decomp_ica(sce, n_components = 10)
  expect_equal(dim(SingleCellExperiment::reducedDim(sce, 'ica')), c(1000, 10))
})

# to avoid throwing even more warnings…
load('data2test/sce.rda')

test_that('decomp_ica jade works', {
  sce <- decomp_ica(sce, n_components = 10, method = 'jade')
  expect_equal(dim(SingleCellExperiment::reducedDim(sce, 'ica')), c(1000, 10))
})
