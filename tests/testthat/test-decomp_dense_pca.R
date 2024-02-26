load('data2test/sce.rda')

test_that("decomp_dense_pca works", {
  sce <- decomp_dense_pca(sce, n_components = 10)
  expect_equal(dim(SingleCellExperiment::reducedDim(sce, 'dense-pca')), c(1000, 10))
})

# to avoid throwing even more warnings…
load('data2test/sce.rda')

test_that('decomp_dense_pca auto works', {
  sce <- decomp_dense_pca(sce, n_components = 10, method = 'auto')
  expect_equal(dim(SingleCellExperiment::reducedDim(sce, 'dense-pca')), c(1000, 10))
})
