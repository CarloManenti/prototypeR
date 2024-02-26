load('data2test/sce.rda')

test_that("decomp_va works", {
  sce <- decomp_va(sce, n_components = 5, accelerator = 'cpu', verbose = TRUE, max_epochs = 10)
  expect_equal(dim(SingleCellExperiment::reducedDim(sce, 'va')), c(1000, 5))
})

# if this test fails, it is much better… :)
# for now we can smile at the most complex batman song to be
# ever written… (for now Q1 2024)
test_that('decomp_dense_pca auto works', {
  expect_error(decomp_va(sce, n_components = 5, accelerator = 'gpu', verbose = TRUE, max_epochs = 10))
})
