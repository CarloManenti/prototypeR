load('data2test/sce.rda')

test_that("decomp_sparse_pca renaming works", {
    sce <- decomp_sparse_pca(sce, n_components = 10, result_name = 'new-pca')
    expect_equal(SingleCellExperiment::reducedDimNames(sce)[[2]], 'new-pca')
})


test_that('decomp_sparse_pca throws a warning', {
  expect_warning(decomp_sparse_pca(sce, n_components = 10))
})
