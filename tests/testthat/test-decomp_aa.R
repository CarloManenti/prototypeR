load('data2test/sce.rda')

test_that('decomp_aa uses the specified reduced dimensionality', {
  sce <- decomp_aa(sce, 5, reduced_representation = 'pca', boostrap_number = 2)
  expect_equal(SingleCellExperiment::reducedDimNames(sce)[2], 'aa')
})


test_that("decomp_aa fails", {
    expect_error(decomp_aa(sce, 5, boostrap_number = 2, reduced_representation = '#404'))
})
