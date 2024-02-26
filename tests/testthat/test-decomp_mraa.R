load('data2test/sce.rda')

test_that('decomp_mraa throws an error', {
  expect_error(decomp_mraa(sce, n_components = 10))
})
