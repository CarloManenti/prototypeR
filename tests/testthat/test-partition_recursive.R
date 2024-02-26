load('data2test/sce.rda')

test_that('partition_recursive fails', {
  expect_error(partition_recursive(sce, mc.cores = 1, n_iter = 1))
})
