load('data2test/sce.rda')

test_that("partition_metacells works", {
  sce <- partition_metacells(sce, target_number_of_metacells = 2, min_umi = 5)
  expect_equal(!is.null(S4Vectors::metadata(sce)['metacells']), TRUE)
})

test_that('partition_metacells fails', {
  expect_error(partition_metacells(sce, target_number_of_metacells = 2, min_umi = 100000))
})
