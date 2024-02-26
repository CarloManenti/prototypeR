load('data2test/sce.rda')

test_that("partition_graph leiden works", {
  sce <- partition_graph(sce)
  expect_equal(!is.null(SummarizedExperiment::colData(sce)['leiden']), TRUE)
})

test_that('partition_graph louvain works', {
  sce <- partition_graph(sce, method = 'louvain', result_name = 'louvain')
  expect_equal(!is.null(SummarizedExperiment::colData(sce)['louvain']), TRUE)
})
