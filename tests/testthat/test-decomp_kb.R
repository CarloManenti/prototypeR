load('data2test/sce.rda')

prior_knowledge = list('0' = list('gene_set_1' = list("KLHL17",
                                                      "CCNL2",
                                                      "ATAD3B",
                                                      "NOL9"),
                                  'gene_set_2' = list("CENPS-CORT",
                                                      "FBXO2",
                                                      "KLHL21")),
                       '1' = list('gene_set_1' = list("NPPA",
                                                      "FHAD1",
                                                      "FBXO42",
                                                      "ATP13A2"),
                                  'gene_set_3' = list("PADI2",
                                                      "AHDC1",
                                                      "PPP1R8",
                                                      "LAPTM5")),
                       # set of genes which will be evaluated for all the
                       # cells in data…
                       'global' = list('gene_set_4' = list("AZIN2",
                                                           "ZSCAN20",
                                                           "ZC3H12A")))

test_that("decomp_kb works", {
  sce <- decomp_kb(sce, metadata_key = 'is.doublet', prior_knowledge = prior_knowledge, n_hvgs = NULL, num_epochs = 100)
  expect_equal(SingleCellExperiment::reducedDimNames(sce)[2], 'kb')
})


test_that('decomp_kb fails', {
  expect_error(decomp_kb(sce, prior_knowledge = NULL))
})
