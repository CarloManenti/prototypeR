store_W <- function(sce,
                    w.matrix,
                    latent_name,
                    result_name
                    ){

  ### Description ###
  # Handy function to store the W matrix genes x cells,
  # in the metadata field of a SingleCellExperiment object
  # as a sparse dgCMatrix with row and column names.

  # example usage
  # store_W(sce, w.matrix, 'pca', 'PC)

  n_latents <- ncol(w.matrix)
  # setting the names to the H matrix
  colnames(w.matrix) <- paste0(rep(latent_name, n_latents),
                               seq_len(n_latents))
  rownames(w.matrix) <- rownames(sce)
  # setting the type of the H matrix
  w.dgCMatrix <- as(w.matrix, 'sparseMatrix')
  # checking for an already in use result_name
  result_name <- change_default_name(result_name, reducedDimNames(sce))

  # storing the actual object
  metadata(sce)[[result_name]] <- w.dgCMatrix
  return(sce)
}
