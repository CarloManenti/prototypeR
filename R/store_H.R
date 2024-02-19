store_H <- function(sce,
                    h.matrix,
                    result_name,
                    latent_name){

  ### Description ###
  # Stores the H matrix into a reduceDim space of the SingleCellExperiment
  # object as a sparse dgCMatrix. It also assigns names to both columns and
  # rows.

  # example usage
  # store_H(sce, pca_h.matrix, 'pca', 'PC')



  # initializing values
  n_latents <- ncol(h.matrix)
  # setting the names to the H matrix
  colnames(h.matrix) <- paste0(rep(latent_name, n_latents),
                               seq_len(n_latents))
  rownames(h.matrix) <- colnames(sce)
  # setting the type of the H matrix
  h.dgCMatrix <- as(h.matrix, 'sparseMatrix')
  # checking for an already in use result_name
  result_name <- change_default_name(result_name, reducedDimNames(sce))

  # storing the actual object
  reducedDim(sce, result_name) <- h.dgCMatrix
  return(sce)
}
