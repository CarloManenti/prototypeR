#' Handy Function to Store the W matrix in a SingleCellExperiment Object
#'
#' This function stores a dense W matrix into a SingleCellExperiment (SCE)
#' object in the metadata slot in a sparse format (dgCMatrix). It also
#' adds the names to the rows and changes the name of the stored W matrix
#' if it is already present that name inside the reducedDimNames slot.
#' Hence you want to use the store_W before the store_H function!
#'
#' @param sce <SingleCellExperiment object> SCE object
#' @param w.matrix <matrix generally dense> w matrix (genes x latents) to be
#' stored in the SCE object.
#' @param result_name <character> Name used to store the result in the
#' SCE object.
#' @param latent_name <character> Name used for the latents. It will be used
#' as column name followed by the number of the latent. (example AA; AA1 AA2 …)
#' @return The SCE object with the new W representation stored in the
#' metadata slot.
#' @examples
#' library(packageX)
#' data(sce)
#' store_W(sce, w.matrix =  Matrix::t(matrix(stats::rnorm(2 * nrow(sce)), nrow = nrow(sce))), result_name = 'aa', latent_name = 'AA')
#' @export
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
    w.dgCMatrix <- methods::as(w.matrix, 'sparseMatrix')
    # checking for an already in use result_name
    result_name <- change_default_name(result_name, SingleCellExperiment::reducedDimNames(sce))

    # storing the actual object
    S4Vectors::metadata(sce)[[result_name]] <- w.dgCMatrix
    return(sce)
}
