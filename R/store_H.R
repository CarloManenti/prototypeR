#' Handy Function to Store the H matrix in a SingleCellExperiment Object
#'
#' This function stores a dense H matrix into a SingleCellExperiment (SCE)
#' object in the reducedDim slot in a sparse format (dgCMatrix). It also
#' adds the names to the rows and changes the name of the stored H matrix
#' if it is already present that name inside reducedDim slot.
#'
#' @param sce <SingleCellExperiment object> SCE object
#' @param h.matrix <matrix generally dense> H matrix (cells x latents) to be
#' stored in the SCE object.
#' @param result_name <character> Name used to store the result in the
#' SCE object.
#' @param latent_name <character> Name used for the latents. It will be used
#' as column name followed by the number of the latent. (example AA; AA1 AA2 …)
#' @return The SCE object with the new H representation stored in the
#' reducedDim slot.
#' @examples
#' store_H(sce, h.matrix = aa_h.matrix, result_name = 'aa', latent_name = 'AA')
#' @export
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
