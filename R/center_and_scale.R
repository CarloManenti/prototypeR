#' Performs Column Wise Centering and Scaling of a Matrix
#'
#' This function taken as input a sparse matrix, generally dgCMatrix, and may
#' center and scale it in a column-wise fashion, or simply return the matrix
#' as it is.
#'
#' @param matrix.dgCMatrix <matrix> Input matrix, possibly in sparse format
#' @param center <bool> default FALSE; Whether to center the input matrix
#' @param scale <bool> default FALSE; Whether to scale the input matrix
#' @return The input matrix in dgCMatrix sparse format after centering,
#' scaling, both or just as it is depending on the flag center and scale.
#' @examples
#' library(packageX)
#' library(SingleCellExperiment)
#' data(sce)
#' center_and_scale(SingleCellExperiment::logcounts(sce), center = TRUE, scale = TRUE)
#' @export
center_and_scale <- function(matrix.dgCMatrix, center=FALSE, scale=FALSE){
    ### Description ###
    # Centers and scales features-wise (row-wise) a sparse matrix.
    # In case both center and scales are false, it returns the original matrix



    if(any(center, scale)){
      # to compute feature wise centering | scaling
      # we need to transpose the matrix
      matrix.dense <- Matrix::t(base::scale(Matrix::t(matrix.dgCMatrix),
                                            center = center,
                                            scale = scale))
      matrix.dgCMatrix <- methods::as(matrix.dense, 'sparseMatrix')
      # it might make no sense to have it sparse, since is filled with non 0s
      # but the type is consistent
    }
    return(matrix.dgCMatrix)
}
