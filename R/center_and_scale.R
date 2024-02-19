#' Performs Column Wise Centering and Scaling of a Matrix
#'
#' This function taken as input a sparse matrix, generally dgCMatrix, and may
#' center and scale it in a colum-wise fashion, or simply return the matrix
#' as it is.
#'
#' @param matrix.dgCMatrix <matrix> Input matrix, possibly in sparse format
#' @param center=FALSE <bool> Whether to center the input matrix
#' @param scale=FALSE <bool> Whether to scale the input matrix
#' @return The input matrix in dgCMatrix sparse format after centering,
#' scaling, both or just as it is depending on the falg center and scale.
#'
#' @examples
#' center_and_scale(logcounts(sce), center = TRUE, scale = TRUE)
#' @export
center_and_scale <- function(matrix.dgCMatrix, center=FALSE, scale=FALSE){
    ### Description ###
    # Centers and scales features-wise (row-wise) a sparse matrix.
    # In case both center and scales are false, it returns the original matrix

    # example usage
    # center_and_scale(counts(sce), center = T, scale = T)

    if(any(center, scale)){
      # to compute feature wise centering | scaling
      # we need to transpose the matrix
      matrix.dense <- t(base::scale(t(matrix.dgCMatrix),
                                        center = center,
                                        scale = scale))
      matrix.dgCMatrix <- as(matrix.dense, 'sparseMatrix')
      # it might make no sense to have it sparse, since is filled with non 0s
      # but the type is consistent
    }
    return(matrix.dgCMatrix)
}
