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
