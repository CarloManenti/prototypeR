#' Splits a Sparse Matrix (without converting it into a dense one)
#'
#' This function splits a sparse matrix (generally dgCMatrix) without coercing
#' it into a dense representation.
#'
#' @param sparse_matrix <matrix, better if in sparse format> The matrix to be
#' split.
#' @param split_vector <vector> vector used to split the sparse_matrix into
#' smaller ones. The split is by default by column (by_col = TRUE); Hence the
#' split_vector must be the same length as the number of columns.
#' @param by_col <bool> default 'X'; Whether to split the matrix by columns or
#' by rows.
#' @param sorted <bool> default TRUE; Whether to sort the output or not.
#' This parameter is useful when dealing with a numeric split_vector, since the
#' if sorted = FALSE, the order of the split vector will define the order of
#' the output matrices.
#' @return A list with many smaller matrices, one for each unique element in
#' the split_vector. The order of the matrices depends on the sorted variable.
#' @examples
#' library(packageX)
#' data(sce)
#' split_sparse(sparse_matrix = SingleCellExperiment::counts(sce), split_vector = sce$cluster)
#' @export
split_sparse <- function(sparse_matrix,
                         split_vector,
                         by_col = TRUE,
                         sorted = TRUE){
    ### Description ###

    # Splits a sparse matrix (dgCMatrix from package Matrix)
    # based on a vector
    # By default it splits by col and the splits are ordered
    # based on the sort(unique(split_vector))
    # If it can not split by col it splits by row.



    # initializing parameters
    transpose <- FALSE

    # do we want to split by col or means?
    split_dimension <- ifelse(by_col, 2, 1)

    # how to split the matrix? if by default we split by col
    # do we need to transpose the sparse matrix to split by column ?
    dimension <- dim(sparse_matrix)[split_dimension]
    transpose <- ifelse(dimension == length(split_vector), FALSE, TRUE)
    # yes we need
    if(transpose){
        sparse_matrix <- Matrix::t(sparse_matrix)
    }

    # getting the 'levels'
    splits <- unique(split_vector)

    # sorting ?
    if(sorted){
        splits <- sort(splits)
    }

    # splitting
    # might be better ot use vapply with S4 FUN.VALUE
    splitted_list <- lapply(splits,
                            function(split_i){
                                split_i.idx <- which(split_vector == split_i)
                                mat <- sparse_matrix[, split_i.idx]
                                # have we transposed the matrix?
                                # YES -> we transpose again
                                # NO -> we leave it as it is
                                if(transpose) mat <- Matrix::t(mat)
                                return(mat)
                            })

    # setting the names
    names(splitted_list) <- splits

    return(splitted_list)
}
