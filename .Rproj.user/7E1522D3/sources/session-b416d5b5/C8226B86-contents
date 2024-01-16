split_sparse <- function(sparse_matrix,
                         split_vector,
                         by_col = TRUE,
                         sorted = T){
    ### Description ###

    # Splits a sparse matrix (dgCMatrix from package Matrix)
    # based on a vector
    # By default it splits by col and the splits are ordered
    # based on the sort(unique(split_vector))
    # If it can not split by col it splits by row.

    # example usage
    # getting a Single Cell Experiment object
    # sce <- readRDS('~/sce.rds')
    # split_sparse(sparse_matrix = counts(sce.obj), split_vector = clusters)



    # initializing parameters
    transpose <- FALSE

    # do we want to split by col or means?
    split_dimension <- ifelse(by_col, 2, 1)

    # how to split the matrix? if by default we split by col
    # do we need to transpose the sparse matrix to split by column ?
    dimension <- dim(sparse_matrix)[split_dimension]
    transpose <- ifelse(dimension == length(split_vector), FALSE, TRUE)
    # yes we need
    if(transpose) sparse_matrix <- t(sparse_matrix)

    # getting the 'levels'
    splits <- unique(split_vector)

    # sorting ?
    if(sorted) splits <- sort(splits)

    # splitting
    # might be better ot use vapply with S4 FUN.VALUE
    splitted_list <- lapply(splits,
                            function(split_i){
                            split_i.idx <- which(split_vector == split_i)
                            mat <- sparse_matrix[, split_i.idx]
                            # have we transposed the matrix?
                            # YES -> we transpose again
                            # NO -> we leave it as it is
                            if(transpose) mat <- t(mat)
                            return(mat)
                            }
                      )

    # setting the names
    names(splitted_list) <- splits

    return(splitted_list)
}
