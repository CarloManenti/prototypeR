list_to_assignment_matrix <- function(input.list){
    ### Description ###
    # Defines a sparse membership matrix (with 1 and 0s) features x observations
    # form a vector (or list).

    # example usage
    # input.list = c('a', 'a', 'b', 'c')
    # list_to_assignment_matrix(input.list)


    elements.vec <- sort(unique(unlist(input.list)))
    membership.matrix <- t(sapply(elements.vec, function(element_i){
                               as.numeric(input.list %in% element_i)
                           }))
    # giving names to the axis
    rownames(membership.matrix) <- elements.vec
    colnames(membership.matrix) <- colnames(input.list)
    # making it sparse
    membership.dgCMatrix <- as(membership.matrix, 'sparseMatrix')
    return(membership.dgCMatrix)
}
