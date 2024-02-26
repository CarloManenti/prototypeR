#' It transforms a List into an Assignment Matrix
#'
#' This function takes a list and returns a sparse assignment matrix.
#'
#' @param input.list <list or vector> the input list which will be transformed
#' into an assignment matrix
#' @return a sparse (dgCMatrix) assignment matrix
#' @examples
#' list_to_assignment_matrix(input.list = c('a', 'a', 'b', 'c'))
#' @export
list_to_assignment_matrix <- function(input.list){
    ### Description ###
    # Defines a sparse membership matrix (with 1 and 0s) features x observations
    # form a vector (or list).



    elements.vec <- sort(unique(unlist(input.list)))
    membership.matrix <- Matrix::t(sapply(elements.vec, function(element_i){
                              as.numeric(input.list %in% element_i)
                           }))
    # giving names to the axis
    rownames(membership.matrix) <- elements.vec
    colnames(membership.matrix) <- colnames(input.list)
    # making it sparse
    membership.dgCMatrix <- methods::as(membership.matrix, 'sparseMatrix')
    return(membership.dgCMatrix)
}
