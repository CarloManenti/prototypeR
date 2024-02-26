#' Hard Assignment of Latents to an Observation (row)
#'
#' This function performs hard assignment of latents (columns) to
#' observations (rows) by using the highest weight as criteria of assignment.
#'
#' @param matrix <dgCMatrix or matrix generally> The matrix H (observations x
#' latents) which will be used for the hard assignment.
#' @param margin <either 1 or 2> dimension used to pick the maximum, 1 to
#' work on the latents, so the rows!
#' @return Returns a vector with the assigned latent for each cell.
#' @examples
#' #hard_assignment(SingleCellExperiment::reducedDim(sce, 'pca'))
#' @export
hard_assignment <- function(matrix, margin = 1){
    ### Description ###

    # Takes a Matrix cells x prototypes
    # it returns a vector in which each cell
    # has the maximum contributing prototype assign to it

    hard_assignment <- apply(matrix,
                             margin,
                             function(col_i){
                                 # to avoid double assingment
                                 # this could be done better!
                                 min(which(col_i == max(col_i)))
                             }
                       )
    return(hard_assignment)
}

