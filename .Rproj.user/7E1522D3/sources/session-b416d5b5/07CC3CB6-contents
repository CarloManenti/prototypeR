hard_assignment <- function(matrix, margin = 1){
    ### Description ###

    # Taken a Matrix cells x prototypes
    # it returns a vector in which each cell
    # has the maximum contributing prototype assign to it

    # example usage
    # get the data.
    # sce <- readRDS('~/sce.rds')
    # normalize and reduce
    # ace <- normalize.ace(sce)
    # ace <- reduce.ace(ace)
    # hard_assignment(ace$ACTION) # on the PCA H matrix




    hard_assignment <- apply(matrix,
                             margin,
                             function(col_i){
                                 # to avoid double assingment
                                 # I select the mimum value among all
                                 min(which(col_i == max(col_i)))
                             }
                       )
    return(hard_assignment)
}

