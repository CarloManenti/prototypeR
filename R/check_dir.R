#' Checks if a Directory Exists, if it does not it creates it
#'
#' This function taken as input a directory and checks if it exists, if it does
#' not it creates it.
#'
#' @param directory2check.dir <character> path to the directory you want to
#' check
#' @return either nothing, if the directory already exists, or it creates
#' the directory and prompt the user with a message of where the directory has
#' been created
#'
#' @examples
#' check_dir('PhD')
#' @export
check_dir <- function(directory2check.dir){
    ### Description ####
    # Checks if the directory exists;
    # if not, it creates the directory

    # example usage
    # check_dir('~/Documents/PhD')



    if(!dir.exists(directory2check.dir)){
        message(paste0('--- new directory defined at : ',
                       directory2check.dir , '---\n'))
        dir.create(directory2check.dir)
    }
}

