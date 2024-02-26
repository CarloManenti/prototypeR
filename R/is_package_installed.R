#' Checks if an R Package is Available
#'
#' Handy function to check if an R package can be used or not.
#' If the specified R package is nowhere to be found it will fail
#' gracefully.
#'
#' @param package_x <character> Name of the needed package
#' @param lib.loc <NULL or character> default NULL; Path to the R library.
#' If NULL, it will use the default R library. Instead if specified it will
#' use the one specified. This could be particular useful when working on
#' a shared R library and a user specific R library in a High Performance
#' Computing (HPC).
#' @return nothing, it is a simple side-effect function. It will break the
#' function if the package is not available. But it will fail gracefully.
#' @export
is_package_installed <- function(package_x, lib.loc=NULL){
    ### Description ###
    # Checks if a package is installed or not in your general .R folder

    if(!requireNamespace(package_x, quietly = TRUE, lib.loc=NULL)){
        stop(package_warning(package_x))
    }
    # no return | side effect function
}
