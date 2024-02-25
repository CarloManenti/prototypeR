#' Handy function for Warnings Messages
#'
#' This function takes the name of an R package we want to install and returns
#' a nice message warning on how to install it. It is a general message, but
#' it makes for a graceful break if the package has not been installed.
#' @param package_x <character> the name of the package need to be installed
#' @return a nice message providing instruction on how to install it
#' @examples
#' package_warning('Matrix')
#' @export
package_warning <- function(package_x){
    ### Description ###
    # Prints a warning message to install a package

    # example usage
    # package_warning('ruler')



    msg <- paste0('Please install ',
                  package_x,
                  ' to your home directory using install.packes(\'',
                  package_x,
                  '\', lib = \'</home/user.name/R>\')')
    return(msg)
}
