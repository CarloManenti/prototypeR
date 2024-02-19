is_package_installed <- function(package_x, lib.loc=NULL){
    ### Description ###
    # Checks if a package is installed or not in your general .R folder

    # example usage
    # is_package_installed('ACTIONet') -> TRUE
    # is_package_installed('ruler') -> …
    # Please install ruler to your home directory
    # using install.packes('ruler', lib = '</home/user.name/R>')
    # Error in is.package_x.install("ruler")



    if(!requireNamespace(package_x, quietly = TRUE, lib.loc=NULL)){
        stop(package_warning(package_x))
    }
    # no return | side effect function
}
