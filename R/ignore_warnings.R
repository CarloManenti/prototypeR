#' Ignore Warnings
#'
#' This function sets the python enviroment to ignore warnings if provided
#' with ignore_warning = TRUE.
#'
#' @param ignore_warnings <bool> Whether to ignore warning related to the
#' python enviroment.
#' @param envname <character> default 'r-decomp';
#' Specify the name of the python virtual
#' environment to be used. If it does not exists it will create one and use it.
#' @param verbose <bool> default FALSE; Whether to be prompted with message
#' for each step of the analysis.
#' @return nothing, it is a simple side-effect function.
#' @examples
#' Ignore_Warnings(TRUE)
#' @export

ignore_warnings <- function(ignore_warnings,
                            envname = 'r-decomp',
                            verbose = FALSE){
    if(ignore_warnings){
        is_package_installed('reticulate')
        reticulate::use_virtualenv(envname)
        warnings <- reticulate::import('warnings', delay_load = TRUE)
        warnings$filterwarnings("ignore")
        if(verbose){
            message('warnings are now suppressed')
        }
    }
}
