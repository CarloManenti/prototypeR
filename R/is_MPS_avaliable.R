#' Checks if MPS Acceleration is Available
#'
#' This function checks if the MPS acceleration is available. It is
#' participatory useful when working with Neural Networks like a Variational
#' AutoEncoder.
#'
#' @param envname <character> default 'r-decomp';
#' Specify the name of the python virtual
#' environment to be used. If it does not exists it will create one and use it.
#' @param verbose <bool> default FALSE; Whether to be prompted with message
#' for each step of the analysis.
#' @return nothing, it is a simple side-effect function. If verbose is set to
#' TRUE, it will print a message to the end user stating whether it could
#' use MPS acceleration or not.
#' @examples
#' #is_MPS_avaliable(envname = 'r-decomp', verbose = TRUE)
#' @export
is_MPS_avaliable <- function(envname, verbose){
    ### Description ###
    # Checks if MPS accelaration is avaliable!

    #is_MPS_avaliable(envname = 'r-decomp', verbose = FALSE)

    is_python_package_installed(envname = envname,
                                packages.vec = c('torch'))
    torch <- reticulate::import('torch', delay_load = TRUE)

    if(torch$backends$mps$is_available()){
        mps_device <- torch$device("mps")
        if(verbose){
            message('--- MPS device avaliable ----')
        }
    }else{
        if(verbose){
            message("--- MPS device not found ---")}
        }
    }
