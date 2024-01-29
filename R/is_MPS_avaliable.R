is_MPS_avaliable <- function(envname){
    ### Description ###
    # Checks if MPS accelaration is avaliable!

    #is_MPS_avaliable(envname = 'r-decomp')

    is_python_package_installed(envname = envname, packages.vec = c('torch'))
    torch <- reticulate::import('torch')

    if(torch$backends$mps$is_available()){
        mps_device = torch$device("mps")
        message('--- MPS device avaliable ----')
    }else{message("--- MPS device not found ---")}

    }
