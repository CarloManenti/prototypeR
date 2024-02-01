is_python_package_installed <- function(packages.vec,
                                        envname='r-reticulate'){

    ### Description ###
    # Checks if you have installed a python package in a given virtual env
    # by default it checks for r-reticulate since it is suite for a SLURM HPC
    # If the package is not present, it directly tries to install it…

    # example usage
    # is_python_package_installed(c('numpy', 'pandas'))



    is_package_installed('reticulate')

    msg.chr <- paste0('Installing in ',envname , ' virtual environment packages : \n')
    # retrieving the specific of the virtual|conda environment

    if(!(envname %in% reticulate::virtualenv_list())){
      message('Defining a new virtual environment called : ', envname)
      reticulate::virtualenv_create(envname = envname)
    }

    env_packages.table <- reticulate::py_list_packages(envname = envname)
    env_packages.vec <- env_packages.table[, 'package']

    installed.boolvec <- packages.vec %in% env_packages.vec
    if(any(!installed.boolvec)){
        packages2install.vec <-  packages.vec[!installed.boolvec]
        message(paste0(msg.chr, packages2install.vec))
        reticulate::virtualenv_install(envname = envname,
                                       packages = packages2install.vec)
    }
    else(message('All python packages are installed'))
}
