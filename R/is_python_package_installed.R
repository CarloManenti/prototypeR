#' Checks if a Python Package is Available in a given Virtual Environment
#'
#' Handy function to check if a Python package is available in a specified
#' virtual environment or not. If the virtual environment does not exists, it
#' will create one with the same name. If the python package is not to be found,
#' it will install it. It is based on reticulate!
#'
#' @param packages.vec <vector of characters> Vector of the needed python
#' packages.
#' @param envname <character> default 'r-reticulate'; name of the virtual
#' enviroment to be used. If the virtual environment already exists, it will
#' use it to check if the packages are available. Otherwise it will create
#' a new one with the specified name. The new environment can be removed with
#' *reticulate::virtualenv_remove(envname = 'my-env')* and the list of all
#' available packages can be seen using *reticulate::virtualenv_list()*.
#' @param forced <bool> default TRUE; Whether to force the use of the specified
#' virtual environment or not. Generally you want to specify the use of the
#' virtual environment if you want to use right away the packages.
#' @return nothing, it is a simple side-effect function. It will break the
#' function if the package is not available. But it will fail gracefully.
#' @examples
#' is_python_package_installed(packages.vec = 'numpy', envname = 'r-decomp')
#' @export
is_python_package_installed <- function(packages.vec,
                                        envname='r-reticulate',
                                        forced=TRUE){

  ### Description ###
  # Checks if you have installed a python package in a given virtual env
  # by default it checks for r-reticulate since it is suite for a SLURM HPC
  # If the package is not present, it directly tries to install it…

  # you should always load python packages by :
  # 1st declaring which virtual env you want to
  #     use via reticulate::use_virtualenv('my-virtualenv')
  # 2nd importing the packages via reticulate::import('my-package')
  # this is why forced is a variable, since it is forcing this execution order

  # example usage
  # is_python_package_installed(c('numpy', 'pandas'), envname = 'r-decomp')

  is_package_installed('reticulate')

  # initializing installing message
  msg.chr <- paste0('Installing in ',
                    envname,
                    ' virtual environment packages : \n')

  # retrieving the specific of the virtual environment
  if(!(envname %in% reticulate::virtualenv_list())){
    message('Defining a new virtual environment called : ', envname)
    reticulate::virtualenv_create(envname = envname)
    # forcing the new virtual environment
    reticulate::use_virtualenv(envname)
  }

  # checking if the packages needed are already installed in the given env.
  # forced is useful to use the correct virtual enviroment BEFORE loading the
  # packages! Otherwise, it may lead to errors in the code.
  if(forced == TRUE){
    reticulate::use_virtualenv(envname)
  }
  env_packages.table <- reticulate::py_list_packages(envname = envname)
  env_packages.vec <- env_packages.table[, 'package']
  installed.boolvec <- packages.vec %in% env_packages.vec

  # installing new packages
  if(any(!installed.boolvec)){
    packages2install.vec <-  packages.vec[!installed.boolvec]
    message(paste0(msg.chr, packages2install.vec))
    reticulate::virtualenv_install(envname = envname,
                                   packages = packages2install.vec)
  }

}
