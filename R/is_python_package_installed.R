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
  }else{
    #(message('All python packages are already installed'))
  }

}
