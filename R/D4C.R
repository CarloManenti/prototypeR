change_default_name <- function(result_name, name_list){
  ### Description ###
  # Changes the result_name to avoid overriding data

  # example usage
  # change_default_name('my_method', reduceDimNames(sce))

  if(result_name %in% name_list){
    # setting the name of the matrix if you have multiple default names
    method_n <- sum(grepl(result_name, name_list))
    result_name <- paste0(result_name, method_n)
  }
  return(result_name)
}

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
  return(TRUE)
}

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

is_python_package_installed <- function(packages.vec, envname='r-reticulate'){

  ### Description ###
  # Checks if you have installed a python package in a given virtual env
  # by default it checks for r-reticulate since it is suite for a SLURM HPC
  # If the package is not present, it directly tries to install it…

  # example usage
  # is_python_package_installed(c('numpy', 'pandas'))



  is_package_installed('reticulate')

  msg.chr <- 'Installing in r-reticulate virtual environment packages : \n'
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


package_warning <- function(package_x){
  ### Description ###
  # Prints a warning message to install a package

  # example usage
  # package_warning('ruler')



  message('Please install ',
          package_x,
          ' to your home directory using install.packes(\'',
          package_x,
          '\', lib = \'</home/user.name/R>\')')
}


load('data/sce.rda')
