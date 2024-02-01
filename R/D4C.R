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



sce2adata_sparse <- function(sce,
                             envname = 'r-metacell',
                             main_layer = 'counts',
                             ...){

  ### Description ###
  # Converts a sce into an anndata, taking care also of the
  # coversion of the CSC to CSR sparse matrix.

  # example usage
  # sce2adata_sparse(sce)

  is_package_installed('sceasy')
  is_package_installed('reticulate')
  is_python_package_installed(envname = envname, packages.vec = c('scipy'))
  reticulate::use_virtualenv(envname)
  scipy <- reticulate::import('scipy')

  # the matrix will be already transposed! cells x features
  adata <- sceasy::convertFormat(sce,
                                 from = "sce",
                                 to = "anndata",
                                 main_layer = main_layer,
                                 drop_single_values = FALSE,
                                 ...)

  # Converting a CSC sparse matrix to CSR sparse matrix
  # to run fuster the model
  if(scipy$sparse$csc$isspmatrix_csc(adata$X)){
    adata$X <- scipy$sparse$csr_matrix(adata$X)
  }

  return(adata)
}



center_and_scale <- function(matrix.dgCMatrix, center=FALSE, scale=FALSE){
  ### Description ###
  # Centers and scales features-wise (row-wise) a sparse matrix.
  # In case both center and scales are false, it returns the original matrix

  # example usage
  # center_and_scale(counts(sce), center = T, scale = T)

  if(any(center, scale)){
    # to compute feature wise centering | scaling
    # we need to transpose the matrix
    matrix.dense <- t(base::scale(t(matrix.dgCMatrix),
                                  center = center,
                                  scale = scale))
    matrix.dgCMatrix <- as(matrix.dense, 'sparseMatrix')
    # it might make no sense to have it sparse, since is filled with non 0s
    # but the type is consistent
  }
  return(matrix.dgCMatrix)
}


load('data/sce.rda')
