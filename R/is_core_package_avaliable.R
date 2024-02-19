is_core_package_avaliable <- function(core_package = 'torch')
python_configuration <- reticulate::py_config()
numpy_version <- python_configuration[['numpy']][['version']]


configuration.msg <- paste0('Please use a virtual environment in which you have an avaliable ',
                            core_package,
                            ' version by using set_virtual_env(\'my-python-env\')')

if(is.null(numpy_version)){
  stop(configuration.msg)
}



python_config <- reticulate::py_discover_config()
avliable_python <- python_config$python
reticulate::use_python(python = avliable_python, require = TRUE)
reticulate::py_config()



