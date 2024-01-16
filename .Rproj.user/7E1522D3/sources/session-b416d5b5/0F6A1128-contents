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
