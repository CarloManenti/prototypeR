#' Changes Default Result Name only if it is already present
#'
#' This function changes the default result name if it is already present
#' in a named list.
#'
#' @param result_name <character> default result_name.
#' @param name_list <list> named list which should not
#' contain the default result_name as one of the names.
#' @return Either the same result_name <character > or a new result_name which
#' accounts for the presence of multiple result_names similar to the one passed
#' to it. example result_name2 or result_name.
#'
#' @examples
#' #result_name <- change_default_name('pca', SingleCellExperiment::reducedDimNames(sce))
#' @export
change_default_name <- function(result_name, name_list){
    ### Description ###
    # Changes the result_name to avoid overriding data

    # example usage
    # change_default_name('my_method', reduceDimNames(sce))

    if(result_name %in% name_list){

        # setting the name of the matrix if you have multiple default names
        method_n <- sum(grepl(result_name, name_list))
        new_name <-  paste0(result_name, method_n)
        warning(paste0('Chaning the result name to ', new_name))
        result_name <- new_name
    }
  return(result_name)
}
