#' Handy Function to Return Multiple Objects
#'
#' This function lets you return 2 different object, the SingleCellExperiment
#' (SCE) object and the model used during the analysis. It may return only the
#' SCE object or both the SCE and the model given the return_model boolean
#' variable.
#'
#' @param sce <SingleCellExperiment object> SCE object
#' @param model <object of various types> model used in a
#' partition or decomposition method.
#' @param return_model <bool> default FALSE; Whether to return also
#' the model and not only the SingleCellExperiment object.
#' @param ... <extra arguments to be returned, must be specified as
#' new_argumnet_name = new_argument>
#' @return either a SCE object or the SCE object and the model used to perform
#' the decompostion or partition method
#' @examples
#' #return_model(sce = c(1, 2 ,3), model = c('a', 'b', 'c'), return_model = TRUE, bool = c(TRUE, FALSE))
#' @export

return_model <- function(sce, model, return_model, ...){
  ### Description ####
  # Handy function to return either only the SingleCellExperiment object
  # or (both) the SingleCellExperiment object AND the model used
  # during the analysis.

  # example usage
  # return_model(sce, model, return_model)



  if(!return_model){
      return(sce)
  }
  return(list(obj = sce, model = model, ...))
}
