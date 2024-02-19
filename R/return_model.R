return_model <- function(sce, model, return_model, ...){
  ### Description ####
  # Handy function to return either only the SingleCellExperiment object
  # or (both) the SingleCellExperiment object AND the model used
  # during the analysis.


  # example usage
  # return_model(sce, model, return_model)

  if(!return_model) return(sce)
  return(list(obj = sce, model = model, ...))
}
