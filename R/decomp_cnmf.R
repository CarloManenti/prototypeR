#' Decomposition Approach : Consensus Non-Negative Matrix Factorization
#'
#' This function performs Consensus Non-Negative Matrix Factorization (cNMF)
#' on an assay of a SingleCellExperiment and stores the results directly in
#' the object.
#'
#' @param sce <SingleCellExperiment object> SCE object
#' @param n_components <integer> number of Non-Negative Factors (NNF) desired.
#' @param levels <vector of integers> sequence of number of components to
#' investigate. It will run cNMF for each of these levels. The n_components
#' must fall inside this interval.
#' @param assay <character> default counts;
#' specifying the assay to use, preferred raw counts!
#' @param num_iterations <integer> default 100;
#' number of iterations used to build
#' the consensus matrices
#' @param num_workers <integer> default 0;
#' number of parallel processes; if 0 it will detect the number of cores
#' available in the machine and will set it to the max_number of cores - 2.
#' Generally, it is best to set it yourself given the dimension of the data set.
#' @param num_hvgenes <integer> default 5000;
#' number of Highly Variable Genes to be selected
#' prior to performing cNMF.
#' @param return_TPM <bool> default FALSE;
#' Whether to return the Transcript Per Million or
#' the actual gene weights. Generally gene weights highlights better difference
#' in the data.
#' @param density_threshold <double> default 0.01;
#' distance threshold (euclidean distance) of
#' a NNF to its K Nearest Neighbors to be consider for building the
#' consensus matrices.
#' @param run_name <character> default Intermediate_results;
#' name of the specific run of cNMF, it will be also
#' the name of the dedicated folder with the intermediate results.
#' @param output_dir <character> default ./cnmf;
#' directory used to store intermediate results.
#' @param ignore_warnings <bool> default TRUE; Whether to ignore warnings
#' @param file_name <character> default 'Intermediate_anndata.h5ad';
#' Name of the file in which is stored the AnnData object processed by cNMF.
#' To pass the data to the cNMF core function, one must save the AnnData
#' object in a .h5ad file. file_name is the name of the .h5ad file.
#' @param result_name <character> default cnmf;
#' name of used to store the result in the
#' SingleCellExperiment object.
#' @param envname <character> default r-decomp;
#' specify the name of the python virtual
#' environment to be used. If it does not exists it will create one and use it.
#' @param return_model <bool> default FALSE;
#' Whether to return also the model and not only
#' the SingleCellExperiment object.
#' @param seed <integer> default 42; to set the seed for reproducibility.
#' @param verbose <bool> default FALSE;
#' Whether to be prompted with message for each step of the analysis.
#' @return either a SingleCellExperiment object with cNMF representation for
#' both genes and cells, or the SingleCellExperiment object and the model
#' used to perform cNMF.
#'
#' @examples
#' library('packageX')
#' data(sce)
#' decomp_cnmf(sce = sce, n_components = 3, levels = 3 , num_iterations = 50, density_threshold = 1)
#' @export
decomp_cnmf <- function(sce,
                        n_components,
                        levels,
                        assay='counts',
                        num_iterations=100,
                        num_workers=4,
                        num_hvgenes=5000,
                        return_TPM=FALSE,
                        density_threshold=0.01,
                        run_name='Intermediate_results',
                        output_dir='~/cnmf',
                        result_name='cnmf',
                        envname='r-decomp',
                        return_model=FALSE,
                        seed=42,
                        verbose=FALSE,
                        ignore_warnings=TRUE,
                        file_name='Intermediate_anndata.h5ad'){
  ### Description ###
  # Performs cNMF
  # Be sure that levels contains the number of components you want
  # Also, take a look at the figures for picking the density threshold!
  # note : It could be nice to suppress warnings!
  # example usage
  # decomp_cnmf(sce)

  # EXT - split this core functions into smaller ones to make use of previous
  #       runs of cNMF.


  if(ignore_warnings){
    # must be fixed in future implementations… :)
    # removing warnings
    warnings <- reticulate::import('warnings')
    warnings$filterwarnings("ignore")

    # to avoid silly warnings from cNMF
    warnings$simplefilter("ignore")

  }

  if(verbose == TRUE){
    message('--- Checking packages ---')
  }
  is_package_installed('reticulate')  # python in r
  is_package_installed('sceasy')      # sce -> <- anndata
  is_python_package_installed(envname = envname,
                              packages.vec = c('cnmf', 'scanpy'))

  cnmf <- reticulate::import('cnmf', delay_load = TRUE)
  sc   <- reticulate::import('scanpy',  delay_load = TRUE)

  # enforcing the type of variable to respect python constraints
  num_iterations <- as.integer(num_iterations)
  num_workers    <- as.integer(num_workers)
  num_hvgenes    <- as.integer(num_hvgenes)
  n_components   <- as.integer(n_components)
  seed           <- as.integer(seed)

  # Setting the levels…
  # if we want a given number of patterns
  if(methods::is(levels, 'numeric')){
    levels <- as.integer(levels)
  }
  # since it does not make sense to extract a single latent (for now)
  if(min(levels) < 2) levels <- seq(2, max(levels))
  density_threshold <- as.double(density_threshold)



  if(verbose == TRUE){
    message('--- Writing intermediate AnnData ---')
  }
  # cNMF needs a file with an AnnData obj to run…
  adata <- sc$AnnData(X = Matrix::t(SummarizedExperiment::assay(sce, assay)))
  file_name <- paste0('./', file_name)
  sc$write(file_name,   adata)



  if(verbose == TRUE){
    message('--- Performing cNMF ---')
  }
  # initializing the cNMF object
  check_dir(paste0(output_dir, '/', run_name))
  cnmf_obj = cnmf$cNMF(output_dir = output_dir,
                       name       = run_name)

  # pre-processing the data to retain only the number of over-dispersed genes
  # that we want
  cnmf_obj$prepare(counts_fn         = file_name,
                   components        = levels,
                   n_iter            = num_iterations,
                   seed              = seed,
                   num_highvar_genes = num_hvgenes)

  # performing parallel cNMF
  cnmf_obj$factorize_multi_process(total_workers = num_workers)

  # defining the consensus
  cnmf_obj$combine()

  # Stability vs Loss
  cnmf_obj$k_selection_plot(close_fig = FALSE)
  cnmf_obj$paths['k_selection_plot']

  # Coherence of the decomposition level
  sapply(levels, function(level_i){ # which decomposition level
    cnmf_obj$consensus(k = level_i,
                       density_threshold = density_threshold,
                       show_clustering = TRUE,
                       close_clustergram_fig = TRUE)
  })

  # cnmf.model data fields…
  # 1 : usage_norm
  # 2 : gep_scores
  # 3 : gep_tpm
  # 4 : topgenes
  cnmf.model <- cnmf_obj$load_results(K = n_components,
                                      density_threshold = density_threshold)



  if(verbose == TRUE){
    message('--- Storing Results ---')
  }
  # gene view
  if(return_TPM){
    data_field <- 3
  }else{
    data_field <- 2
  }
  cnmf_w <- cnmf.model[[data_field]]
  sce <- store_W(sce = sce,
                 w.matrix = cnmf_w,
                 latent_name = 'NNF',
                 result_name = result_name)

  # cell view
  cnmf_h <- cnmf.model[[1]]
  sce <- store_H(sce = sce,
                 h.matrix = cnmf_h,
                 result_name = result_name,
                 latent_name = 'NNF')



  # this is the return
  return_model(sce = sce, model = cnmf.model, return_model = return_model)
}
