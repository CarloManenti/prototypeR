decomp_cnmf <- function(sce,
                        assay='counts',
                        n_components,         # number of prototypes to use
                        levels=seq(13, 15),      # range of prototypes to be extracted
                        ## running multiple times NMF and defining a consensus
                        num_iterations=100,      # number of iterations used to build the consensus
                        num_workers=4,           # number of parallel processes
                        num_hvgenes=5000,        # number of over-dispersed genes to consider
                        ## Quality metrics plots : Coherence of the decomposition level
                        density_threshold=0.01, # threshold to filter out prototypes|latents (max euclidean distance between replicates)
                        ## Folders to store the results of the runs
                        run_name='Intermediate_results', # name of the specific run
                        output_dir='./cnmf',
                        envname='r-cnmf',
                        file_name='Intermediate_anndata.h5ad', # stored in data
                        seed=42,
                        result_name='cnmf',
                        return_model=FALSE,
                        ...){
    ### Description ###
    # Performs cNMF
    # Be sure that levels contains the number of components you want
    # Also, take a look at the figures for picking the density threshold!

    # note : It could be nice to suppress warnings!

    # example usage
    # decomp_cnmf(sce)



    message('--- Checking packages ---')
    is_package_installed('reticulate')
    is_python_package_installed(envname = envname, packages.vec = c('cnmf', 'scanpy'))
    reticulate::use_virtualenv(envname)

    # to avoid silly warnings from cNMF
    warnings <- reticulate::import('warnings')
    warnings$simplefilter('ignore')

    cnmf <- reticulate::import('cnmf')
    sc <- reticulate::import('scanpy')

    num_iterations <- as.integer(num_iterations)
    num_workers    <- as.integer(num_workers)
    num_hvgenes    <- as.integer(num_hvgenes)
    n_components <- as.integer(n_components)
    seed <- as.integer(seed)
    # since it does not make sense to extract a single latent (for now)
    if(min(levels) < 2) levels <- seq(2, max(levels))
    density_threshold <- as.double(density_threshold)



    message('--- Writing intermediate AnnData ---')
    data.h5ad <- sc$AnnData(X = t(assay(sce, assay)))
    file_name <- paste0('./', file_name)
    sc$write(file_name,   data.h5ad)



    message('--- Performing cNMF ---')
    # initializing the cNMF object
    cnmf_obj = cnmf$cNMF(output_dir = output_dir,
                         name = run_name)

    # pre-processing the data to retain only the number of over-dispersed genes
    # that we want
    cnmf_obj$prepare(counts_fn = file_name,
                     components = levels,
                     n_iter = num_iterations,
                     seed = seed,
                     num_highvar_genes = num_hvgenes)

    # performing parallel cNMF
    cnmf_obj$factorize_multi_process(total_workers = num_workers)

    # defining the consensus
    cnmf_obj$combine()

    # Stability vs Loss
    cnmf_obj$k_selection_plot(close_fig=FALSE)
    cnmf_obj$paths['k_selection_plot']

    # Coherence of the decomposition level
    sapply(levels, function(level_i){ # which decomposition level we want to look at
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



    message('--- Storing Results ---')
    patterns_names <- paste0(rep('NNF_', n_components), seq_len(n_components))

    # Cell view
    cnmf_h.dgCMatrix <- as(cnmf.model[[1]], 'sparseMatrix')
    colnames(cnmf_h.dgCMatrix) <- patterns_names
    rownames(cnmf_h.dgCMatrix) <- colnames(sce)
    result_name <- change_default_name(result_name, reducedDimNames(sce))
    reducedDims(sce)[[result_name]] <- cnmf_h.dgCMatrix

    # Gene view
    cnmf_w.dgCMatrix <- as(cnmf.model[[2]], 'sparseMatrix')
    colnames(cnmf_w.dgCMatrix) <- patterns_names
    rownames(cnmf_w.dgCMatrix) <- rownames(sce)
    metadata(sce)[[result_name]] <- cnmf_w.dgCMatrix

    if(!return_model) return(sce)
    return(list(obj = sce, model = cnmf.model))
}



