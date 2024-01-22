decomp_aa <- function(sce,
                      n_components=14,
                      reduction_method='pca',  # dimensionality reduction method
                      reduced_representation=NULL,
                      n_dimensions=10,
                      assay='counts',          # assay to be reduced in dimensions
                      delta=0,                 #
                      conv_crit=1e-04,         #
                      seed=42,                 #
                      parallel_type='m',       # to run locally or on a single node
                      volume_ratio='t_ratio',  # use none for quicker and less precise computations
                      sample_proportion=NULL,  # 0.25 | with NULL it runs multiple times on the whole data set, to evaluate local optima.
                      envname='r-decomp',
                      method='pcha',
                      normalise_var=TRUE,
                      verbose=FALSE,
                      return_model=F,
                      boostrap_number=2,
                      result_name='aa'){
    ### Description ###
    # Performs Acrchetypal Analysis with bootstrap
    # bootstrap can be avoided!

    # example usage
    # decomp_aa(sce)



    message('--- Checking packages ---')
    is_package_installed('ParetoTI')

    n_components <- as.integer(n_components)
    n_dimensions <- as.integer(n_dimensions)
    conv_crit    <- as.double(conv_crit)
    seed         <- as.integer(seed)
    boostrap_number <- as.integer(boostrap_number)
    delta        <- as.double(delta)
    if(!is.null(sample_proportion)) as.double(sample_proportion)



    ## Note : Look at select_aas
    message(paste0('--- Reducing ', assay, ' matrix via ', reduction_method, ' ---'))
    sce <- switch(reduction_method,
                  'pca' = decomp_pca(sce,
                                     n_dimensions,
                                     assay = assay,
                                     result_name = reduction_method),
                  'nmf' = decomp_nmf(sce,
                                     n_dimensions,
                                     assay = assay,
                                     result_name = reduction_method),
                  'none' = sce)



    if(reduction_method == 'none'){
      reduced_matrix <- assay(sce, assay)
    }else{
      # we need to transpose the matrix to have it features x cells
      reduced_matrix <- t(reducedDims(sce)[[reduction_method]])}



    message(paste0('--- Computing archetypal analysis ---'))
    # Performing boostrap with the chosen k number of prototypes (archetypes)
    aa.model <- ParetoTI::fit_pch_bootstrap(data            = reduced_matrix, # matrix features x cells
                                            noc             = n_components,              # number of archetypes to find
                                            n               = boostrap_number,           # how many bootstraps do we want to do?
                                            sample_prop     = sample_proportion,         # fraction of the data set to consider when finding archetypes
                                            check_installed = T,                         # check if python module py_pcha is found.
                                            type            = parallel_type,             # do we want to parallelize? then use 'm' for multicore
                                            seed            = seed,                      # reproducibility seed
                                            delta           = delta,                     # inflates original polytope (simplex) fit such that it may contain more points of the dataset
                                            conv_crit       = conv_crit,                 # convergence criteria; more accuracy with 1e-6 but more computational time.
                                            method          = method,                    # pcha by default
                                            normalise_var   = normalise_var,             # normalize the variability of archetypes by the variance in the data
                                            verbose         = verbose,                   # whether you want more information during the run
                                            volume_ratio    = volume_ratio)              # this might be quite computational expensive for more than 8 dimensions



    message('--- Storing Results ---')
    patterns_names <- paste0(rep('A_', n_components), seq_len(n_components))

    # storing results from the best run among all!
    best_run <- as.integer(which(aa.model$pch_fits$SSE == min(aa.model$pch_fits$SSE)))

    # Cell view
    aa_h <- t(aa.model[['pch_fits']][['S']][[best_run.int]])
    aa_h.dgCmatrix <- as(aa_h, 'sparseMatrix')
    colnames(aa_h.dgCmatrix) <- patterns_names
    rownames(aa_h.dgCmatrix) <- colnames(sce)
    result_name <- change_default_name(result_name, reducedDimNames(sce))
    reducedDims(sce)[[result_name]] <- aa_h.dgCmatrix

    # Note: The Gene View can not be access since it returns C which has unexpected
    #       dimensions (n_archetypes x n_archetypes); instead of cells x archetypes

    # we should fix this later (even via a linear regression or whatever)



    if(!return_model) return(sce)
    return(list(obj = sce, model = aa.model))
}
