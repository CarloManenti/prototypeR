decomp_nmf <- function(sce,
                       n_components,
                       method='Lsnmf',
                       assay='counts',
                       envname='r-decomp',
                       initialization='nndsvd',
                       result_name='nmf',
                       return_model=FALSE,
                       ...){
    ### Description ###
    # Performs various implementations of NFM via nimfa


    # example usage
    # decomp_nmf(sce, n_components = 14)



    message('--- Checking packages ---')
    is_package_installed('reticulate')
    is_python_package_installed(envname = envname, packages.vec = c('nimfa'))
    reticulate::use_virtualenv(envname)
    # to avoid environment releted warnings with nimfa and MATLAB
    warnings <- reticulate::import('warnings')
    warnings$simplefilter('ignore')
    nimfa    <- reticulate::import('nimfa')


    n_components <- as.integer(n_components)



    # we need to transpose the matrix to work with nimfa
    ## WARNING
    # the packages does not work with sparse matrix or it is super slow!
    # this should be fixed in next versions of my wrapper…
    message('--- Dense coercion :C ---')
    matrix.dense <- as.matrix(t(assay(sce, assay)))



    message('--- Performing ', method,  ' NMF ---')
    nfm.model <- switch(method,
                        'bayesian'  = nimfa$Bd(matrix.dense,
                                               rank = n_components,
                                               seed = initialization,
                                               ...),
                        'separable' = nimfa$methods$sepnmf(matrix.dense,
                                                           rank = n_components,
                                                           seed = initialization,
                                                           ...),
                        'default'  = nimfa$nmf(matrix.dense,
                                               rank = n_components,
                                               seed = initialization,
                                               ...),
                        'Lsnmf'    = nimfa$Lsnmf(matrix.dense,
                                                 rank = n_components,
                                                 seed = initialization,
                                                 ...))
    # fitting the model
    nmf.model = nfm.model()



    message('--- Storing Results ---')
    patterns_names <- paste0(rep('NNF_', n_components), seq_len(n_components))

    # Cell view
    nmf_h.dgCmatrix <- as(nmf.model$basis(), 'sparseMatrix')
    colnames(nmf_h.dgCmatrix) <- patterns_names
    rownames(nmf_h.dgCmatrix) <- colnames(sce)
    result_name <- change_default_name(result_name, reducedDimNames(sce))
    reducedDims(sce)[[result_name]] <- nmf_h.dgCmatrix

    # Gene view
    nmf_w.dgCmatrix <- as(t(nmf.model$coef()), 'sparseMatrix')
    colnames(nmf_w.dgCmatrix) <- patterns_names
    rownames(nmf_w.dgCmatrix) <- rownames(sce)
    metadata(sce)[[result_name]] <- nmf_w.dgCmatrix



    if(!return_model) return(sce)
    return(list(obj = sce, model = nmf.model))
}
