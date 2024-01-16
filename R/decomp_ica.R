decomp_ica <- function(sce,
                       n_components,
                       assay='counts',
                       method='fast', # fast, imax, jade
                       center=TRUE, # centering is generally requested!
                       scale=FALSE,
                       seed=42,
                       result_name='ica',
                       return_model=FALSE,
                       ...){
    ### Description ###
    # Computes ICA via three main implementations
    # InfoMax
    # JADE
    # FastICA


    # example usage
    #load(file = 'data/simulated_data.sce.rda')
    #sce <- simulated_data.sce
    #n_components <- 14



    message('--- Checking packages ---')
    is_package_installed('ica') # implements many variants of ICA

    # enforcing the type to avoid crashes with python functions
    n_components <- as.integer(n_components)
    seed <- as.integer(seed)



    message(paste0('--- Performing ', method ,' ICA ---'))
    assay.dgCMatrix <- assay(sce, assay)

    centered_assay.dgCMatrix <- center_and_scale(assay.dgCMatrix,
                                                 center = center,
                                                 scale = scale)

    ica.model <- ica::ica(X = t(centered_assay.dgCMatrix), # cells x features
                          nc = n_components, # number of Independent components
                          method = method,
                          center = FALSE,
                          ...) # other arguments which can be taken from ica



    message(paste0('--- Storing results ---'))

    # cell view
    ica_h.matrix <- ica.model$S
    colnames(ica_h.matrix) <- paste0(rep('IF_', n_components),
                                     seq_len(n_components))
    ica_h.dgCMatrix <- as(ica_h.matrix, 'sparseMatrix')
    result_name <- change_default_name(result_name, reducedDimNames(sce))
    reducedDims(sce)[[result_name]] <- ica_h.dgCMatrix

    # gene view
    ica_w.matrix <- ica.model$M
    colnames(ica_w.matrix) <- colnames(ica_h.matrix)
    rownames(ica_w.matrix) <- rownames(sce)
    ica_w.dgCMatrix <- as(ica_w.matrix, 'sparseMatrix')
    metadata(sce)[[result_name]] <- ica_w.dgCMatrix

    if(!return_model) return(sce)
    return(list(obj = sce, model = ica.model))
}

