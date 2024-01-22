decomp_pca <- function(sce,
                       n_components,
                       assay='counts',
                       method='randomized',
                       center=TRUE,
                       scale=FALSE,
                       seed=42,
                       envname='r-decomp',
                       result_name='pca',
                       return_model=FALSE){
    ### Description ###
    # Computes pca on a single cell object assay
    # note : by setting the method you can select the type of pca
    #        you want to perform. For big matrix is best to use randomized
    # note : since it is a svd based pca, it is required to at least center
    #        the data before performing the pca.

    # example usage
    # load(file='data/simulated_data.sce.rda')
    # simulated_data.sce <- decomp_pca(simulated_data.sce, 12)



    message('--- Checking packages ---')
    is_package_installed('reticulate')

    is_python_package_installed(packages.vec = 'scikit-learn',
                                envname = envname)

    reticulate::use_virtualenv(envname)

    # importing 'only' the decomposition approaches of sklrn
    sklr.decomposition <- reticulate::import('sklearn.decomposition')

    # enforcing the type to avoid crashes with python functions
    n_components <- as.integer(n_components)

    seed <- as.integer(seed)



    message(paste0('--- Performing ', method, ' PCA ---'))
    assay.dgCMatrix <- assay(sce, assay)

    # feature-wise centering and scaling
    scaled_assay.matrix <- center_and_scale(matrix.dgCMatrix = assay.dgCMatrix,
                                            center = center,
                                            scale = scale)

    pca.model <- sklr.decomposition$PCA(n_components = n_components,
                                              svd_solver   = method,
                                              random_state = seed)

    # it requires a dense matrix to fit the model
    pca.model$fit(as.matrix(scaled_assay.matrix))




    message(paste0('--- Storing results ---'))

    # cell view
    pca_h.matrix <- t(pca.model$components_)

    colnames(pca_h.matrix) <- paste0(rep('PC_', n_components),
                                     seq_len(n_components))

    rownames(pca_h.matrix) <- colnames(sce)

    pca_h.dgCMatrix <- as(pca_h.matrix, 'sparseMatrix')

    result_name <- change_default_name(result_name, reducedDimNames(sce))

    reducedDims(sce)[[result_name]] <- pca_h.dgCMatrix

    # gene view not present for pca

    if(!return_model) return(sce)
    return(list(obj = sce, model = pca.model))
}






