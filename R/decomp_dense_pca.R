decomp_dense_pca <- function(sce,
                       n_components=50,
                       assay='logcounts',
                       method='randomized',
                       center=TRUE,
                       scale=FALSE,
                       seed=42,
                       envname='r-decomp',
                       result_name='pca',
                       return_model=FALSE,
                       verbose=FALSE){
    ### Description ###
    # Computes pca on a single cell object assay !converts it to a dense matrix!
    # note : by setting the method you can select the type of pca
    #        you want to perform. For big matrix is best to use randomized
    # note : since it is a svd based pca, it is required to at least center
    #        the data before performing the pca.

    # example usage
    # load(file='data/simulated_data.sce.rda')
    # simulated_data.sce <- decomp_pca(simulated_data.sce, 12)



    if(verbose == TRUE){
        message('--- Checking packages ---')
    }
    is_package_installed('reticulate')
    is_python_package_installed(packages.vec = 'scikit-learn',
                                envname = envname)
    # importing 'only' the decomposition approaches of sklrn
    sklr.decomposition <- reticulate::import('sklearn.decomposition',
                                             delay_load = TRUE)
    # enforcing the type to avoid crashes with python functions
    n_components <- as.integer(n_components)
    seed <- as.integer(seed)



    if(verbose == TRUE){
        message(paste0('--- Performing dense ', method, ' PCA ---'))
    }
    assay.dgCMatrix <- assay(sce, assay)
    # feature-wise centering and scaling
    scaled_assay.matrix <- center_and_scale(matrix.dgCMatrix = assay.dgCMatrix,
                                            center = center,
                                            scale  = scale)
    # defining the model parameters
    pca.model <- sklr.decomposition$PCA(n_components = n_components,
                                              svd_solver   = method,
                                              random_state = seed)

    # it requires a dense matrix to fit the model
    pca.model$fit(as.matrix(scaled_assay.matrix)) # double check this passage




    if(verbose == TRUE){
        message(paste0('--- Storing results ---'))
    }
    # just for this case, since we want to store also %VarianceExplained
    result_name <- change_default_name(result_name, reducedDimNames(sce))

    # cell view
    pca_h.matrix <- t(pca.model$components_)
    # storing the percentage of variance explained
    sce <- store_H(sce        = sce,
                  h.matrix    = pca_h.matrix,
                  result_name = result_name,
                  latent_name = 'PC')

    # storing the %VarianceExplained by the PCA method
    metadata_name <- paste0(result_name,'_%VarianceExplained')
    metadata(sce)[[metadata_name]] <- pca.model$explained_variance_ratio_



    if(!return_model) return(sce)
    return(list(obj = sce, model = pca.model))
}

