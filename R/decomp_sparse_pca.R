decomp_sparse_pca <- function(sce,
                              n_components=50,
                              assay='logcounts',
                              center=TRUE,
                              scale=TRUE,
                              seed=42,
                              result_name='pca',
                              return_model=FALSE,
                              verbose=FALSE){
### Description ###
# Computes PCA on a sparse matrix specified stored in a specific assays
# of a SingleCellExperiment Object.

# example usage
# decomp_sparse_pca(sce, 100)


    if(verbose == TRUE){
      message('--- Checking packages ---')
    }
    is_package_installed('irlba') # fast and sparse PCA implementation in R
    # to make results reproducible
    set.seed(seed)

    if(verbose == TRUE){
      message('--- Performing sparse Implicitly Restarted Lanczos Bidiagonalization Algorithm (IRLBA) PCA ---')
    }
    assay.dgCMatrix <- assay(sce, assay)
    # performing IRLBA PCA on a sparse matrix
    pca.model <- irlba::prcomp_irlba(assay.dgCMatrix,
                                     n = n_components,
                                     center = center,
                                     scale. = scale)



    if(verbose == TRUE){
      message(paste0('--- Storing results ---'))
    }
    # just for this case, since we want to store also %VarianceExplained
    result_name <- change_default_name(result_name, reducedDimNames(sce))

    # gene view
    pca_w.matrix <- pca.model$x
    sce <- store_W(sce         = sce,
                   w.matrix    = pca_w.matrix,
                   result_name = result_name,
                   latent_name = 'PC')

    # cell view
    pca_h.matrix <- pca.model$rotation
    sce <- store_H(sce         = sce,
                   h.matrix    = pca_h.matrix,
                   result_name = result_name,
                   latent_name = 'PC')



    # storing the %VarianceExplained (pve) by the PCA method
    metadata_name <- paste0(result_name,'_%VarianceExplained')
    pve <- pca.model$sdev ^ 2 / sum (pca.model$sdev ^ 2)
    metadata(sce)[[metadata_name]] <- pve

    return_model(sce, model, return_model)
}
