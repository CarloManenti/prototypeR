#' Decomposition Approach : Sparse-Matrix Principal Component Analysis
#'
#' This function performs Principal Component Analysis (PCA) keeping a
#' sparse representation for the matrix and leveraging irlba::prcomp_irlba
#' implementation.
#'
#' @param sce <SingleCellExperiment object> SCE object
#' @param n_components <integer> number of Principal Components desired.
#' Must be n_components == min(n_cells, n_genes)
#' @param assay <character> default logcounts;
#' specifying the assay to use, preferred raw counts!
#' @param center <bool> default TRUE;
#' Whether to center in a column wise fashion
#' the assay prior to performing PCA.
#' @param scale <bool> default TRUE;
#' Whether to scale in a column wise fashion
#' the assay prior to performing PCA.
#' @param result_name <character> default 'pca';
#' Name of used to store the result in the SingleCellExperiment object.
#' @param return_model <bool> default FALSE; Whether to return also
#' the model and not only
#' the SingleCellExperiment object.
#' @param seed <integer> default 42; to set the seed for reproducibility.
#' @param verbose <bool> default FALSE; Whether to be prompted with message
#' for each step of the analysis.
#' @return either a SingleCellExperiment object with PCA representation for
#' only genes, or the SingleCellExperiment object and the model
#' used to perform PCA.
#' @examples
#' decomp_sparse_pca(sce, n_components = 50)
#' @export
decomp_sparse_pca <- function(sce,
                              n_components=50,
                              assay='logcounts',
                              center=TRUE,
                              scale=TRUE,
                              result_name='pca',
                              return_model=FALSE,
                              seed=42,
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
