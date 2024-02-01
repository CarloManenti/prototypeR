sce2adata_sparse <- function(sce,
                             envname = 'r-metacell',
                             main_layer = 'counts',
                             ...){

    ### Description ###
    # Converts a sce into an anndata, taking care also of the
    # conversion of the CSC to CSR sparse matrix.

    # example usage
    # sce2adata_sparse(sce)

    is_package_installed('sceasy')
    is_package_installed('reticulate')
    is_python_package_installed(envname = envname, packages.vec = c('scipy'))
    reticulate::use_virtualenv(envname)
    scipy <- reticulate::import('scipy')

    # the matrix will be already transposed! cells x features
    adata <- sceasy::convertFormat(sce,
                                   from = "sce",
                                   to = "anndata",
                                   main_layer = main_layer,
                                   drop_single_values = FALSE,
                                   ...)

    # Converting a CSC sparse matrix to CSR sparse matrix
    # to run faster the models
    if(scipy$sparse$csc$isspmatrix_csc(adata$X)){
      adata$X <- scipy$sparse$csr_matrix(adata$X)
    }

    return(adata)
}
