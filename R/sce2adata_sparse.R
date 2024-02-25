#' Handy Function to Convert a SingleCellExperiment Object into AnnData Object
#'
#' This function converts a SingleCellExperiment (SCE) object into AnnData
#' object and it keeps most of the fields intact! It is also very useful in
#' converting the matrix types in a proper way.
#'
#' @param sce <SingleCellExperiment object> SCE object
#' @param result_name <character> default 'X';
#' Name used to store the result in the SingleCellExperiment object.
#' @param envname <character> default 'r-decomp';
#' Specify the name of the python virtual
#' environment to be used. If it does not exists it will create one and use it.
#' @param main_layer <character> default counts. Main layer to be converted.
#' @return an AnnData object with most of the fields of the original SCE
#' @examples
#' sce2adata_sparse(sce)
#' @export
sce2adata_sparse <- function(sce,
                             envname='r-decomp',
                             main_layer='counts',
                             ...){
    ### Description ###
    # Converts a sce into an anndata, taking care also of the
    # coversion of the CSC to CSR sparse matrix.

    # example usage
    # sce2adata_sparse(sce)

    is_package_installed('sceasy')
    is_package_installed('reticulate')
    is_python_package_installed(envname = envname, packages.vec = 'scipy')
    reticulate::use_virtualenv(envname)
    scipy <- reticulate::import('scipy')

    # the matrix will be already transposed! cells x features
    adata <- sceasy::convertFormat(sce,
                                   from="sce",
                                   to="anndata",
                                   main_layer=main_layer,
                                   drop_single_values=FALSE,
                                   ...)

    # Converting a CSC sparse matrix to CSR sparse matrix
    # to run fuster the model
    if(scipy$sparse$csc$isspmatrix_csc(adata$X)){
        adata$X <- scipy$sparse$csr_matrix(adata$X)
    }

    return(adata)
}
