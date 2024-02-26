#' Decomposition Approach : Dense-Matrix Principal Component Analysis
#'
#' This function performs Principal Component Analysis (PCA) by cohering
#' the matrix to be dense and using sklr.decomposition.PCA() with a
#' method of choice.
#'
#' @param sce <SingleCellExperiment object> SCE object
#' @param n_components <integer> number of Principal Components desired.
#' Must be n_components == min(n_cells, n_genes)
#' @param assay <character> specifying the assay to use, preferred raw counts!
#' @param method <character> by default 'randomized'.
#' If **auto** : The solver is selected by a default policy based on X.shape
#' and n_components: if the input data is larger than 500x500 and the
#' number of components to extract is lower than 80% of the smallest
#' dimension of the data, then the more efficient ‘randomized’
#' method is enabled. Otherwise the exact full SVD
#' is computed and optionally truncated afterwards.
#' If **full** : run exact full SVD calling the standard LAPACK solver
#' via scipy.linalg.svd and select the components by postprocessing
#' If **arpack** : run SVD truncated to n_components calling ARPACK
#' solver via scipy.sparse.linalg.svds. It requires strictly
#' 0 < n_components < min(X.shape).
#' If **randomized** : run randomized SVD by the method of Halko et al.
#' (from *sklearn.decomposition.PCA*)
#' @param center <bool> default TRUE; Whether to center in a column wise fashion
#' the assay prior to performing PCA.
#' @param scale <bool> default FALSE; Whether to scale in a column wise fashion
#' the assay prior to performing PCA.
#' @param result_name <character> default 'dense-pca';
#' Name used to store the result in the SingleCellExperiment object.
#' @param envname <character> default 'r-decomp';
#' Specify the name of the python virtual
#' environment to be used. If it does not exists it will create one and use it.
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
#' #decomp_dense_pca(sce, n_components = 50)
#' @export
decomp_dense_pca <- function(sce,
                             n_components=50,
                             assay='logcounts',
                             method='randomized',
                             center=TRUE,
                             scale=FALSE,
                             result_name='dense-pca',
                             envname='r-decomp',
                             return_model=FALSE,
                             seed=42,
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
    result_name <- change_default_name(result_name, SingleCellExperiment::reducedDimNames(sce))

    # cell view
    pca_h.matrix <- Matrix::t(pca.model$components_)
    # storing the percentage of variance explained
    sce <- store_H(sce        = sce,
                  h.matrix    = pca_h.matrix,
                  result_name = result_name,
                  latent_name = 'PC')

    # storing the %VarianceExplained by the PCA method
    metadata_name <- paste0(result_name,'_%VarianceExplained')
    S4Vectors::metadata(sce)[[metadata_name]] <- pca.model$explained_variance_ratio_



    # this is the return
    return_model(sce = sce, model = pca.model, return_model = return_model)
}

