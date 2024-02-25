#' Decomposition Approach : Independent Component Analysis
#'
#' This function performs Independent Component Analysis (ICA) on a
#' SingleCellExperiment object and stores the results directly in it.
#' It is based on the ica package.
#'
#' @param sce <SingleCellExperiment object> SCE object
#' @param n_components <integer> number of Independent Component desired.
#' @param assay <character> assay used as input for ICA.
#' @param method <character> default **'fast'**; method used to perform
#' ICA. It can be one of (**"fast"**, **"imax"**, **"jade"**).
#' **'fast'** : Hyvarinen’s (1999) FastICA algorithm.
#' **'imax'** : Bell and Sejnowski’s (1995) Information-Maximization algorithm.
#' **'jade'** : Cardoso and Souloumiac’s (1993, 1996) Joint Approximate
#' Diagonalization of Eigenmatrices (JADE) algorithm.
#' @param center <bool> default TRUE; Whether to center the data or not;
#' centering is performed column wise, so on the cells (observations).
#' @param scale <bool> default FALSE; Whether to scale the data or not;
#' scaling is performed column wise, so on the cells (observations).
#' @param result_name <character> default 'ica';
#' Name of used to store the result in the SingleCellExperiment object.
#' @param return_model <bool> default FALSE; Whether to return also
#' the model and not only the SingleCellExperiment object.
#' @param seed <integer> default 42; to set the seed for reproducibility.
#' @param verbose <bool> default FALSE; Whether to be prompted with message
#' for each step of the analysis.
#' @return either a SingleCellExperiment object with PCA representation for
#' only genes, or the SingleCellExperiment object and the model
#' used to perform PCA.
#' @examples
#' decomp_ica(sce, 6, logcounts)
#' @export
decomp_ica <- function(sce,
                       n_components,
                       assay='logcounts',
                       method='fast', # fast, imax, jade
                       center=TRUE, # centering is generally requested!
                       scale=FALSE,
                       result_name='ica',
                       return_model=FALSE,
                       seed=42,
                       verbose=FALSE,
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
    #decomp_ica(sce, 6, 'logcounts')


    if(verbose == TRUE){
        message('--- Checking packages ---')
    }
    is_package_installed('ica') # implements many algorithms to perform ICA
    # enforcing the type
    n_components <- as.integer(n_components)
    seed <- as.integer(seed)



    if(verbose == TRUE){
        message(paste0('--- Performing ', method ,' ICA ---'))
    }
    # pre-processing the data
    assay.dgCMatrix <- assay(sce, assay)
    centered_assay.dgCMatrix <- center_and_scale(assay.dgCMatrix,
                                                 center = center,
                                                 scale = scale)
    # performing ICA
    ica.model <- ica::ica(X = t(centered_assay.dgCMatrix), # cells x features
                          nc = n_components, # number of Independent components
                          method = method,
                          center = FALSE, # center may be performed before…
                          ...) # other arguments which can be taken from ica



    if(verbose == TRUE){
        message('--- Storing results ---')
    }

    # gene view
    ica_w.matrix <- ica.model$M
    sce <- store_W(sce = sce,
                   w.matrix = ica_w.matrix,
                   latent_name = 'IC',
                   result_name = result_name)
    # cell view
    ica_h.matrix <- ica.model$S
    sce <- store_H(sce = sce,
                   h.matrix = ica_h.matrix,
                   result_name = result_name,
                   latent_name = 'IC')

    # this is the return
    return_model(sce = sce, model = ica.model, return_model = return_model)
}

