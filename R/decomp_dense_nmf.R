#' Decomposition Approach : Dense matrix NMF
#'
#' This function performs NMF with a dense
#' coercion of the input matrix. One can choose between different NMF
#' implementations and initialization. The results are directly stored in
#' the input SingleCellExperiment (SCE) object.
#'
#' @param sce <SingleCellExperiment object> SCE object
#' @param n_components <integer> number of Non Negative Factors desired.
#' @param assay <character> default 'logcounts'; name of the assay to be used
#' during the downstream analysis.
#' @param method <character> default 'Lsnmf'; Type of algorithm used to perform
#' NMF. It can be one of…
#' *Lsnmf* : Alternating Non Negative Least Squares Matrix Factorization Using
#' Projected Gradient (bound constrained optimization) method for each
#' sub-problem (LSNMF) (Lin2007). (From nimfa)
#' *bayesian* : Bayesian Decomposition (BD) - Bayesian Non Negative Matrix
#' Factorization Gibbs sampler (bayesian) (Schmidt2009). (From nimfa)
#' *separable*: Separable Nonnegative Matrix Factorization (SepNMF)
#' (Damle2014), (Benson2014), (Kumar2013), (Gillis2014), (Tepper2015),
#' (Kapralov2016). (From nimfa)
#' *default*: Standard NMF with Euclidean /
#' Kullback-Leibler update equations and Frobenius /
#' divergence / connectivity cost functions (Lee2001), (Brunet2004).
#' (from nimfa)
#' @param initialization <character> default nndsvd; Type of initialization used
#' for performing NMF. Can be one of:
#' *nndsvd* : Non Negative Double Singular Value Decomposition
#' (NNDSVD) (Boutsidis2007) (from nimfa).
#' *randomc* Random C (Albright2006) is an inexpensive initialization
#' method for NMF. (From nimfa).
#' *randomv* Random V (Albright2006) is an inexpensive initialization
#' method for NMF. (From nimfa).
#' *random* Random is the simplest MF initialization method. (From nimfa)
#' @param ignore_warnings <bool> default FALSE; Whether to avoid warnings.
#' note: using verbose = FALSE and ignore warnings = TRUE may result in
#' the suppression of warnings due to suppressing most of the messages from
#' the function.
#' @param result_name <character> default 'dense-nmf';
#' Name of used to store the result in the SingleCellExperiment object.
#' @param envname <character> default 'r-decomp';
#' Specify the name of the python virtual
#' environment to be used. If it does not exists it will create one and use it.
#' @param return_model <bool> default FALSE; Whether to return also
#' the model and not only
#' the SingleCellExperiment object.
#' @param seed <integer> default 42; to set the seed for reproducibility.
#' @param verbose <bool> default FALSE; Whether to be prompted with message
#' for each step of the analysis.
#' @param ... <extra arguments for the specific implementation of nmf used>
#' @return either a SingleCellExperiment object with both H and W
#' representations or the SingleCellExperiment object and the model
#' used to perform NMF.
#' @examples
#' library('packageX')
#' data(sce)
#' decomp_dense_nmf(sce, n_components = 3)
#' decomp_dense_nmf(sce, n_components = 3, method = 'bayesian')
#' decomp_dense_nmf(sce, n_components = 3, method = 'separable')
#' decomp_dense_nmf(sce, n_components = 3, method = 'default')
#' @export
decomp_dense_nmf <- function(sce,
                       n_components,
                       assay='logcounts',
                       method='Lsnmf',
                       initialization='nndsvd',
                       ignore_warnings=FALSE,
                       result_name='dense-nmf',
                       envname='r-decomp',
                       return_model=FALSE,
                       seed=42,
                       verbose=FALSE,
                       ...){
    ### Description ###
    # Performs various implementations of NFM via nimfa
    # example usage
    # decomp_nmf(sce, n_components = 14)


    if(verbose){
        message('--- Checking Packages ---')
    }
    is_package_installed('reticulate')
    is_python_package_installed(envname = envname,
                                packages.vec = c('nimfa'))
    # to avoid environment related warnings with nimfa and MATLAB
    ignore_warnings(ignore_warnings = ignore_warnings,
                    envname = envname,
                    verbose = verbose)
    # delayed laod for CRAN assement
    nimfa    <- reticulate::import('nimfa', delay_load = TRUE)
    # enforcing the type of variables to avoid python crashesh
    n_components <- as.integer(n_components)



    # we need to transpose the matrix to work with nimfa
    ## WARNING
    # the packages does not work with sparse matrix or it is super slow!
    # this should be fixed in next versions of my wrapper…
    if(verbose){
        warning('--- Dense coercion ---')
    }
    matrix.dense <- invisible(as.matrix(Matrix::t(assay(sce, assay))))



    if(verbose){
       message('--- Performing ', method,  ' NMF ---')
    }
    nfm.model <- switch(method,
                        'bayesian'  = nimfa$Bd(matrix.dense,
                                               rank = n_components,
                                               seed = initialization,
                                               ...),
                        'separable' = nimfa$methods$sepnmf(matrix.dense,
                                                           rank = n_components,
                                                           seed = initialization,
                                                           ...),
                        'default'  = nimfa$methods$nmf(matrix.dense,
                                               rank = n_components,
                                               seed = initialization,
                                               ...),
                        'Lsnmf'    = nimfa$Lsnmf(matrix.dense,
                                                 rank = n_components,
                                                 seed = initialization,
                                                 ...))
    # fitting the model
    nmf.model <- nfm.model()



    if(verbose){
        message('--- Storing Results ---')
    }
    # gene view
    nmf_w.matrix <- Matrix::t(nmf.model$coef())
    sce <- store_W(sce = sce,
                   w.matrix = nmf_w.matrix,
                   result_name = result_name,
                   latent_name = 'NNF')
    # cell view
    nmf_h.matrix<- nmf.model$basis()
    sce <- store_H(sce = sce,
                   h.matrix = nmf_h.matrix,
                   result_name = result_name,
                   latent_name = 'NNF')



    # this is the return
    return_model(sce = sce, model = nmf.model, return_model = return_model)
}
