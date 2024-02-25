#' Decomposition Approach : MultiResolution Archetypal Analysis
#'
#' This function performs MultiResolution Archetypal Analysis as implemented by
#' the ACTIONet framework. It stores directly in the SingleCellExperiment object
#' the results. It can work in parallel using multiple cores.
#'
#' @param sce <SingleCellExperiment object> SCE object
#' @param n_components <integer or NULL> default NULL; number of archetypes
#' desired. If specified it will force a single resolution, the number of
#' specified archetypes. Otherwise, it will work at multiple resolution and will
#' prune useless or repeated archetypes found at multiple resolution, keeping
#' only informative ones.
#' @param norm.assay <character or NULL> default NULL; name of the assay which
#' stores the normalized counts. If NULL it will be obtained via shifted
#' logarithm normalization.
#' @param reduced.assay <character or NULL> default NULL; name of the assay
#' which stores the dimensionality reduction representation of H
#' (cells x latents). If NULL it will be obtained via randomized PCA.
#' @param store_all_results <bool> default FALSE; Whether to store all the
#' intermediate results obtained via ACTIONet in the SingleCellExperiment
#' object.
#' @param result_name <character> default 'mraa';
#' Name used to store the result in the SingleCellExperiment object.
#' @param return_model <bool> default FALSE; Whether to return also
#' the model and not only
#' the SingleCellExperiment object.
#' @param seed <integer> default 42; to set the seed for reproducibility.
#' @param verbose <bool> default FALSE; Whether to be prompted with message
#' for each step of the analysis.
#' @return either a SingleCellExperiment object with H and W representation
#' obtained via MultiResolution Archetypal Analysis, or the SingleCellExperiment
#' and the ACTIONet object used to perform the various steps
#' @examples
#' decomp_mraa(sce, n_components = 14, verbose = T)
#' @export

decomp_mraa <- function(sce,
                        n_components=NULL,
                        norm.assay=NULL,
                        reduced.assay=NULL,
                        store_all_results=FALSE,
                        result_name='mraa',
                        return_model=FALSE,
                        seed=42,
                        verbose=FALSE,
                        ...){
    ### Description ###
    # Performs MultiResolution Archetypal Analysis (mraa)
    # on a Single Cell Experiment Object

    # example usage
    # decomp_mraa(sce)

    if(verbose){
        message('--- Checking packages ---')
    }
    is_package_installed('ACTIONet')
    is_package_installed('ACTIONetExperiment')
    is_package_installed('Matrix')
    library('ACTIONet')
    set.seed(seed)

    ace <- ACTIONetExperiment::as.ACTIONetExperiment(sce)



    #normalizing the data only if they are not already normalized
    if(is.null(norm.assay)){
        if(verbose){
          message('--- Normalizing Raw Counts ---')
        }
        ace <- ACTIONet::normalize.ace(ace,
                                       assay_name = "counts",
                                       assay_out = "logcounts")
        norm.assay <- 'logcounts'
    }

    # reducing the sparse log counts into a dense matrix of 50 PCs
    if(is.null(reduced.assay)){
        if(verbose){
          message('--- Reducing Dimensions ---')
        }
        ace <- ACTIONet::reduce.ace(ace)
    }

    if(verbose){
        message('--- Running MultiResolution Archeytpal Analysis ---')
    }

    # running Multi-Resolution Archetypal Anlysis
    if(is.null(n_components)){
        ace <- ACTIONet::runACTIONet(ace,
                                     assay_name = norm.assay,
                                     reduction_slot = reduced.assay,
                                     ...)
    }else{ # forcing single resolution
        ace <- ACTIONet::runACTIONet(ace,
                                     k_min = n_components,
                                     k_max = n_components,
                                     assay_name = norm.assay,
                                     reduction_slot = reduced.assay,
                                     ...)
    }



    if(verbose){
        message('--- Storing Results ---')
    }
    # General settings to store results without store_H and store_W
    result_name <- change_default_name(result_name, reducedDimNames(sce))
    n_mraa <- ncol(colMaps(ace)[['H_unified']])
    patterns_names <- paste0(rep('MRAA_', n_mraa), seq_len(n_mraa))

    # setting the names
    colnames(colMaps(ace)[['H_unified']]) <- patterns_names
    colnames(rowMaps(ace)[['unified_feature_profile']]) <- patterns_names

    colMaps(ace)[[result_name]] <- colMaps(ace)[['H_unified']]
    metadata(ace)[[result_name]] <- rowMaps(ace)[['unified_feature_profile']]

    if(store_all_results==TRUE){
        sce <- ACTIONetExperiment::as.SingleCellExperiment(ace)
    }else{
        reducedDims(sce)[[result_name]] <- colMaps(ace)[[result_name]]
        metadata(sce)[[result_name]] <- metadata(ace)[[result_name]]
    }



    # this is the return
    return_model(sce = sce, model = ace, return_model = return_model)
}









