decomp_mraa <- function(sce,
                        result_name='mraa',
                        store_all_results=FALSE,
                        return_model=FALSE,
                        ...){
    ### Description ###
    # Performs MultiResolution Archetypal Analysis (mraa)
    # on a Single Cell Experiment Object

    # example usage
    # decomp_mraa(sce)

    message('--- Checking packages ---')
    is_package_installed('ACTIONet')
    is_package_installed('ACTIONetExperiment')
    is_package_installed('Matrix')

    library('ACTIONetExperiment')
    library('Matrix')

    ace <- ACTIONetExperiment::as.ACTIONetExperiment(sce)




    message('--- Running MultiResolution Archeytpal Analysis ---')
    #normalizing the data
    ace <- ACTIONet::normalize.ace(ace)
    # reducing the sparse log counts into a dense matrix of 50 PCs
    ace <- ACTIONet::reduce.ace(ace)
    # running Multi-Resolution Archetypal Anlysis
    ace <- ACTIONet::runACTIONet(ace, ...)



    message('--- Storing Results ---')
    # General settings to store results
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

    if(!return_model) return(sce)
    return(list(obj = sce, model = ace))
}









