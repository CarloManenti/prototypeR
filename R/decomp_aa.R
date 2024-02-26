#' Decomposition Approach : Archetypal Analysis
#'
#' This function take has input a SingleCellExperiment object and performs
#' Archetypal Analysis (AA) leveraging (based on ParetoTI) on a
#' specified compressed representation.One can perform AA on the whole data
#' set, or on subsets of it via bootstrapping to account for possible outliers.
#'
#' @param sce <SingleCellExperiment object> with a your data
#' @param n_components <integer> number of archetypes desired
#' @param reduction_method <character> default *pca*;
#'  one of : *'pca'* (sparse pca), *'nmf'* (dense nmf) or *'none'*
#'  specifying the type of dimensionality reduction to be performed prior to AA.
#' @param assay <character> default 'logcounts'; assay used as input for the
#' reduction_method. If reduction_method = 'none, it will be the input for AA.
#' Ignored if reduced_representation is specified
#' @param n_dimensions <integer> default 50; number of dimensions used
#' to perform the reduction_method; ignored if reduced_representation
#' is specified.
#' @param reduced_representation <character> default NULL; It specifies an
#' already existing compressed representation to be used prior to AA.
#' If specified, reduction_method, n_dimensions and assay will be ignored.
#' @param delta <integer> default 0;  Inflates the original polytope (simplex)
#' fit such that it may contain more points of the data set. (from ParetoTI)
#' @param conv_crit <double> default 1e-04; convergence criteria; a lower value
#' (like 1e-6) will result in higher accuracy costing a longer computational
#' time. (from ParetoTI)
#' @param parallel_type <character> default *m*; one of *s*, *m*, *cmq*.
#' Do we want to parallelize? then use *m* for multicore.
#' *s* means single core processing using lapply.
#' *m* means multi-core parallel procession using parLapply.
#' *cmq* means multi-node parallel processing on a computing cluster
#' using clustermq package. (from ParetoTI)
#' @param volume_ratio <character> default *t_ratio*;
#' *t_ratio* : defines the volume of the convex hull of the data and the
#' t-ratio.
#' *variance_ratio* : defined the ratio of  variance of archetype positions to
#' variance of data; .
#' *none* :  do nothing; .
#' Caution! Geometric figure should be at least simplex:
#' qhull algorithm will fail to find convex hull of flat 2D shapes
#' in 3D, 3D shapes in 4D and so on. So, for calculating this dimensions
#' are selected based order of rows in data. Makes sense for principal
#' components but not for original data. Caution 2!
#' Computation time and memory use increase very quickly with dimensions.
#' Do not use for more than 7-8 dimensions. (from ParetoTI)
#' @param sample_proportion <double or NULL> default NULL; either NULL or the
#' proportion of data set that should be included in each sample.
#' If NULL the polytope fitting algorithm is run n times on
#' the same data which is useful for evaluating how often the
#' algorithm gets stuck in local optima. (from ParetoTI)
#' @param method <character> default *pcha*; method for archetypal analysis:
#' *pcha* : Principal Convex Hull,
#' *kmeans* : KMeans
#' *louvain* : Louvain (should be avoided)
#' @param normalise_var <bool> default TRUE; either TRUE or FALSE,
#' whether to normalise variance in position of archetypes by variance
#' in data in each dimensions.(from ParetoTI).
#' @param boostrap_number <integer or NA> default 20; number of bootstrap
#' samples on random data to measure variability in archetype positions.
#' Set to NA fit once to all data instead of bootstraping.
#' When this option is positive seed bootstrap_seed and
#' sample_prop must be provided. (from ParetoTI)
#' @param result_name <character> default aa; name of used to store the
#' result in the SingleCellExperiment object.
#' @param envname <character> default r-decomp; specifies the name of the
#' python virtual environment to be used. If it does not exists it will
#' create one and use it.
#' @param return_model <bool> default FALSE; Whether to return also the model
#' and not only the SingleCellExperiment object.
#' @param seed <integer> default 42; to set the seed for reproducibility.
#' @param verbose <bool> default FALSE; Whether to be prompted with message
#' for each step of the analysis.
#' @return either a SingleCellExperiment object with AA representation for both
#' genes and cells, or the SingleCellExperiment object and the model used to
#' perform archetypal analysis.
#' @examples
#' library('packageX')
#' data(sce)
#' decomp_aa(sce = sce, n_components = 5, reduced_representation = 'pca')
#' @export
decomp_aa <- function(sce,
                      n_components,
                      reduction_method='pca',
                      assay='logcounts',
                      n_dimensions=50,
                      reduced_representation=NULL,
                      delta=0,
                      conv_crit=1e-04,
                      parallel_type='m',
                      volume_ratio='t_ratio',
                      sample_proportion=NULL,
                      method='pcha',
                      normalise_var=TRUE,
                      boostrap_number=20,
                      result_name='aa',
                      envname='r-decomp',
                      return_model=FALSE,
                      seed=42,
                      verbose=FALSE){
    ### Description ###
    # Performs Acrchetypal Analysis with bootstrap
    # bootstrap can be avoided!

    # example usage
    # decomp_aa(sce)


    if(verbose == TRUE){
        message('--- Checking packages ---')
    }

    is_python_package_installed(packages.vec = c('py-pcha', 'geosketch'), envname = envname)
    is_package_installed('ParetoTI')

    n_components    <- as.integer(n_components)
    n_dimensions    <- as.integer(n_dimensions)
    conv_crit       <- as.double(conv_crit)
    seed            <- as.integer(seed)
    boostrap_number <- as.integer(boostrap_number)
    delta           <- as.double(delta)

    if(!is.null(sample_proportion)){
        as.double(sample_proportion)
    }
    ## Note : Look at select_aas
    if(verbose == TRUE){
        message(paste0('--- Reducing ', assay, ' matrix via ', reduction_method, ' ---'))
    }

    # checking if the reduced dimension representation is already present
    if(!is.null(reduced_representation)){
      if(reduced_representation %in% SingleCellExperiment::reducedDimNames(sce)){
        reduced_matrix <- Matrix::t(SingleCellExperiment::reducedDim(sce, reduced_representation))
      }else{
        stop('Please provide a reduced representation present in the SingleCellExperiment object')
      }
    }else{

        sce <- switch(reduction_method,
                      'pca' = decomp_sparse_pca(sce,
                                                n_dimensions,
                                                assay = assay,
                                                result_name = reduction_method),
                      'nmf' = decomp_dense_nmf(sce,
                                               n_dimensions,
                                               assay = assay,
                                               result_name = reduction_method),
                      'none' = sce)

        # we need to transpose the matrix to have it features x cells
        reduced_matrix <- Matrix::t(SingleCellExperiment::reducedDim(sce, reduction_method))
    }

    if(reduction_method == 'none'){
       reduced_matrix <- SummarizedExperiment::assay(sce, assay)
    }




    if(verbose == TRUE){
        message(paste0('--- Computing archetypal analysis ---'))
    }
    # Performing boostrap with the chosen k number of prototypes (archetypes)
    aa.model <- ParetoTI::fit_pch_bootstrap(data            = reduced_matrix, # matrix features x cells
                                            noc             = n_components,              # number of archetypes to find
                                            n               = boostrap_number,           # how many bootstraps do we want to do?
                                            sample_prop     = sample_proportion,         # fraction of the data set to consider when finding archetypes
                                            check_installed = TRUE,                      # check if python module py_pcha is found.
                                            type            = parallel_type,
                                            seed            = seed,                      # reproducibility seed
                                            delta           = delta,                     # inflates original polytope (simplex) fit such that it may contain more points of the dataset
                                            conv_crit       = conv_crit,                 # convergence criteria; more accuracy with 1e-6 but more computational time.
                                            method          = method,                    # pcha by default
                                            normalise_var   = normalise_var,             # normalize the variability of archetypes by the variance in the data
                                            verbose         = verbose,                   # whether you want more information during the run
                                            volume_ratio    = volume_ratio)              # this might be quite computational expensive for more than 8 dimensions



    if(verbose == TRUE){
        message('--- Storing Results ---')
    }

    if(is.null(sample_proportion)){
        # storing results from the best run among all!
        best_run <- as.integer(which(aa.model$pch_fits$SSE == min(aa.model$pch_fits$SSE)))
        # cell view
        aa_h <- Matrix::t(aa.model[['pch_fits']][['S']][[best_run]])
        sce <- store_H(sce, h.matrix = aa_h, result_name = result_name, latent_name = 'A')

        # Note: The Gene View can not be access since it returns C which has unexpected
        #       dimensions (n_archetypes x n_archetypes); instead of cells x archetypes
        # we should fix this later (even via a linear regression or whatever)
    }else{
        # Storing the reduced H given the sample_proportion specified
        # still looking only at the best iteration out of all with respect to the
        # SSE.
        result_name <- change_default_name(result_name = result_name, name_list = SingleCellExperiment::reducedDimNames((sce)))
        n_latents <- ncol(aa_h)
        colnames(aa_h) <- paste0(rep('A', n_latents), seq_len(n_latents))
        # rownames are already conserved and correct given that there only the
        # amount defined by sample_proportion.
        h.dgCMatrix <- methods::as(aa_h, 'sparseMatrix')
        S4Vectors::metadata(sce)[[result_name]] <- h.dgCMatrix
    }



    # this is the return
    return_model(sce = sce, model =  aa.model, return_model = return_model)
}
