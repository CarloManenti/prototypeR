#' Handy Function to Select an Appropriate Number of Archetypes
#'
#' This function take has input a SingleCellExperiment object and performs
#' Archetypal Analysis (AA) leveraging (based on ParetoTI) on a
#' specified compressed representation. The core idea of this function is
#' to help the user to select an appropriate number of archetypes for the
#' downstream analysis. It sinergises with decomp_aa.
#' note: one can perform AA on the whole data set,
#'  or on subsets of it via bootstrapping to account for possible outliers.
#'
#' @param sce <SingleCellExperiment object> with a your data
#' @param interval_of_k <vector of integers or NULL> default NULL;
#' Sequence of number of latents to try on the data. Generally it never includes
#' 1 in it. If NULL, the sequence will range form 2 up to max_k.
#' note : If both max_k and interval_of_k are NULL, the function will fail!
#' @param max_k <integer> default NULL; Maximum number of latents to find in
#' the data. If specified and the interval_of_k is NULL, it will define the
#' higher limit of the sequence of latents to be found in the data. This
#' range starts from 2 and goes up to max_k (included).
#' If NULL it will use the specified interval_of_k.
#' note : If both max_k and interval_of_k are NULL, the function will fail!
#' @param reduction_method <character> default *pca*;
#'  one of : *'pca'* (sparse pca), *'nmf'* (dense nmf) or *'none'*
#'  specifying the type of dimensionality reduction to be performed prior to AA.
#' @param assay <character> default 'logcounts'; assay used as input for the
#' reduction_method. If reduction_method = 'none, it will be the input for AA.
#' Ignored if reduced_representation is specified.
#' @param n_dimensions <integer> default 50; number of dimensions used
#' to perform the reduction_method; ignored if reduced_representation
#' is specified.
#' @param reduced_representation <character> default NULL; It specifies an
#' already existing compressed representation to be used prior to AA.
#' If specified, reduction_method, n_dimensions and assay will be ignored.
#' @param max_iter <integer> default 10; Number of runs to used to find a given
#' number of latents.
#' note: Total amount of runs will be :
#' max_iter ° bootstrap_n ° length(interval_of_k).
#' @param boostrap_n <integer or NA> default 20; number of bootstrap
#' samples on random data to measure variability in archetype positions.
#' Set to NA fit once to all data instead of bootstraping.
#' When this option is positive seed bootstrap_seed and
#' sample_prop must be provided. (from ParetoTI)
#' @param bootstrap <bool> default TRUE; Whether to perform bootstrap.
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
#' @param order_type <character> default align; how to allign latents found in
#' multiple runs. One of:
#' *align*  : aligns the latents by matching a reference shape with the one they
#' define.
#' *cosine* : aligns the latents computing the cosine similarity between them
#' *side*   : aligns the latents matching the sides of the shape they define.
#' @param sample_proportion <double or NULL> default NULL; either NULL or the
#' proportion of data set that should be included in each sample.
#' If NULL the polytope fitting algorithm is run n times on
#' the same data which is useful for evaluating how often the
#' algorithm gets stuck in local optima. (from ParetoTI)
#' @param point_size <double> default 0.5; dimensions of the plotted point.
#' for quality metrics plots
#' @param line_size <double> default 0.5; width of the plotted lines.
#' for quality metrics plots
#' @param text_axis_size <integer or double> default 10; dimension of the text
#' on the axis. for quality metrics plots
#' @param plot_title_size <integer or double> default 15; dimension of
#' the title. for quality metrics plots
#' @param hjust <double> default 0.5; Position of the title (by default in the
#' central). for quality metrics plots
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
#' #select_aas(sce, interval_of_k = seq(3, 4), boostrap_n = 2, reduced_representation = 'pca', verbose = TRUE)
#' @export
select_aas <- function(sce,
                       interval_of_k=NULL,
                       max_k=NULL,
                       reduction_method='pca',
                       assay='logcounts',
                       n_dimensions=50,
                       reduced_representation=NULL,
                       max_iter=10,
                       boostrap_n=20,
                       bootstrap=TRUE,
                       delta=0,
                       conv_crit=1e-04,
                       parallel_type='m',
                       volume_ratio='t_ratio',
                       order_type="align",
                       sample_proportion=NULL,
                       point_size=0.5,
                       line_size=0.5,
                       text_axis_size=5,
                       plot_title_size=10,
                       hjust=0.5,
                       envname='r-decomp',
                       return_model=FALSE,
                       seed=42,
                       verbose=FALSE){
    ### Description ###
    # Runs Boostraped AA on multiple ks (number of archetypes)
    # to pick the best number of AAs for the dataset
    # It also plots a few quality metrics to select the number of AAs.

    # example usage
    # select_aas(sce, seq(2, 5))

    # EXT - silence reticulate virtual env install


    if(verbose){
        message('--- Checking packages ---')
    }
    is_package_installed('ParetoTI')
    is_package_installed('reticulate')
    # installing useful python packages for AA analysis with ParetoTI
    is_python_package_installed(envname = 'r-decomp',
                                packages.vec = c('py_pcha',
                                                 'datetime',
                                                 'scipy',
                                                 'numpy'))

    # enforcing bootstrapping parameters
    n_dimensions         <- as.integer(n_dimensions)
    max_k                <- as.integer(max_k)
    max_iter             <- as.integer(max_iter)
    boostrap_n           <- as.integer(boostrap_n)
    delta                <- as.integer(delta)
    conv_crit            <- as.double(conv_crit)
    seed                 <- as.integer(seed)
    # enforcing plotting variables
    point_size      <- as.double(point_size)
    line_size       <- as.double(line_size)
    text_axis_size  <- as.integer(text_axis_size)
    plot_title_size <- as.integer(plot_title_size)
    hjust           <- as.double(hjust)

    # setting the range of archetypes to test on the data
    if(is.null(interval_of_k) & is.null(max_k)){
          stop('Please provide either interval_of_k or max_k')
    }
    # falling back to max_k
    if(is.null(interval_of_k) | length(interval_of_k) < 2){
        interval_of_k <- seq(2, max_k)
        if(max_k <= 2){
             stop('Please use a max_k higher then 2 to get a comparison across resolutions!')
        }
    }



    # cheking if the reduced dimension representation is already present
    if(!is.null(reduced_representation)){
        if(reduced_representation %in% SingleCellExperiment::reducedDimNames(sce)){
            reduced_matrix <- Matrix::t(SingleCellExperiment::reducedDim(sce, reduced_representation))
        }else{
          stop('Please provide a reduced representation present in the SingleCellExperiment object')
        }
    }else{
        # performing dimensionality reduction
        # since there is not one already specified
        if(verbose){
            message(paste0('--- Reducing ', assay, ' matrix via ', reduction_method, ' ---'))
        }
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
    }
    # using the raw assay if we do not want to perform dimensionality reduction
    if(reduction_method == 'none'){
        reduced_matrix <- assay(sce, assay)
    }else{
        # we need to transpose the matrix to have it features x cells
        reduced_matrix <- Matrix::t(SingleCellExperiment::reducedDims(sce)[[reduction_method]])
    }



    if(verbose){
        message('--- Bootstraping Archetypes Across Resolutions ---')
    }
    aas.model = ParetoTI::k_fit_pch(data            = reduced_matrix,
                                    ks              = interval_of_k,
                                    check_installed = TRUE,
                                    bootstrap       = bootstrap,
                                    bootstrap_N     = boostrap_n,
                                    maxiter         = max_iter,
                                    bootstrap_type  = parallel_type,
                                    seed            = seed,
                                    volume_ratio    = volume_ratio,
                                    delta           = delta,
                                    conv_crit       = conv_crit,
                                    order_type      = order_type,
                                    sample_prop     = sample_proportion)



    ## Quality plots to pick the number of prototypes (archetypes / latents)
    # Percentage of Variance Explained Plot
    if(verbose){
        message('--- Plotting Quality Metrics ---')
    }
    variance_explained_plot <- plot_aas_feature(aas.model,
                                                "varexpl",
                                                'Percentage of Variance Explained',
                                                point_size = point_size,
                                                line_size = line_size,
                                                text_axis_size = text_axis_size,
                                                plot_title_size = plot_title_size,
                                                hjust = hjust)

    variance_gain_plot <- plot_aas_feature(aas.model,
                                           'res_varexpl',
                                           'Variance Gain',
                                           point_size = point_size,
                                           line_size = line_size,
                                           text_axis_size = text_axis_size,
                                           plot_title_size = plot_title_size,
                                           hjust = hjust)

    total_variance_plot <- plot_aas_feature(aas.model,
                                            'total_var',
                                            'Archetypes Variability',
                                            point_size = point_size,
                                            line_size = line_size,
                                            text_axis_size = text_axis_size,
                                            plot_title_size = plot_title_size,
                                            hjust = hjust)

    t_ratio_plot <- plot_aas_feature(aas.model,
                                     't_ratio',
                                     'T Ratio',
                                     point_size = point_size,
                                     line_size = line_size,
                                     text_axis_size = text_axis_size,
                                     plot_title_size = plot_title_size,
                                     hjust = hjust)

    quality_plot <- list(variance_explained = variance_explained_plot,
                         variance_gain = variance_gain_plot,
                         total_variance = total_variance_plot,
                         t_ratio = t_ratio_plot)



    # slightly ugly return but it is more suitable than return_model…
    if(!return_model) return(quality_plot)
    return(list(plot = quality_plot, model = aas.model))
}







