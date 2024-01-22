select_aas <- function(sce,
                       n_dimensions=10,         # number of dimensions for the compressed representation
                       reduction_method='pca',  # dimensionality reduction method
                       assay='counts',          # assay to be reduced in dimensions
                       max_k=3,                 # maximum k number
                       maxiter=10,              # maximum number of iterations to find the archetypes in a given run
                       boostrap_number=5,       #
                       delta=0,                 #
                       conv_crit=1e-04,         #
                       seed=42,                 #
                       bootstrap=TRUE,          # whether to perform bootstrap for more robust results
                       parallel_type='m',       # to run locally or on a single node
                       volume_ratio='t_ratio',  # use none for quicker and less precise computations
                       order_type="align",      #
                       sample_proportion=NULL,  # 0.25 | with NULL it runs multiple times on the whole data set, to evaluate local optima.
                       interval_of_k=NULL,
                       compressed_representation=NULL,
                       envname='r-decomp',
                       point_size=0.5,          # plotting parameters
                       line_size=0.5,
                       text_axis_size=5,
                       plot_title_size=10,
                       hjust=0.5,
                       return_model=F){
  ### Description ###
  # Runs Boostraped AA on multiple ks (number of archetypes)
  # to pick the best number of AAs for the dataset
  # It also plots a few quality metrics to select the number of AAs.

  # example usage
  # select_aas(sce)

  # note: the default values are just for a toy example!
  #       change them depending on the data set and the purpose
  #       of the analysis.




    message('--- Checking packages ---')
    is_package_installed('ParetoTI')
    is_package_installed('reticulate')
    is_package_installed('patchwork')

    is_python_package_installed(envname = 'r-decomp',
                                packages.vec = c('py_pcha', 'datetime',
                                                 'scipy', 'numpy'))
    reticulate::use_virtualenv(envname)



    # Bootstrapping parameters
    n_dimensions         = as.integer(n_dimensions)           # number of components
    max_k                = as.integer(max_k)           # maximum k number
    maxiter              = as.integer(maxiter)         # maximum number of iterations to find the archetypes in a given run
    boostrap_number      = as.integer(boostrap_number) # number of trial via bootstrap
    delta                = as.integer(delta)
    conv_crit            = as.double(conv_crit)        # convergence threshold
    seed                 = as.integer(seed)

    # plotting variables
    point_size      = as.double(point_size)
    line_size       = as.double(line_size)
    text_axis_size  = as.integer(text_axis_size)
    plot_title_size = as.integer(plot_title_size)
    hjust           = as.double(hjust)

    # setting the range of archetypes to test on the data
    if(is.null(interval_of_k)) interval_of_k <- seq(2, max_k)
    if(max_k <= 2) warning('Please use a k higher then 2 get a comparison')


    # Setting the compressed representation

    ## ADD ALSO ICA!
    ## test for nmf and none

    ## also add an option to specify where to pick the decomp matrix!
    ## instead of running it again
    ## make this a function on its own!
    message(paste0('--- Reducing ', assay, ' matrix via ', reduction_method, ' ---'))
    sce <- switch(reduction_method,
                  'pca' = decomp_pca(sce,
                                     n_dimensions,
                                     assay = assay,
                                     result_name = reduction_method),
                  'nmf' = decomp_nmf(sce,
                                     n_dimensions,
                                     assay = assay,
                                     result_name = reduction_method),
                  'none' = sce)


    if(reduction_method == 'none'){
        reduced_matrix <- assay(sce, assay)
    }else{
        # we need to transpose the matrix to have it features x cells
        reduced_matrix <- t(reducedDims(sce)[[reduction_method]])}


    message('--- Bootstraping Archetypes ---')
    aas.model = ParetoTI::k_fit_pch(data            = reduced_matrix,
                                    ks              = interval_of_k,
                                    check_installed = TRUE,
                                    bootstrap       = bootstrap,
                                    bootstrap_N     = boostrap_number,
                                    maxiter         = maxiter,
                                    bootstrap_type  = parallel_type,
                                    seed            = seed,
                                    volume_ratio    = volume_ratio,
                                    delta           = delta,
                                    conv_crit       = conv_crit,
                                    order_type      = order_type,
                                    sample_prop     = sample_proportion)



    ## Quality plots to pick the number of prototypes (archetypes / latents)
    # Percentage of Variance Explained Plot
    message('--- Plotting Quality Metrics ---')
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

    library(patchwork)
    quality_plot <- variance_explained_plot /
                    variance_gain_plot      /
                    total_variance_plot     /
                    t_ratio_plot

    if(!return_model) return(quality_plot)
    return(list(plot = quality_plot, model = aas.model))
}







