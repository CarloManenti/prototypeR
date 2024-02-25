#' Decomposition Approach : Kernel Archetypal Analysis
#'
#' This function performs Kernel Archetypal Analysis leveraging SEACells
#' implementation and store directly in the SingleCellExperiment object the
#' results.
#'
#' @param sce <SingleCellExperiment object> SCE object
#' @param n_components <integer or NULL> number of archetypes desired.
#' Can be set to NULL to fallback to the
#' suggested number of archetypes (n_cells / suggested_cells4aa).
#' @param n_highly_variable_genes <integer> default 5000; number of
#'  Highly Variable Genes to be selected and then used to find the archetypes.
#' note: setting the HVGs number too low will cause the method to fail,
#' empirically you want to use 1500 genes or more as HVGs.
#' For smaller data set one may use roughly half of the total number of genes.
#' @param n_pcs <integer> default 50; number of Principal Components
#' used for downstream analysis.
#' @param kernel_dimension <charter> default 'X_pca'; representation of
#' the data used to build the kernel (it refers to the AnnData object).
#' @param min_iter <integer> default 10; minimum number of iterations to be
#' performed before convergence of the AA algorithm. If convergence is reached
#' before the minimum number of iterations it will not be considered.
#' @param max_iter <integer> default 10000; maximum number of iterations to be
#' performed before convergence of the AA algorithm. Once reached the method
#' will simply fail.
#' @param suggested_cells4aa <integer> default 75; Number of suggested cells
#' for a single archetype. 75 is the number suggested by the authors of the
#' core method (SEACells). If n_components is specified, suggested_cells4aa will
#' be ignored.
#' @param n_waypoint_eigs <integer> default NULL; if NULL it will fallback
#' to the suggested number of way points, 15. Way points are the number of
#' eigen-values to consider when initializing archetypes.
#' note: the user should never lower the number of way points further 3!
#' @param convergence_epsilon <double> default 1e-5; Convergence criteria.
#' A lower convergence_epsilon will result in more accurate solutions, but
#' also in a longer computation time.
#' @param minimum_weight <double> default 0.05; minimal threshold to consider
#' the contribution of an archetype to a cell. If the contribution is lower, it
#' will be set to 0!
#' @param use_highly_variable_genes <bool> default TRUE; Whether to use only
#' Highly Variable Genes to compute Archetypes.
#' @param summarize_layer <character> default 'raw'; Starting data matrix used
#' for computing archetypes. It refers to the AnnData object.
#' @param low_dim_embedding <character> default 'X_pca'; low dimensional
#' representation used to compute the compactness.
#' @param plots <bool> default FALSE; whether to plot quality metrics plot.
#' @param bins <integer> default 5; Number of bins in the plot of
#' kernel archetypes contributions to cells.
#' @param celltype_label <character> default NULL; if specified it will compute
#' archetype purity based on the defined field of the object. Generally is
#' cell_type or cluster… . Used for plotting quality metrics.
#' @param nth_nearest_neighbor <integer> default 1;
#' computes the diffusion distance from its nth nearest neighbor, used for
#' plotting quality metrics.
#' @param max_n_comps_pca <integer> default 80; maximum number of Principal
#' Components computed while performing PCA and visualized in the elbow plot.
#' @param kernel_matrix_plot_dim <integer vectors> default seq(200). Number of
#' cells to plot in the kernel matrix (it will be a squared matri). Used for
#' plotting quality metrics.
#' @param figures_dir <character> default 'kaa_plots'. Name of the directory
#' used to store all the quality metrics plots. If it does not exists it will
#' be created.
#' @param elbow_plot_name <character> default '1_kaa_elbow_plot.pdf';
#' Name of the elbow plot saved in figures_dir.
#' @param kernel_matrix_plot_name <character> default
#' '2_portion_of_kernel_matrix.pdf';
#' Name of the kernel matrix plot saved in figures_dir.
#' @param initialization_plot_name <character> default
#' '3_initialization_plot.pdf';
#' Name of the initialization plot saved in figures_dir.
#' @param convergence_plot_name <character> default '4_convergence_plot.pdf';
#' Name of the convergence plot saved in figures_dir.
#' @param assignment_plot_name, <character> default '6_cell_assignment_plot.pdf'
#' Name of the cell assignment plot saved in figures_dir.
#' @param purity_plot_name <character> default '7_archetypes\'_purity_plot.pdf'
#' Name of the archetypes' purity plot saved in figures_dir.
#' @param compactness_plot_name default '8_archetypes\'_compactness_plot.pdf';
#' Name of the archetypes' compactness plot saved in figures_dir.
#' @param separation_plot_name default '9_archetypes\'_separation_plot.pdf';
#' Name of the archetypes' separation plot saved in figures_dir.
#' @param result_name <character> default 'kaa';
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
#' decomp_kaa(sce, n_components = 5, n_highly_variable_genes = 250, n_pcs = 10,
#' n_waypoint_eigs = 3, verbose = TRUE, plots = TRUE)
#' @export
decomp_kaa <- function(sce,
                       n_components,
                       n_highly_variable_genes=500,
                       n_pcs=50,
                       kernel_dimension='X_pca',
                       min_iter=10,
                       max_iter=10000,
                       suggested_cells4aa=75,
                       n_waypoint_eigs=NULL,
                       convergence_epsilon=1e-5,
                       minimum_weight=0.05,
                       use_highly_variable_genes=TRUE,
                       summarize_layer='raw',
                       low_dim_embedding='X_pca',
                       plots=FALSE,
                       bins=5,
                       celltype_label=NULL,
                       nth_nearest_neighbor=1,
                       max_n_comps_pca=80,
                       kernel_matrix_plot_dim=seq(200),
                       figures_dir='kaa_plots',
                       elbow_plot_name='1_kaa_elbow_plot.pdf',
                       kernel_matrix_plot_name='2_portion_of_kernel_matrix.pdf',
                       initialization_plot_name='3_initialization_plot.pdf',
                       convergence_plot_name='4_convergence_plot.pdf',
                       assignment_plot_name='6_cell_assignment_plot.pdf',
                       purity_plot_name='7_archetypes\'_purity_plot.pdf',
                       compactness_plot_name='8_archetypes\'_compactness_plot.pdf',
                       separation_plot_name='9_archetypes\'_separation_plot.pdf',
                       result_name='kaa',
                       envname='r-decomp',
                       return_model=FALSE,
                       seed=42,
                       verbose=TRUE
){

    if(verbose){
        message('--- Checking packages ---')
    }

    packages.vec <- c('SEACells',
                    'jupyter',
                    'ipywidgets',
                    'fastcluster',
                    'alabaster',
                    'anndata',
                    'cmake',
                    'Cython',
                    'h5py',
                    'joblib',
                    'kiwisolver',
                    'legacy-api-wrap',
                    'leidenalg',
                    'llvmlite',
                    'louvain',
                    'matplotlib',
                    'munkres',
                    'ncls',
                    'numba',
                    'numpy',
                    'palantir',
                    'pandas',
                    'psutil',
                    'pyranges',
                    'pyrle',
                    'scanpy',
                    'scikit-learn',
                    'scipy',
                    'seaborn',
                    'six',
                    'sorted-nearest',
                    'statsmodels',
                    'tabulate',
                    'tqdm',
                    'tzdata',
                    'tzlocal',
                    'umap-learn')

    is_python_package_installed(envname = envname,
                                packages.vec = packages.vec)
    #.rs.restartR() in case there are problems with the loading metacells
    # enforcing the use of the correct environment
    seacells <- reticulate::import('SEACells', delay_load = TRUE)      # kernel archetypal analysis
    sc <- reticulate::import('scanpy', delay_load = TRUE)              # pre-processing for kernela
    ad <- reticulate::import('anndata', delay_load = TRUE)             # passing from SCE to anndat
    plt <- reticulate::import('matplotlib.pyplot', delay_load = TRUE)  # visualizing results
    np <- reticulate::import('numpy', delay_load = TRUE)               # basic handling of data
    sns <- reticulate::import('seaborn', delay_load = TRUE)            # for quality plotting



    # enforcing types on variables
    if(!is.null(n_waypoint_eigs)){
        n_waypoint_eigs <- as.integer(n_waypoint_eigs)
    }

    if(!is.null(n_components)){
      n_components <- as.integer(n_components)
    }

    n_highly_variable_genes <- as.integer(n_highly_variable_genes)
    max_n_comps_pca <- as.integer(max_n_comps_pca)
    n_pcs <- as.integer(n_pcs)
    # selecting archetypes parameters
    min_iter <- as.integer(min_iter)
    max_iter <- as.integer(max_iter)
    # aggregating cells into Archetypes parameters
    seed <- as.integer(seed)
    # fixed parameters parameters
    suggested_cells4aa <- as.integer(suggested_cells4aa)
    convergence_epsilon <- as.double(convergence_epsilon)
    minimum_weight <- as.double(minimum_weight)
    bins <- as.integer(bins)
    nth_nearest_neighbor <- as.integer(nth_nearest_neighbor)
    figures_dir <- paste0('./', figures_dir)
    check_dir(figures_dir)
    figures_dir <- paste0(figures_dir, '/')



    if(verbose){
        message('--- Converting SCE into AnnData ---')
    }
    matrix.dgCMatrix <- counts(sce)
    data.h5ad <- ad$AnnData(t(matrix.dgCMatrix), dtype = 'float32')
    data.h5ad$obs_names <- colnames(matrix.dgCMatrix)
    data.h5ad$var_names <- rownames(matrix.dgCMatrix)
    data.h5ad$obs       <- reticulate::r_to_py(as.data.frame(colData(sce)))
    data.h5ad$obs$index <- colnames(sce)
    # enforcing unique feature names
    data.h5ad$var_names_make_unique()
    sc.obj <- data.h5ad
    # Counts must be stored in .raw
    # so we need to provide a AnnData object with col and row names.
    h5ad_raw.obj <- sc$AnnData(sc.obj['X'])
    # note: we need to pass a vector such that reticulate will handle it
    h5ad_raw.obj$obs_names <- np$array(sc.obj$obs_names)
    h5ad_raw.obj$var_names <- np$array(sc.obj$var_names)
    # storing the counts with row and col names in the scanpy object
    sc.obj$raw <- h5ad_raw.obj



    if(verbose){
        message('--- Preprocessing ---')
    }
    # From here onward we will work on the X matrix overwriting it at each step
    # It might be useful to store each layer individually…
    # 1. Normalizing
    sc$pp$normalize_total(sc.obj, inplace = TRUE)
    # normalize_cell is deprecated!
    # 2. Log conversion
    sc$pp$log1p(sc.obj)
    # 3. Selection of highly variable genes
    sc$pp$highly_variable_genes(sc.obj,
                                n_top_genes = n_highly_variable_genes)
    # 4. Computing PCA and selecting the main components
    sc$tl$pca(sc.obj,
              n_comps = max_n_comps_pca,
              use_highly_variable = use_highly_variable_genes)
    if(plots){
        # visualizing the variance explained by the PCs
        # saving the plot
        plot_name.pdf <- paste0(figures_dir, elbow_plot_name)  # .pdf file
        pdf(plot_name.pdf)
        sc$pl$pca_variance_ratio(sc.obj, n_pcs = max_n_comps_pca)
        dev.off()
    }
    # setting the number of PCs wanted by sub-setting
    # (this can not be performed from terminal)
    sc.obj$obsm['X_pca'] <- sc.obj$obsm['X_pca'][, seq(0, n_pcs)]



    if(verbose){
        message('--- Computing Kernel Based Archetypes --- \n')
    }
    # The developers of Kernel aa generally suggest 1 prototype each 75 cells
    if(is.null(n_components)){
        n_components <- round(sc.obj$shape[[1]] / suggested_cells4aa)
        n_components <- as.integer(n_components)
    }
    # Setting the number of way points
    if(is.null(n_waypoint_eigs)) n_waypoint_eigs = as.integer(15)
    # this is quite tricky since it will crash the code!
    # Furthermore, using more than 15 waypoints is a one way trip to crash town.
    # Ideally we need to test further when the code crash or not.
    # Would be even better an actual documentation…

    # Designing the model
    kernel_aa.model = seacells$core$SEACells(ad = sc.obj,
                                             build_kernel_on     = kernel_dimension,
                                             n_SEACells          = n_components,
                                             n_waypoint_eigs     = n_waypoint_eigs,
                                             convergence_epsilon = convergence_epsilon)

    # Performing the Kernel 'Trick'
    kernel_aa.model$construct_kernel_matrix()



    ## This will raise error if run from terminal
    if(plots){
        if(verbose){
           message('--- Plotting the Kernel Matrix ---')
        }
        is_package_installed('ComplexHeatmap')
        kernel_matrix.dgRMatrix <- kernel_aa.model$kernel_matrix
        kernel_matrix_plot <- kernel_matrix.dgRMatrix[kernel_matrix_plot_dim,
                                                      kernel_matrix_plot_dim]
        # making it dense
        kernel_matrix_plot <- as.matrix(kernel_matrix_plot)
        kernel_matrix.plot <- ComplexHeatmap::Heatmap(kernel_matrix_plot,
                                                      col = blues9)
        # saving the plot
        plot_name.pdf <- paste0(figures_dir, kernel_matrix_plot_name)
        pdf(file = plot_name.pdf)
        ComplexHeatmap::draw(kernel_matrix.plot)
        dev.off() # free memory closing plots!
    }

    # Computing Archetypes
    kernel_aa.model$initialize_archetypes()
    # this **could** break if the number of way points is higher
    # than the number of archetypes!

    kernel_aa.model$fit(min_iter = min_iter,
                        max_iter = max_iter)
    # this could break easily!
    # better put a very high threshold!!!
    # note: one could manually force the model to run
    # for more iterations by looping
    # over … model.step()



    if(verbose){
        message('--- Storing Results ---')
    }
    # General settings to store results
    ## Hard Assignment of cells to prototypes|archetypes
    # We need to remove 'SEACell-' preceding
    # the actual value of archetype!
    # then we rename the vector,
    # otherwise we will lose the names
    # with the  conversion to integer.
    membership.matrix <- kernel_aa.model$get_hard_assignments()
    # getting the actual assignment chr
    membership.vec <- membership.matrix[, 1]
    # removing SEACell-
    membership.vec <- as.integer(gsub(".*-","", membership.vec))
    # assigning cell names to the membership vector
    names(membership.vec) <- rownames(membership.matrix)
    colData(sce)[[result_name]] <- membership.vec

    # gene view
    kernel_aa_W.obj <- seacells$core$summarize_by_soft_SEACell(sc.obj,
                                                               np$matrix(kernel_aa.model$A_),
                                                               summarize_layer = summarize_layer,
                                                               minimum_weight  = minimum_weight)
    ## storing  raw counts
    kernel_aa_W.obj$layers['counts'] <- kernel_aa_W.obj$X
    ## then we can normalize and log transform W_raw
    ## overwriting X
    sc$pp$normalize_total(kernel_aa_W.obj, inplace = TRUE)
    ## storing log-normalized counts in log1p
    kernel_aa_W.obj$layers['log1p'] <- sc$pp$log1p(kernel_aa_W.obj['X'],
                                                   copy = TRUE)
    # finally storing the gene view
    kaa_w.matrix <- t(kernel_aa_W.obj$layers['log1p'])
    sce <- store_W(sce = sce,
                   w.matrix = kaa_w.matrix,
                   latent_name = 'KA',
                   result_name = result_name)

    # cell view
    kaa_h.matrix <- t(kernel_aa.model$A_)
    sce <- store_H(sce = sce,
                   h.matrix = kaa_h.matrix,
                   result_name = result_name,
                   latent_name = 'KA')



    if(plots){
      if(verbose){
          message('--- Quality metrics Plots ---')
      }
      # Computing a UMAP for visualization purposes
      sc$tl$umap(sc.obj)
      # saving the initialization plot…
      # Are my initial archetypes representative of my data?
      # are we happy?
      plot_name.pdf <- paste0(figures_dir, initialization_plot_name)
      seacells$plot$plot_initialization(ad = sc.obj,
                                        model = kernel_aa.model,
                                        save_as = plot_name.pdf)
      # saving the convergence plot
      # Does the model follow a smooth convergence curve or not?
      # are we happy?
      plot_name.pdf <- paste0(figures_dir, convergence_plot_name)
      kernel_aa.model$plot_convergence(save_as = plot_name.pdf)

      plot_name.pdf <- paste0(figures_dir,
                              '5_cell_assignment_plot.pdf')

      seacells$plot$plot_2D(sc.obj,
                            key ='X_umap',
                            colour_metacells = TRUE,
                            save_as = plot_name.pdf)

      # Cell Contribution to Archetypes
      ## it uses a deprecated function of seabron!
      plot_name.pdf <- paste0(figures_dir,
                              assignment_plot_name)

      seacells$plot$plot_SEACell_sizes(sc.obj,
                                       bins = bins,
                                       title = 'Distribution of Kernel based Archetypes Sizes',
                                       save_as = plot_name.pdf)
      # Archetypes' Purity
      if(!is.null(celltype_label)){
          purity.vec <- seacells$evaluate$compute_celltype_purity(sc.obj,
                                                                  col_name = celltype_label)[, 2]

          # plotting purity across all archetypes
          purity_plot <- sns$violinplot(purity.vec)
          plt$title('Archetypes\' Purity')
          plt$gca()$set_xticklabels(list()) # to set to null the x ticks
          purity_plot$axes$set(xlabel = NULL)

          # saving the plot
          plt$show()
          plot_name.pdf <- paste0(figures_dir,
                                  purity_plot_name)
          plt$savefig(plot_name.pdf)
          plt$close() # remember to close plots to free memory!
      }

      # Archetypes' Compactness
      # Computes the variance in diffusion components in a selected compressed
      # dimensional representation. (double check that…)
      # Lower values of compactness suggest more compact/lower
      # variance archetypes
      compactness.vec <- seacells$evaluate$compactness(sc.obj,
                                                       low_dim_embedding = low_dim_embedding)

      # plotting compactness across all archetypes
      compactness_plot <- sns$violinplot(compactness.vec)
      plt$title(paste0('Archetypes\' Compactness in ', low_dim_embedding, '\n reduced space'))
      plt$gca()$set_xticklabels(list()) # to set to null the x ticks
      compactness_plot$axes$set(xlabel = NULL)

      # saving the plot
      plt$show()
      plot_name.pdf <- paste0(figures_dir,
                              compactness_plot_name)
      plt$savefig(plot_name.pdf)
      plt$close() # remember to close plots to free memory!

      # Archetypes' Separation
      ## For each cell it computes the diffusion distance between itself and its
      ## nth_nearest_neighbor in a reduced space (like PCA space).
      ## I do not know how diffusion is performed, double check it!
      separation.vec <- seacells$evaluate$separation(ad = sc.obj,
                                                     low_dim_embedding = low_dim_embedding,
                                                     nth_nbr = nth_nearest_neighbor)
      ## Specifying the cluster variables as 'cell_types' could be interestingly
      ## It would force to compute the distances between cells belonging
      ##to the same cell type… . But does not seem to work as expected!

      ## plotting separation across all archetypes
      separation_plot <- sns$violinplot(separation.vec)
      plt$title(paste0('Archetypes\' Separation in ',low_dim_embedding ,'\n reduced space'))
      plt$gca()$set_xticklabels(list()) # to set to null the x ticks
      separation_plot$axes$set(xlabel = NULL)

      ## saving the plot
      plt$show()
      plot_name.pdf <- paste0(figures_dir,
                              separation_plot_name)

      plt$savefig(plot_name.pdf)
      plt$close()
    }



    # this is the return
    return_model(sce = sce, model = kernel_aa.model, return_model = return_model)
}




