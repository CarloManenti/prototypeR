decomp_kaa <- function(sce,
                       envname='r-reticulate',
                       result_name='kaa',
                       n_highly_variable_genes=500,
                       # setting the HVGs number too low will crash the code 100%
                       # better stick to more tha 1500 genes
                       # For smaller dataset use roughly half of the total number of genes.
                       ## PCA and selection of PCs
                       use_highly_variable_genes=TRUE,
                       max_n_comps_pca=10,
                       n_comps_pca=10,
                       ## Applying the Kernel
                       n_components=5, # if null it follow default setting of n_cells / 75
                       kernel_dimension='X_pca', # representation of the data used to build the kernel
                       ## Selecting Archetypes
                       min_iter=10,
                       max_iter=min_iter *1000,
                       ## Aggregating cells into Archetypes
                       celltype_label='cluster',
                       seed=42,
                       figures_dir='kaa_plots',
                       return_model=FALSE,
                       ## Fixed Parameters
                       suggested_cells4aa=75, # suggested number of cells for each archetype.
                       n_waypoint_eigs=3, # number of eigen-values to consider when initializing archetypes. Falling back to deafult with NULL.
                       convergence_epsilon=1e-5,  # stopping criteria.
                       kernel_matrix_plot_dim=seq(200),  # dimension of the kernel matrix plot.
                       summarize_layer='raw',# type of data to aggregate.
                       minimum_weight=0.05,  # minimal threshold to consider the contribution of an archetype.
                       bins=5,  # for the plot of the distribution of cell contributions to kernel archetypes
                       low_dim_embedding='X_pca', # low dimensional representation used to compute the compactness
                       nth_nearest_neighbor=1,  # computes the diffusion distance from its first neighbour!
                       plots=FALSE,
                       elbow_plot_name='1_kaa_elbow_plot.pdf',
                       kernel_matrix_plot_name='2_portion_of_kernel_matrix.pdf',
                       initialization_plot_name='3_initialization_plot.pdf',
                       convergence_plot_name='4_convergence_plot.pdf',
                       assignment_plot_name='6_cell_assignment_plot.pdf',
                       purity_plot_name='7_arcehtypes\'_purity_plot.pdf',
                       compactness_plot_name='8_arcehtypes\'_compactness_plot.pdf',
                       separation_plot_name='9_arcehtypes\'_separation_plot.pdf'
){

    message('--- Checking packages ---')

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
                    #'tables',
                    'tabulate',
                    'tqdm',
                    'tzdata',
                    'tzlocal',
                    'umap-learn')





    is_python_package_installed(envname = envname,
                                packages.vec = packages.vec)
    #.rs.restartR() in case there are problems with the loading metacells
    # enforcing the use of the correct environment
    reticulate::use_virtualenv(envname)

    seacells <- reticulate::import('SEACells')      # kernel archetypal analysis
    sc <- reticulate::import('scanpy')              # pre-processing for kernela
    ad <- reticulate::import('anndata')             # passing from SCE to anndat
    plt <- reticulate::import('matplotlib.pyplot')  # visualizing results
    np <- reticulate::import('numpy')               # basic handling of data
    sns <- reticulate::import('seaborn')            # for quality plotting



    ## Enforcing types on variables
    if(!is.null(n_waypoint_eigs)){
        n_waypoint_eigs <- as.integer(n_waypoint_eigs)
    }

    if(!is.null(n_components)){
      n_components <- as.integer(n_components)
    }

    n_highly_variable_genes <- as.integer(n_highly_variable_genes)
    max_n_comps_pca <- as.integer(max_n_comps_pca)
    n_comps_pca <- as.integer(n_comps_pca)
    # Selecting Archetypes
    min_iter <- as.integer(min_iter)
    max_iter <- as.integer(max_iter)
    # Aggregating cells into Archetypes
    seed <- as.integer(seed)
    # Fixed Parameters
    suggested_cells4aa <- as.integer(suggested_cells4aa)
    convergence_epsilon <- as.double(convergence_epsilon)
    minimum_weight <- as.double(minimum_weight)
    bins <- as.integer(bins)
    nth_nearest_neighbor <- as.integer(nth_nearest_neighbor)
    figures_dir <- paste0('./', figures_dir)
    check_dir(figures_dir)
    figures_dir <- paste0(figures_dir, '/')



    message('--- Converting SCE into AnnData ---')
    matrix.dgCMatrix <- counts(sce)
    data.h5ad <- ad$AnnData(t(matrix.dgCMatrix), dtype = 'float32')
    data.h5ad$obs_names <- colnames(matrix.dgCMatrix)
    data.h5ad$var_names <- rownames(matrix.dgCMatrix)
    data.h5ad$obs <- reticulate::r_to_py(as.data.frame(colData(sce)))
    data.h5ad$obs$index <- colnames(sce)
    # enforcing unique feature names
    data.h5ad$var_names_make_unique()

    # file_name <- paste0('./', file_name)
    # sc$write(file_name,   data.h5ad)
    # sc.obj <- sc$read_h5ad(file_name)
    sc.obj <- data.h5ad



    # Counts must be stored in .raw
    # so we need to provide a AnnData object with col and row names.
    h5ad_raw.obj <- sc$AnnData(sc.obj['X'])
    # note: we need to pass a vector such that reticulate will handle it
    h5ad_raw.obj$obs_names <- np$array(sc.obj$obs_names)
    h5ad_raw.obj$var_names <- np$array(sc.obj$var_names)
    # storing the counts with row and col names in the scanpy object
    sc.obj$raw <- h5ad_raw.obj



    message('--- Preprocessing ---')
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

    # setting the number of PCs wanted by subsetting
    # (this can not be performed from terminal)
    sc.obj$obsm['X_pca'] <- sc.obj$obsm['X_pca'][, seq(0, n_comps_pca)]



    cat('--- Computing Kernel Based Archetypes --- \n')
    # The developers of Kernel aa generally suggest 1 prototype each 75 cells
    if(is.null(n_components)){
        n_components <- round(sc.obj$shape[[1]] / suggested_cells4aa)
        n_components <- as.integer(n_components)
    }
    # Setting the number of waypoints
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

    # Performing a Kernel 'Trick'
    kernel_aa.model$construct_kernel_matrix()

    ## This will raise error if run from terminal
    if(plots){
        message('--- Plotting the Kernel Matrix ---')
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
    # this could break if the number of waypoints is higher
    # than the number of archetypes ? probably!


    kernel_aa.model$fit(min_iter = min_iter,
                        max_iter = max_iter)
    # this could break easily!
    # better put a very high threshold!!!

    # note: one could manually force the model to run
    # for more iterations by looping
    # over … model.step()



    message('--- Storing Results ---')
    # General settings to store results
    n_ka <- length(unique(membership.vec))
    patterns_names <- paste0(rep('KA_', n_ka), seq_len(n_ka))
    result_name <- change_default_name(result_name, reducedDimNames(sce))

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
    # assinging cell names to the membership vector
    names(membership.vec) <- rownames(membership.matrix)

    colData(sce)[[result_name]] <- membership.vec



    # Cell View
    kaa_H.dgCMatrix <- Matrix::Matrix(t(kernel_aa.model$A_), sparse = T)
    colnames(kaa_H.dgCMatrix) <- patterns_names
    rownames(kaa_H.dgCMatrix) <- colnames(sce)
    reducedDims(sce)[[result_name]] <- kaa_H.dgCMatrix



    # Gene View
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

    # storing  result
    ## Gene View
    kaa_W.dgRMatrix <- t(kernel_aa_W.obj$layers['log1p'])
    rownames(kaa_W.dgRMatrix) <- rownames(sce)
    colnames(kaa_W.dgRMatrix) <- patterns_names
    metadata(sce)[[result_name]] <- kaa_W.dgRMatrix



    if(plots){
      message('--- Quality metrics Plots ---')
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
      # Lower values of compactness suggest more compact/lower variance archetypes
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
      ## Specifying the cluster variables as 'cell_types' could be interestingly.
      ## It would force to compute the distances between cells belonging to the same
      ## cell type… . But does not seem to work as expected!

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



    if(!return_model) return(sce)
    return(list(obj = sce, model = kernel_aa.model))
}




