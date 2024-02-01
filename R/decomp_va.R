decomp_va <- function(sce,
                      n_components,
                      max_epochs=400,
                      n_highly_variable_genes=5000, # nHVGs
                      layer='counts',
                      categorical_covariate_keys = c(), # to account for possible categorical batches
                      continuous_covariate_keys = c(),  # to account for possible continous batches
                      target_sum=1e4, # for normalizing (could be the meadian of the dataset or else)
                      accelerator='cpu', # use GPU for TRAINING (must be installed and avaliable!)
                      seed=42,
                      envname='r-decomp',
                      result_name='va',
                      return_model=FALSE,
                      ...){
    ### Description ###
    # Fits a Variational Autoencoder to extract latents|prototypes

    # example usage
    #decomp_va(sce, n_components = 14)

    message('--- Checking Packages ---')
    is_package_installed('reticulate') # using python code in R
    is_package_installed('sceasy')     # moving from and to sce and anndata

    packages.vec = c('scvi-tools', # variational autoencoder
                     'scanpy',     # single cell analysis framework
                     'scikit-misc',
                     'scipy')     # data handlying

    is_python_package_installed(envname = envname,
                                packages.vec = packages.vec)

    reticulate::use_virtualenv(envname)

    sc <- reticulate::import('scanpy')
    scvi <- reticulate::import('scvi')
    os <- reticulate::import('os')
    scipy <- reticulate::import('scipy')

    # enforcing the type
    n_components <- as.integer(n_components)
    n_highly_variable_genes <- as.integer(n_highly_variable_genes)
    max_epochs <- as.integer(max_epochs)
    seed <- as.integer(seed)

    # setting the seed for reproducibility
    scvi$settings$seed = seed # scvi-tools seed

    #is_MPS_avaliable(envname = envname)



    message('--- Converting the SCE into AnnData ---')
    # the matrix will be already transposed! cells x features
    adata <- sceasy::convertFormat(sce,
                                   from = "sce",
                                   to = "anndata",
                                   main_layer = "counts",
                                   drop_single_values = FALSE)

    # Converting a CSC sparse matrix to CSR sparse matrix
    # to run fuster the model
    if(scipy$sparse$csc$isspmatrix_csc(adata$X)){
        adata$X <- scipy$sparse$csr_matrix(adata$X)
    }



    #message('--- Preprocessing ---') #TBD



    message('--- Feature Selection ---')
    # Normalization
    adata$layers[layer] = adata$X$copy()  # preserve counts
    sc$pp$normalize_total(adata$copy(), target_sum = target_sum)
    sc$pp$log1p(adata$copy())
    adata$raw = adata  # freeze the state in `.raw`

    # Selction of Higly Variable Genes
    sc$pp$highly_variable_genes(adata,
                                n_top_genes = n_highly_variable_genes,
                                subset = TRUE, # modifying the object keeping only HVGs
                                layer = layer,
                                flavor = "seurat_v3")



    message('--- Fitting Variational Autoencoder ---')
    scvi$model$SCVI$setup_anndata(adata,  # make it work with sparse CSR matrix
                                  layer = layer,
                                  categorical_covariate_keys = categorical_covariate_keys,
                                  continuous_covariate_keys = continuous_covariate_keys)

    model = scvi$model$SCVI(adata, n_latent = n_components, ...)

    # os$environ['PYTORCH_ENABLE_MPS_FALLBACK'] = '1' # Make the fall back work!

    model$train(accelerator = accelerator,
                max_epochs = max_epochs,
                ...)

    va_h.dense = model$get_latent_representation()
    # note: to sample the distribution we can use give_mean = FALSE



    message('--- Storing Results ---')
    # General info on the results
    patterns_names <- paste0(rep('VA_', n_components), seq_len(n_components))
    result_name <- change_default_name(result_name, reducedDimNames(sce))

    # Cell View
    va_h.dgCmatrix <- as(va_h.dense, 'sparseMatrix')
    colnames(va_h.dgCmatrix) <- patterns_names
    rownames(va_h.dgCmatrix) <- colnames(sce)
    reducedDims(sce)[[result_name]] <- va_h.dgCmatrix

    if(!return_model) return(sce)
    return(list(obj = sce, model = model))
}



