#' Non Linear Decomposition Approach : Variational AutoEncoder (VAE)
#'
#' This function performs non linear decomposition leveraging a Variational
#' AutoEncoder. It stores the results directly inside the SingleCellExperiment
#' (SCE) input object. Also It tries to leverage acceleration (GPU, or others…)
#' if it is available. note: on Mac ARM architectures it will fail to use 'gpu'
#' acceleration even if it is available and is detected!
#'
#' @param sce <SingleCellExperiment object> SCE object
#' @param n_components <integer> number of latents (VAE) desired.
#' @param result_name <character> default 'va';
#' Name used to store the result in the SingleCellExperiment object.
#' @param max_epochs <integer> default 400; Number of epochs to fit the VAE.
#' The higher the number of epochs, the better (hopefully) the results.
#' note : it is not always the case given that we are working with Neural
#' Networks. A very naive approach would be to fit the model with a moderate
#' number of epochs and try increasing the number later on…
#' @param n_highly_variable_genes <integer> default 5000; number of Highly
#' Variable Genes to be selected for downstream analysis and fitting the
#' VAE.
#' @param assay <character> default 'counts'; Assay storing the raw counts.
#' @param categorical_covariate_keys < vector of characters or NULL> default
#' NULL; states the fields of categorical covariates to correct for.
#' @param continuous_covariate_keys < vector of characters or NULL> default
#' NULL; states the fields of continuous covariates to correct for.
#' @param target_sum <integer or NULL> default NULL; target sum used to
#' normalize the raw counts. If NULL it uses the mean of all the cells in the
#' raw count matrix. Otherwise is generally set to 1e4.
#' @param accelerator <character> default 'cpu'. Type of accelerations to be
#' used to fit the model. It can be one of :
#' *cpu* : Central Processing Unit (CPU), vanilla use of the CPU's processors.
#' *gpu* : Graphics Processing Unit (GPU), preferable acceleration when
#' note: gpu acceleration, due to inherent variability in GPU processing of data
#' available will introduce slight differences every time is run, thus making
#' the process less reproducible, but much faster!
#' Other types of acceleration are *tpu*, *ipu*, *hpu*, *mps*, *auto*
#' as well as custom accelerator instances.
#' note: Only *cpu* is available for MAC ARCH architectures, for now (2024)
#' @param ignore_warnings <bool> default FALSE; Whether to avoid warnings.
#' note: using verbose = FALSE and ignore warnings = TRUE may result in
#' the suppression of warnings due to suppressing most of the messages from
#' the function.
#' @param envname <character> default 'r-decomp';
#' Specify the name of the python virtual
#' environment to be used. If it does not exists it will create one and use it.
#' @param return_model <bool> default FALSE; Whether to return also
#' the model and not only
#' the SingleCellExperiment object.
#' @param seed <integer> default 42; to set the seed for reproducibility.
#' @param verbose <bool> default FALSE; Whether to be prompted with message
#' for each step of the analysis.
#' @param ... <extra parameters for the function train>
#' @return either a SingleCellExperiment object with both H and W
#' representations or the SingleCellExperiment object and the model
#' used to perform X
#' @examples
#' library(packageX)
#' data(sce)
#' decomp_va(sce, n_components = 5, accelerator = 'cpu', verbose = TRUE, max_epochs = 10)
#' @export
decomp_va <- function(sce,
                      n_components,
                      max_epochs=400,
                      n_highly_variable_genes=5000,
                      assay='counts',
                      categorical_covariate_keys=c(),
                      continuous_covariate_keys=c(),
                      target_sum=NULL,
                      accelerator='cpu',
                      ignore_warnings=TRUE,
                      result_name='va',
                      envname='r-decomp',
                      return_model=FALSE,
                      seed=42,
                      verbose=FALSE,
                      ...){
    ### Description ###
    # Fits a Variational Autoencoder to extract latents|prototypes

    # example usage
    #decomp_va(sce, n_components = 14)

    if(verbose){
        message('--- Checking Packages ---')
    }
    is_package_installed('reticulate') # using python code in R
    is_package_installed('sceasy')     # moving from and to sce and anndata

    packages.vec <- c('scvi-tools',  # variational autoencoder
                      'scanpy',      # single cell analysis framework
                      'scikit-misc', # handy functions
                      'scipy')       # data handlying

    is_python_package_installed(envname = envname,
                                packages.vec = packages.vec)
    # to avoid warnigns
    ignore_warnings(ignore_warnings = ignore_warnings, verbose = verbose)

    # delay_load to pass CRAN testing…
    sc    <- reticulate::import('scanpy', delay_load = TRUE)
    scvi  <- reticulate::import('scvi', delay_load = TRUE)
    os    <- reticulate::import('os', delay_load = TRUE)
    scipy <- reticulate::import('scipy', delay_load = TRUE)

    # enforcing the type to avoid errors in python
    n_components            <- as.integer(n_components)
    n_highly_variable_genes <- as.integer(n_highly_variable_genes)
    max_epochs              <- as.integer(max_epochs)
    seed                    <- as.integer(seed)

    # setting the seed for reproducibility
    scvi$settings$seed <- seed # scvi-tools seed

    # checking for possible accelerations via MPS
    is_MPS_avaliable(envname = envname, verbose = verbose)


    if(verbose){
        message('--- Converting the SCE into AnnData ---')
    }
    # the matrix will be already transposed! cells x features
    adata <- sceasy::convertFormat(sce,
                                   from = "sce",
                                   to = "anndata",
                                   main_layer = "counts",
                                   drop_single_values = FALSE)

    # Converting a CSC sparse matrix to CSR sparse matrix
    # to run faster the model
    if(scipy$sparse$csc$isspmatrix_csc(adata$X)){
        adata$X <- scipy$sparse$csr_matrix(adata$X)
    }



    #message('--- Preprocessing ---') #EXT



    if(verbose){
        message('--- Feature Selection ---')
    }
    # Normalization
    adata$layers[assay] <- adata$X$copy()  # preserve counts
    sc$pp$normalize_total(adata$copy(), target_sum = target_sum)
    sc$pp$log1p(adata$copy())
    adata$raw <- adata  # freeze the state in `.raw`

    # Selection of Highly Variable Genes
    sc$pp$highly_variable_genes(adata,
                                n_top_genes = n_highly_variable_genes,
                                subset = TRUE, # modifying the object keeps only HvGs
                                layer = assay,
                                flavor = "seurat_v3")



    if(verbose){
        message('--- Fitting Variational Autoencoder ---')
    }
    scvi$model$SCVI$setup_anndata(adata,  # make it work with sparse CSR matrix
                                  layer = assay,
                                  categorical_covariate_keys = categorical_covariate_keys,
                                  continuous_covariate_keys = continuous_covariate_keys)

    model <- scvi$model$SCVI(adata, n_latent = n_components)

    # EXT…
    # os$environ['PYTORCH_ENABLE_MPS_FALLBACK'] = '1' # Make the fall back work!

    # fitting the model
    model$train(accelerator = accelerator,
                max_epochs = max_epochs,
                ...)

    va_h.dense <- model$get_latent_representation()
    # note: to sample the distribution we can use give_mean = FALSE



    if(verbose){
        message('--- Storing Results ---')
    }

    sce <- store_H(sce = sce,
                   h.matrix = va_h.dense,
                   result_name = result_name,
                   latent_name = 'VAE')



    # this is the return
    return_model(sce = sce, model = model, return_model = return_model)
}



