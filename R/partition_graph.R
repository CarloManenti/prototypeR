partition_graph <- function(sce,
                            method='louvain',
                            resolution=0.8,
                            result_name='leiden', # or louvain
                            envname='r-scanpy',
                            n_pcs=50,
                            n_neighbors=10,
                            svd_solver='arpack',
                            seed=42,
                            return_model=FALSE){



    ### Description ###
    # Performs Leiden or Louvain Clustering via scanpy


    # example usage
    # partition_graph(sce)


    message('--- Checking packages ---')
    is_package_installed('reticulate')
    is_python_package_installed(envname = envname,
                                packages.vec = c('scanpy',
                                                 'anndata',
                                                 'leidenalg',
                                                 'louvain'))
    reticulate::use_virtualenv(virtualenv = envname)

    # importing 'only' the clustering approaches of sklrn
    sc       <- reticulate::import('scanpy')
    anndata  <- reticulate::import('anndata')
    #warnings <- reticulate::import('warnings')

    # to avoid future warnings of pandas
    # warnings$simplefilter(action='ignore', category='FutureWarning')

    resolution  <- as.integer(resolution)
    n_pcs       <- as.integer(n_pcs)
    n_neighbors <- as.integer(n_neighbors)
    seed        <- as.integer(seed)

    counts.dgCMatrix <- counts(sce)
    data.h5ad <- sc$AnnData(X = t(counts.dgCMatrix))$copy()
    data.h5ad$obs_names <- colnames(counts.dgCMatrix)
    data.h5ad$var_names <- rownames(counts.dgCMatrix)


    message('--- Preprocessing ---')

    # Features normalization and selection
    sc$pp$normalize_total(data.h5ad)
    sc$pp$log1p(data.h5ad)
    sc$pp$highly_variable_genes(adata = data.h5ad)
    sc$pp$scale(data.h5ad)
    # Dimensionality reduction
    sc$tl$pca(data.h5ad, n_comps = n_pcs, svd_solver = svd_solver)
    # KNN graph
    sc$pp$neighbors(data.h5ad, n_neighbors = n_neighbors, n_pcs = n_pcs)



    message(paste0('--- Clustering with ', method, ' ---'))
    switch(method,
           'louvain' = sc$tl$louvain(data.h5ad,
                                     resolution = resolution,
                                     random_state = seed),
           'leiden'  = sc$tl$leiden(data.h5ad,
                                    resolution = resolution,
                                    random_state = seed))



    message('--- Storing results ---')
    result_name <- change_default_name(result_name, colnames(colData(sce)))
    membership.vec <- data.h5ad[['obs']][method]
    names(membership.vec) <- result_name
    colData(sce) <- cbind(colData(sce), membership.vec)

    if(!return_model) return(sce)
    return(list(obj = sce, model = data.h5ad))
}








