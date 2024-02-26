#' Partition Approach : Graph Clustering
#'
#' This function performs Graph Clustering via either Louvain or Leiden on a
#' SingleCellExperiment (SCE) object. It stores the memebership vector,
#' the result of the function, directly in the colData slot of the SCE Object.
#'
#' @param sce <SingleCellExperiment object> SCE object
#' @param resolution <double> default 1; resolution used for clustering.
#' The higher it is, the finer the clusters; The lower it is, the more coarser
#' the clusters.
#' @param method <character> default *leiden*. Method used for clustering, can
#' be one of:
#' *louvain* : modularity based clustering;
#' *leiden*  : modularity based clustering with significant improvement over
#' louvain.
#' @param n_pcs <integer> default 50; number of PCs to be used after Principal
#' Component Analysis as reduction method.
#' @param n_neighbors <integer> default 10; number of Nearest Neighbors to be
#' used for the K-Nearest Neighbors (KNN) Graph.
#' @param svd_solver <character> default *auto*; Type of solver to perform
#' Principal Component Analysis. Can be one of :
#' *arpack*    : default solver form the ARPACK wrapper in SciPy;
#' *randomized*: for a randomized PCA (Halko - 2009);
#' *auto*      : the solver is picked depending on the type of data used;
#' *lobpcg*    : experimental solver to process quickly a sparse matrix.
#' note : *lobpcg* currently (Q1-2024) only works with the ‘arpack’
#' or ‘lobpcg’ solvers.
#' @param ignore_warnings <bool> Whether to ignore warning or not.
#' @param result_name <character> default 'leiden';
#' Name used to store the result in the SingleCellExperiment object.
#' @param envname <character> default 'r-partition';
#' Specify the name of the python virtual
#' environment to be used. If it does not exists it will create one and use it.
#' @param return_model <bool> default FALSE; Whether to return also
#' the model and not only
#' the SingleCellExperiment object.
#' @param seed <integer> default 42; to set the seed for reproducibility.
#' @param verbose <bool> default FALSE; Whether to be prompted with message
#' for each step of the analysis.
#' @return either a SingleCellExperiment object with both H and W
#' representations or the SingleCellExperiment object and the model
#' used to perform X
#' @examples
#' #partition_graph(sce)
#' @export
partition_graph <- function(sce,
                            resolution=1,
                            method='leiden',
                            n_pcs=50,
                            n_neighbors=10,
                            svd_solver='auto',
                            ignore_warnings=TRUE,
                            result_name='leiden',
                            envname='r-partition',
                            return_model=FALSE,
                            seed=42,
                            verbose=FALSE){



    ### Description ###
    # Performs Leiden or Louvain Clustering via scanpy



    if(verbose){
        message('--- Checking packages ---')
    }
    is_package_installed('reticulate')
    is_python_package_installed(envname = envname,
                                packages.vec = c('scanpy',
                                                 'anndata',
                                                 'leidenalg',
                                                 'louvain',
                                                 'numpy'))

    ignore_warnings(ignore_warnings = ignore_warnings,
                    envname = envname,
                    verbose = verbose)

    # importing 'only' the clustering approaches of sklrn
    sc       <- reticulate::import('scanpy', delay_load = TRUE)
    anndata  <- reticulate::import('anndata', delay_load = TRUE)
    np       <- reticulate::import('numpy', delay_load = TRUE)
    # enforcing the type of variables to avoid crashes in python
    resolution  <- as.integer(resolution)
    n_pcs       <- as.integer(n_pcs)
    n_neighbors <- as.integer(n_neighbors)
    seed        <- as.integer(seed)

    # from SCE to AnnData
    data.h5ad <- sce2adata_sparse(sce = sce,
                                  envname = envname,
                                  main_layer = 'counts')



    if(verbose){
        message('--- Preprocessing ---')
    }
    # Features normalization and selection
    sc$pp$normalize_total(data.h5ad$copy())
    sc$pp$log1p(data.h5ad$copy())
    sc$pp$highly_variable_genes(adata = data.h5ad)
    sc$pp$scale(data.h5ad)
    # Dimensionality reduction
    sc$tl$pca(data.h5ad, n_comps = n_pcs, svd_solver = svd_solver)
    # KNN graph
    sc$pp$neighbors(data.h5ad, n_neighbors = n_neighbors, n_pcs = n_pcs)



    if(verbose){
        message(paste0('--- Clustering with ', method, ' ---'))
    }
    switch(method,
           'louvain' = sc$tl$louvain(data.h5ad,
                                     resolution = resolution,
                                     random_state = seed),
           'leiden'  = sc$tl$leiden(data.h5ad,
                                    resolution = resolution,
                                    random_state = seed))



    if(verbose){
        message('--- Storing results ---')
    }
    # EXT - this could be a function on its own!
    result_name <- change_default_name(result_name, colnames(SummarizedExperiment::colData(sce)))
    membership.vec <- np$array(data.h5ad[['obs']][method])
    names(membership.vec) <- colnames(sce)
    SummarizedExperiment::colData(sce) <- cbind(SummarizedExperiment::colData(sce), membership.vec)
    colnames(SummarizedExperiment::colData(sce))[colnames(SummarizedExperiment::colData(sce)) == 'membership.vec'] <- result_name



     # this is the return
    return_model(sce = sce, model = data.h5ad, return_model = return_model)
}
