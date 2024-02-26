#' Decomposition Approach : Knowledge Based Decomposition
#'
#' This function performs Knowledge Based Decomposition leveraging Spectra
#' implementation. It stores the results directly inside the
#' SingleCellExperiment object.
#'
#' @param sce <SingleCellExperiment object> SCE object
#' @param metadata_key <'character or NULL'> default NULL; Specified prior
#' knowledge to be used, specifically groups of cells (like clusters, types,
#' states, …). If is set to NULL, prior_knowledge should not have any metadata
#' specific gene set!
#' @param prior_knowledge <list of lists of named gene sets> default NULL;
#' List of lists, where each sub list is a named gene set to be tested.
#' Gene sets can be metadata specific or global. If global the gene set will
#' be tested on all the cells, regardless of the grouping due to the metadata.
#' The gene sets should be grouped in a named list; the name of the list must
#' match a name of the metadata_key.
#' look at the example to get an idea of the structure…
#' @param n_components <integer or NULL> default NULL. Generally should match
#' the number of gene sets (authors suggestion).
#' @param assay <character> default 'logcounts'; assay used for the analysis.
#' @param lam <double> default 0.1; trust in prior knowledge.
#' @param n_hvgs <integer or NULL> default 5000; Number of Highly Variable Genes
#' to be selected and used for downstream analysis. If set to 0 or NULL, it will
#' use all the available genes.
#' @param n_markers <integer> default 50; number of genes for each gene set.
#' It accounts also for imputed genes during the analysis based on the data.
#' @param num_epochs <integer> default 10000; number of epocs used to fit
#' the model. 10.000 is the recomended default.
#' @param min_gs_num <integer> default 3; minimum number of gene in a gene set.
#' @param use_metadata <bool> default TRUE; Whether to use the metadata_key
#' or not. If set to FALSE, the user should not specify metadata specific gene
#' sets.
#' @param use_weights <bool> default TRUE; edge weights are estimated based
#' on graph structure and used throughout training (from Spectra).
#' @param delta <double> default 0.001; Scaling parameter to make the values
#' comparable. It varies depending on data and gene sets.
#' It can range from 0.5 up to 0.001.
#' @param gene_set_confidence_assignment <double> default 0.02;
#' minimum fraction of overlapping genes to assign the gene set to a
#' given metadata group of cells.
#' @param additional_results <bool> default TRUE; Whether to store in the
#' metadata field of the SingleCellExperiment object also Markers of each Set,
#' Interpretation of Factors, Influence of the Prior on the Factors and the
#' Gene-Gene graph.
#' @param ignore_warnings <bool> default TRUE; Whether to silence warnings.
#' @param ggg_metadata <character or NULL>; Specifying a metadata
#' group (like a given cell type) to obtain a Metadata-Specific Gene-Gene Graph,
#' thus the name ggg_metadata.
#' @param result_name <character> default 'kb';
#' Name of used to store the result in the SingleCellExperiment object.
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
#' # prior_knowledge = list('0' = list('gene_set_1' = list("KLHL17",
#' #                                                      "CCNL2",
#' #                                                      "ATAD3B",
#' #                                                      "NOL9"),
#' #                                  'gene_set_2' = list("CENPS-CORT",
#' #                                                      "FBXO2",
#' #                                                      "KLHL21")),
#' #                        '1' = list('gene_set_1' = list("NPPA",
#' #                                                       "FHAD1",
#' #                                                       "FBXO42",
#' #                                                       "ATP13A2"),
#' #                                   'gene_set_3' = list("PADI2",
#' #                                                       "AHDC1",
#' #                                                       "PPP1R8",
#' #                                                       "LAPTM5")),
#' #                  # set of genes which will be evaluated for all the
#' #                  # cells in data…
#' #                  'global' = list('gene_set_4' = list("AZIN2",
#' #                                                      "ZSCAN20",
#' #                                                      "ZC3H12A")))
#' # decomp_kb(sce,
#' #           metadata_key = 'is.doublet',
#' #           prior_knowledge = prior_knowledge,
#' #           n_hvgs = NULL,
#' #           num_epochs = 100)
#' @export
decomp_kb <- function(sce,
                      metadata_key,
                      prior_knowledge,
                      n_components=NULL,
                      assay='logcounts',
                      lam=0.10,
                      n_hvgs=5000,
                      n_markers=50,
                      num_epochs=10000,
                      min_gs_num=3,
                      use_metadata=TRUE,
                      use_weights=TRUE,
                      delta=0.001,
                      gene_set_confidence_assignment=0.2,
                      additional_results=TRUE,
                      ignore_warnings=TRUE,
                      ggg_metadata=NULL,
                      result_name='kb',
                      envname='r-decomp',
                      return_model=FALSE,
                      seed=42,
                      verbose=FALSE
){

    ### Description ###
    # Performs Knowledge Based Decomposition via Spectra

    # example usage

    #  we need to define some gene sets, both global and cell specific…
    #  Gene sets must be defined via nested lists…
    # '0' and '1' are the labels assigned to cell types,
    #  in this case 0 is for a high quality cell,
    #  '1' for a doublet!
    #  The name of the nested lists must match with the names assigned
    #  to each cell in the metadata of the sce.

    # #example… with the metadata field 'is.doublet'
    # metadata_key = 'is.doublet'
    # prior_knowledge = list('0' = list('gene_set_1' = list("KLHL17",
    #                                                       "CCNL2",
    #                                                       "ATAD3B",
    #                                                       "NOL9"),
    #                                   'gene_set_2' = list( "CENPS-CORT",
    #                                                        "FBXO2",
    #                                                        "KLHL21")),
    #                        '1' = list('gene_set_1' = list("NPPA",
    #                                                       "FHAD1",
    #                                                       "FBXO42",
    #                                                       "ATP13A2"),
    #                                   'gene_set_3' = list("PADI2",
    #                                                       "AHDC1",
    #                                                       "PPP1R8",
    #                                                       "LAPTM5")),
    #                        # set of genes which will be evauated for all the
    #                        # cells in data…
    #                        'global' = list('gene_set_4' = list("AZIN2",
    #                                                            "ZSCAN20",
    #                                                            "ZC3H12A")))

    ### Running Knowledge Based decomposition…
    ## fixing the number of patterns|prototypes to find
    # decomp_kb(sce,
    #           n_components = 14,
    #           metadata_key = NULL,
    #           prior_knowledge = prior_knowledge$global,
    #           use_metadata = FALSE)



    if(verbose){
        message('--- Checking Packages ---')
    }
    is_package_installed('reticulate')
    is_package_installed('scran')      # selecting HvGs
    is_package_installed('sceasy')     # moving from and to sce and anndata
    packages.vec <- c('scSpectra', 'pandas', 'pyvis', 'numpy')
    is_python_package_installed(envname      = envname,
                                packages.vec = packages.vec)

    spectra <- reticulate::import('Spectra', delay_load = TRUE)
    pd      <- reticulate::import('pandas', delay_load = TRUE)
    np      <- reticulate::import('numpy', delay_load = TRUE)

    ignore_warnings(ignore_warnings = ignore_warnings,
                    envname = envname,
                    verbose = verbose)

    # enforcing the type on the variables to avoid problems in python
    if(!is.null(n_components)){
        n_components <- as.integer(n_components)
    }
    n_markers <- as.integer(n_markers)
    num_epochs <- as.integer(num_epochs)
    if(min_gs_num == 0){
        min_gs_num <- 1
        warning('Forcing the min_gs_num to 1; must be a positive integer')
    }
    min_gs_num <-as.integer(min_gs_num)

    metadata.in_use <- !is.null(metadata_key)
    if(metadata.in_use){
      SummarizedExperiment::colData(sce)[[metadata_key]] <- as.character(SummarizedExperiment::colData(sce)[[metadata_key]])
    }

    # converting SCE into AnnData (with sparse matrix)
    adata <- sce2adata_sparse(sce = sce, envname = envname, main_layer = assay)



    ###### this could be a function on its own #################################
    # to match correctly the metadata_key with prior knowledge we need to
    # transform into characters the chosen filed for the metadata_key
    # S is for string
    if(metadata.in_use){
        metadata.dtype_NOT_string <- as.logical(toupper(adata$obs[[metadata_key]]$dtype != 'S'))

        # enforcing the type of the medata obs
        if(metadata.dtype_NOT_string){
            adata$obs[[metadata_key]] <- adata$obs[[metadata_key]]$astype('int')
            adata$obs[[metadata_key]] <- adata$obs[[metadata_key]]$astype('string')
        }

        # checking if prior_knowledge is found in the data
        check_gene_set_dictionary <- spectra$Spectra_util$check_gene_set_dictionary
        # to convert from r to python correctly the variable type.
        prior_knowledge <- reticulate::r_to_py(prior_knowledge)
        prior_knowledge <- check_gene_set_dictionary(adata       = adata,
                                                    annotations = prior_knowledge,
                                                    obs_key     = metadata_key,
                                                    global_key  = 'global')
    }
    ############################################################################



    ############################################################################
    # Finding HVGs if we want to use them
    # this can be improved! but should work fine for now…
    if(!is.null(n_hvgs)){
        if(n_hvgs != 0){
           hvgs <- scran::getTopHVGs(sce, n = n_hvgs)
           hvgs.bool <- rownames(sce) %in% hvgs
           adata$var['highly_variable'] <- hvgs.bool
           # and setting the flag to use HVGs
           use_highly_variable <- TRUE
        }else{
            use_highly_variable <- FALSE
        }
    }else{
        use_highly_variable <- FALSE
    }
    ############################################################################



    ############################################################################
    # fitting Spectra model; Performing Knowledge Based Decomposition
    if(verbose){
        message('--- Running Knowledge Based Decomposition ---')
    }
    model <- spectra$est_spectra(adata               = adata,
                                  L                   = n_components,
                                  gene_set_dictionary = prior_knowledge,
                                  use_highly_variable = use_highly_variable,
                                  use_cell_types      = use_metadata,
                                  cell_type_key       = metadata_key,
                                  use_weights         = use_weights,
                                  lam                 = lam,
                                  delta               = delta,
                                  n_top_vals          = n_markers,
                                  overlap_threshold   = gene_set_confidence_assignment,
                                  min_gs_num          = min_gs_num,
                                  num_epochs          = num_epochs)
    ############################################################################



    ############################################################################
    if(verbose){
        message('--- Storing Results ---')
    }
    # general parameters
    result_name <- change_default_name(result_name, SingleCellExperiment::reducedDimNames(sce))
    # pattern names will be defined later…

    # Cell View
    kb_h.dense <- model$cell_scores
    kb_h.dgCmatrix <- methods::as(kb_h.dense, 'sparseMatrix')
    # Setting the number of patterns if they were selected by default
    if(is.null(n_components)){
      n_components <- ncol(kb_h.dgCmatrix)
    }
    patterns_names <- paste0(rep('KB_', n_components), seq_len(n_components))
    colnames(kb_h.dgCmatrix) <- patterns_names #since we want to keep the model
    # assignment of cells to specific-patterns
    rownames(kb_h.dgCmatrix) <- colnames(sce)
    SingleCellExperiment::reducedDims(sce)[[result_name]] <- kb_h.dgCmatrix

    # Gene View
    # Spectra factors the patterns of the W matrix
    kb_w.dense <- Matrix::t(as.matrix(adata$uns['SPECTRA_factors']))
    kb_w.dgCmatrix <- methods::as(kb_w.dense, 'sparseMatrix')
    colnames(kb_w.dgCmatrix) <- patterns_names #since we might be interested
    #in the assignment provided by the model of a gene set to a specific gene.
    if(use_highly_variable){
        imputed_genes <- as.vector(np$asarray(adata$var['spectra_vocab']))
        rownames(kb_w.dgCmatrix) <- rownames(sce)[imputed_genes]
    }else{
        rownames(kb_w.dgCmatrix) <- rownames(sce)
    }
    S4Vectors::metadata(sce)[[result_name]] <- kb_w.dgCmatrix
    ############################################################################



    ############################################################################
    # Additional Results
    if(additional_results){
        ## Markers of each Set
        if(use_highly_variable){
            marker_order <- paste0(rep('top_', ), seq_len(sum(imputed_genes)), rep('_marker'))
        }else{
            marker_order <- paste0(rep('top_', ), seq_len(n_markers), rep('_marker'))
        }
        markers4set <- Matrix::t(as.matrix(adata$uns['SPECTRA_markers']))
        colnames(markers4set) <- patterns_names
        rownames(markers4set) <- marker_order
        S4Vectors::metadata(sce)[[paste0(result_name, '_markers')]] <- markers4set

        ## Interpretation of Factors
        # 'index' + '-X-' + 'cell type specificity' + '-X-' + 'assigned label', ...]
        factors_interpretation <- adata$uns['SPECTRA_overlap']
        fact_interpretation <- paste0(result_name, '_factors_interpretation')
        S4Vectors::metadata(sce)[[fact_interpretation]] <- factors_interpretation

        ## Influence of the Prior on the Factors
        prior_influence <- model$return_eta_diag()
        names(prior_influence) <- patterns_names
        S4Vectors::metadata(sce)[[paste0(result_name, '_prior_influence')]] <- prior_influence

        ## Gene-Gene graph
        if(!is.null(ggg_metadata)){ # do we want a cell-type specific graph?
          soft_graph <- as.matrix(model$return_graph(ct = ggg_metadata))
        } else {
          soft_graph <- as.matrix(model$return_graph())
        }
        # storing the graph
        S4Vectors::metadata(sce)[[paste0(result_name, '_GeneGeneGraph')]] <- soft_graph
    }
    ############################################################################



    # this is the return
    return_model(sce = sce, model = model, return_model = return_model)
}


# More notes on this function…
# you can run it without cell type specific patterns
# by using a simple dict and turning NULL the flag of cell_type_key

# nested list of lists gene_set_dict | cell type_i | gene_set_i
# gene sets are not cell specific
# Note that one key in the dictionary must be 'global' with the corresponding
# value being a dictionary of gene sets which apply to all cells

# The cell type labels have to match with the cell type
# labels in the gene set dictionary

# if use_metadata == False then maps gene set names to gene sets ;
# must contain "global" key in addition to every unique cell type under .obs.
# Note that if you have a cell type in your adata for which
# you do not have any gene sets
# in your gene set annotation dictionary you must include an
# empty dictionary under that cell type key.

# Consider looking at the Colab notebook
# to leverage a representation of the network

