decomp_kb <- function(sce,
                      metadata_key,
                      prior_knowledge,
                      n_components=NULL, # number of gene sets + 1 (thumb rule)
                      assay='logcounts',
                      # or a dictionary that maps cell type
                      # to an integer per cell type
                      lam=0.10, # trust in the prior (the higher, the higher the trust in it!)
                      use_highly_variable=TRUE,
                      n_markers=50,
                      num_epochs=1000,
                      min_gs_num=3,
                      use_cell_types=TRUE,
                      use_weights=TRUE, # edge weights are estimated based on graph
                      # structure and used throughout training
                      delta=0.001, # scaling to make the values comparable
                      # varies depending on data and gene sets, try between 0.5 and 0.001
                      gene_set_confidence_assignment = 0.2, # fraction of overlapping genes
                      result_name='kb',
                      additional_results=TRUE,
                      return_model=FALSE,
                      envname='r-decomp',
                      ggg_cell_type=NULL,
                      ignore_warnings=TRUE){

    ### Description ###
    # Performs Knowledge Based Decomposition via Spectra

    # example usage

    # we need to define some gene sets, both global and cell specific…
    # Gene sets must be defined via nested lists…
    # '0' and '1' are the labels assigned to cell types,
    #  in this case 0 is for a high quality cell,
    #  '1' for a doublet!
    #  The name of the nested lists must match with the names assigned
    #  to each cell in the metadata of the sce.

    # example… with the metadata field 'is.doublet'
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
    #           n_components = 14
    #           metadata_key = NULL,
    #           prior_knowledge = prior_knowledge$global,
    #           use_cell_types = FALSE)

    ## Picking the same number of gene set used
    # decomp_kb(sce,
    #           metadata_key = 'is.doublet',
    #           prior_knowledge = prior_knowledge)


    message('--- Checking Packages ---')
    is_package_installed('reticulate')
    is_package_installed('sceasy')     # moving from and to sce and anndata

    packages.vec = c('scSpectra', 'pandas', 'pyvis')

    is_python_package_installed(envname      = envname,
                                packages.vec = packages.vec)

    spectra <- reticulate::import('Spectra')
    pd      <- reticulate::import('pandas')

    if(ignore_warnings){
        warnings <- reticulate::import('warnings')
        warnings$filterwarnings("ignore")
    }

    # enforcing integer type on variables
    if(!is.null(n_components)){
        n_components <- as.integer(n_components)
    }
    n_markers = as.integer(n_markers)
    num_epochs = as.integer(num_epochs)

    if(min_gs_num == 0){
        min_gs_num <- 1
        warning('> Forcing the min_gs_num to 1; must be a positive integer <')
    }
    min_gs_num = as.integer(min_gs_num)

    metadata.in_use <- !is.null(metadata_key)
    if(metadata.in_use){
        colData(sce)[[metadata_key]] <- as.character(colData(sce)[[metadata_key]])
    }


    # converting SCE into AnnData (with sparse matrix)
    adata <- sce2adata_sparse(sce = sce, envname = envname, main_layer = assay)

    ###### this could be a function on its own ######

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
        prior_knowledge <- reticulate::r_to_py(prior_knowledge)
        print('ok')
        prior_knowledge = check_gene_set_dictionary(adata       = adata,
                                                    annotations = prior_knowledge,
                                                    obs_key     = metadata_key,
                                                    global_key  = 'global')
    }

    ############################################################################


    message('--- Running Knowledge Based Decomposition ---')
    model = spectra$est_spectra(adata               = adata,
                                L                   = n_components,
                                gene_set_dictionary = prior_knowledge,
                                use_highly_variable = use_highly_variable,
                                use_cell_types      = use_cell_types,
                                cell_type_key       = metadata_key,
                                use_weights         = use_weights,
                                lam                 = lam,
                                delta               = delta,
                                n_top_vals          = n_markers,
                                overlap_threshold   = gene_set_confidence_assignment,
                                min_gs_num          = min_gs_num,
                                num_epochs          = num_epochs)



    message('--- Storing Results ---')
    # General Parameters
    result_name <- change_default_name(result_name, reducedDimNames(sce))
    # pattern names will be defined later…

    # Cell View
    kb_h.dense <- model$cell_scores
    kb_h.dgCmatrix <- as(kb_h.dense, 'sparseMatrix')
    # Setting the number of patterns if they were selected by default
    if(is.null(n_components)){
      n_components <- ncol(kb_h.dgCmatrix)
    }
    patterns_names <- paste0(rep('KB_', n_components), seq_len(n_components))
    colnames(kb_h.dgCmatrix) <- patterns_names
    rownames(kb_h.dgCmatrix) <- colnames(sce)
    reducedDims(sce)[[result_name]] <- kb_h.dgCmatrix

    # Gene View
    # Spectra factors the patterns of the W matrix
    kb_w.dense <- t(as.matrix(adata$uns['SPECTRA_factors']))
    kb_w.dgCmatrix <- as(kb_w.dense, 'sparseMatrix')
    colnames(kb_w.dgCmatrix) <- patterns_names
    rownames(kb_w.dgCmatrix) <- rownames(sce)
    metadata(sce)[[result_name]] <- kb_w.dgCmatrix

    # Additional Results
    if(additional_results){
        ## Markers of each Set
        marker_order <- paste0(rep('top_', ), seq_len(n_markers), rep('_marker'))
        markers4set <- t(as.matrix(adata$uns['SPECTRA_markers']))
        colnames(markers4set) <- patterns_names
        rownames(markers4set) <- marker_order
        metadata(sce)[[paste0(result_name, '_markers')]] <- markers4set

        ## Interpretation of Factors
        # 'index' + '-X-' + 'cell type specificity' + '-X-' + 'assigned label', ...]
        factors_interpretation <- adata$uns['SPECTRA_overlap']
        fact_interpretation <- paste0(result_name, '_factors_interpretation')
        metadata(sce)[[fact_interpretation]] <- factors_interpretation

        ## Influence of the Prior on the Factors
        prior_influence <- model$return_eta_diag()
        names(prior_influence) <- patterns_names
        metadata(sce)[[paste0(result_name, '_prior_influence')]] <- prior_influence

        ## Gene-Gene graph
        if(!is.null(ggg_cell_type)){ # do we want a cell-type specific graph?
          soft_graph <- as.matrix(model$return_graph(ct = ggg_cell_type))
        } else {
          soft_graph <- as.matrix(model$return_graph())
        }
        # storing the graph
        metadata(sce)[[paste0(result_name, '_GeneGeneGraph')]] <- soft_graph
    }

if(!return_model) return(sce)
return(list(obj = sce, model = model, adata = adata))
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


# if use_cell_types == False then maps gene set names to gene sets ;
# must contain "global" key in addition to every unique cell type under .obs.
# Note that if you have a cell type in your adata for which
# you do not have any gene sets
# in your gene set annotation dictionary you must include an
# empty dictionary under that cell type key.

# Consider looking at the Colab notebook
# to leverage a representation of the network





